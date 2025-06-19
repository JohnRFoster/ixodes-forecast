## main hindcast workflow script
## 1 - setup
## 3 - state data intake
## 4 - weather data intake
## 5 - create informative priors
## 6 - get initial conditions
## 7 - forecast step
## 7a - mice
## 7b - ticks
## 8 - analysis step
## 9 - save

# =========================================== #
#       1 - setup
# =========================================== #

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(lubridate)
library(nimble)
library(parallel)

tests <- expand_grid(
  use_nmme = TRUE,
  use_mice = TRUE,
  remove = c(FALSE, TRUE)
)

test <- 2

tests |>
  slice(test)

source("R/data_assimilation.R")
source("R/function_run_hindcast_nimble.R")
source("R/function_scale_met_forecast.R")
source("R/functions_mice.R")
source("R/functions_misc.R")

n_slots <- 3 # 3 chains
# use_nmme <- TRUE
use_nmme <- tests |> slice(test) |> pull(use_nmme)
production <- TRUE
n_iter <- 20000
nmc <- 2000
horizon <- 365
# use_mice <- TRUE
use_mice <- tests |> slice(test) |> pull(use_mice)

dir_data <- "data"
dir_nmme <- file.path(dir_data, "NMME")
dir_save <- "out"

# field sites
sites <- c("Green", "Henry", "Tea")

# all possible combinajobstions for hindcast experiments
all_jobs <- expand.grid(
  paramsFrom = sites,
  ticksFrom = sites,
  remove = c(TRUE, FALSE)
)

job_num <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# for testing locally
if (is.na(job_num)) {
  job_num <- 9
}

# remove <- all_jobs$remove[job_num]
remove <- tests |> slice(test) |> pull(remove)
params_job <- all_jobs$paramsFrom[job_num]
ticks_job <- all_jobs$ticksFrom[job_num]
mice_job <- ticks_job
model_job <- if_else(use_mice, "WithWeatherAndMiceGlobal", "Weather")
ua_from <- if_else(
  use_mice,
  "mice_ic_parameter_process",
  "ic_parameter_process"
)
calibration_mice <- paste(mice_job, "Control")
calibration_ticks <- paste(ticks_job, "Control")
experiment <- "base"
if (use_mice) {
  experiment <- paste(experiment, "mna", sep = "_")
}
if (remove) {
  experiment <- paste(experiment, "remove", sep = "_")
}
if (use_nmme) {
  experiment <- paste(experiment, "nmme", sep = "_")
}

# =========================================== #
#       1 - state data intake
# =========================================== #

df_tick <- read_csv(file.path(dir_data, "Ticks2006_2021.csv"))
df_mice <- read_csv(file.path(dir_data, "Mice2006_2021.csv")) # currently tick data only through 2020

data_tick <- df_tick |> filter(Grid == calibration_ticks)
data_mice <- df_mice |> filter(Grid == calibration_mice)

dates_tick <- data_tick$Date
dates_mice <- data_mice$Date

start_date <- if_else(use_mice, min(dates_mice), min(dates_tick) - 31) # start one month before first tick observation
end_date <- max(dates_tick) # end at last tick observation
hindcast_period <- seq.Date(start_date, end_date, by = 1)

# hindcast_period represents every single day in the time series
# for a majority of the hindcasts, we only need to run the simulation
# after each observation because the weather is treated as known, then
# also only need to run the hindcasts during the sampling season (not winter)

# =========================================== #
#       1 - weather data intake
# =========================================== #

# we will need the means/SDs from the calibration period
hist_means <- scale_met_forecast()

# get the observed data from Cary and center to historical mean
cary_met <- read_csv(file.path(dir_data, "Cary_Met_Data_Daily.csv")) |>
  suppressWarnings()
met_pre_2021 <- cary_met |>
  slice(-c((nrow(cary_met) - 364):nrow(cary_met))) |>
  mutate(Date = mdy(DATE)) |>
  filter(Date >= "1995-01-01")
met_2021 <- cary_met |>
  slice((nrow(cary_met) - 364):nrow(cary_met)) |>
  mutate(Date = ymd(DATE))

weather_corrected <- bind_rows(met_pre_2021, met_2021) |>
  select(Date, MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) |>
  arrange(Date) |>
  mutate(year = year(Date)) |>
  group_by(year) |>
  mutate(
    gd = if_else(MAX_TEMP > 10, MAX_TEMP - 10, 0),
    cgdd = cumsum(gd),
    MAX_TEMP = (MAX_TEMP - hist_means$means["MAX_TEMP"]) /
      hist_means$sds["MAX_TEMP"],
    MAX_RH = (MAX_RH - hist_means$means["MAX_RH"]) / hist_means$sds["MAX_RH"],
    MIN_RH = (MIN_RH - hist_means$means["MIN_RH"]) / hist_means$sds["MIN_RH"],
    TOT_PREC = (TOT_PREC - hist_means$means["TOT_PREC"]) /
      hist_means$sds["TOT_PREC"]
  )

mu_precip_year <- cary_met |>
  mutate(Date = mdy(DATE)) |>
  filter(
    Date <= "2005-12-31",
    Date >= "1995-01-01",
    TOT_PREC > min(TOT_PREC, na.rm = TRUE)
  ) |>
  mutate(year = year(Date)) |>
  group_by(year) |>
  summarise(sum.precip = sum(TOT_PREC)) |>
  pull(sum.precip) |>
  mean()

# get NOAA NMME forecast dates
nmme_files <- list.files(file.path(dir_nmme))
dates_nmme <- ymd(nmme_files)

if (remove) {
  data_tick <- data_tick |>
    mutate(Larvae = NA)
}

if (use_mice) {
  # get the mice data for the calibration period
  mna <- mna_hindcast(calibration_mice) |>
    suppressMessages()
}

# =========================================== #
#       get initial conditions
# =========================================== #

calibration_timeseries <- "dormantNymphTimeSeries.csv"
fname <- file.path(dir_data, calibration_timeseries)
df_latent <- read_csv(fname)
data_latent <- df_latent |>
  mutate(model = gsub("DormantNymph", "", model)) |>
  filter(
    model == model_job,
    type == "latent",
    statistic == "conf_50",
    ticksFrom == as.character(ticks_job),
    paramsFrom == as.character(params_job),
    ua == ua_from,
    month(DATE) == month(start_date)
  ) |>
  group_by(lifeStage) |>
  summarise(
    mu = mean(value),
    prec = 1 / var(value)
  ) |>
  pivot_wider(names_from = lifeStage, values_from = c(mu, prec))

IC <- tibble(
  mu = c(
    pull(data_latent, mu_larvae),
    pull(data_latent, mu_dormant),
    pull(data_latent, mu_nymphs),
    pull(data_latent, mu_adults)
  ),
  prec = c(
    pull(data_latent, prec_larvae),
    pull(data_latent, prec_dormant),
    pull(data_latent, prec_nymphs),
    pull(data_latent, prec_adults)
  )
) |>
  as.matrix()

# =========================================== #
#       get informative priors
# =========================================== #
calibration_parameters <- "dormantNymphParams.csv"
fname <- file.path(dir_data, calibration_parameters)
df_params <- read_csv(fname)
params_stats <- df_params |>
  filter(
    model == model_job,
    site == as.character(params_job)
  ) |>
  select(parameter, value) |>
  group_by(parameter) |>
  summarise(
    mu = mean(value),
    tau = 1 / var(value)
  )

phi_l <- get_prior("phi.l.mu")
phi_n <- get_prior("phi.n.mu")
phi_a <- get_prior("phi.a.mu")
theta_l2n <- get_prior("theta.ln")
theta_n2a <- get_prior("theta.na")
repro <- get_prior("repro.mu")
repro_mu <- repro[1]

n_beta <- params_stats |>
  filter(grepl("beta", parameter)) |>
  nrow()

pr_beta <- matrix(NA, n_beta, 2)
for (i in seq_len(n_beta)) {
  pr_beta[i, ] <- get_prior(paste0("beta[", i, "]"))
}

# get invgamma parameters
pr_sig <- df_params |>
  filter(
    model == model_job,
    site == params_job,
    grepl("sig", parameter)
  ) |>
  select(parameter, value) |>
  group_by(parameter) |>
  summarise(
    alpha = inv_gamma_mm(value)[1],
    beta = inv_gamma_mm(value)[2]
  )

# =========================================== #
#       run hindcast
# =========================================== #

hindcast_seq_tick <- sort(c(start_date, dates_tick, end_date))

if (use_nmme) {
  hindcast_seq_tick <- hindcast_seq_tick[hindcast_seq_tick >= first(dates_nmme)]
}

if (!production) {
  hindcast_seq_tick <- hindcast_seq_tick[1:10]
}

if (use_mice) {
  hindcast_seq_tick <- hindcast_seq_tick[hindcast_seq_tick < last(mna$time)]
}

# iterate ======================================================================================
t <- 1
# t <- 2
# t <- 3
for (t in seq_along(hindcast_seq_tick)) {
  fx_start_date <- hindcast_seq_tick[t]
  message("---------------------------------------------------")
  mm <- paste0(
    fx_start_date,
    " (",
    round(t / length(hindcast_seq_tick) * 100, 2),
    "%)"
  )
  message("\t", mm, "\n")

  # initialize nimble lists
  constants <- data <- list()

  if (t == 1) {
    fx_sequence <- seq.Date(fx_start_date, by = 1, length.out = horizon)
    y <- matrix(NA, 4, horizon)
  } else {
    # can just read from the last fx (t-1)
    dir_read <- file.path(
      dir_save,
      site_params,
      model_job,
      exp,
      as.character(hindcast_seq_tick[t - 1])
    )
    last_params <- read_csv(file.path(dir_read, "parameterSamples.csv")) |>
      suppressMessages()
    last_fx <- read_csv(file.path(dir_read, "tickSamples.csv")) |>
      suppressMessages()

    # get parameter posterior summary
    params_stats <- last_params |>
      pivot_longer(cols = -c(site, paramsFrom), names_to = "parameter") |>
      select(parameter, value) |>
      group_by(parameter) |>
      summarise(
        mu = mean(value),
        tau = 1 / var(value)
      )

    phi_l <- get_prior("phi.l.mu")
    phi_n <- get_prior("phi.n.mu")
    phi_a <- get_prior("phi.a.mu")
    theta_l2n <- get_prior("theta.ln")
    theta_n2a <- get_prior("theta.na")

    pr_beta <- matrix(NA, n_beta, 2)
    for (i in seq_len(n_beta)) {
      pr_beta[i, ] <- get_prior(paste0("beta[", i, "]"))
    }

    # get invgamma parameters
    pr_sig <- last_params |>
      pivot_longer(cols = -c(site, paramsFrom), names_to = "parameter") |>
      filter(grepl("sig", parameter)) |>
      select(parameter, value) |>
      group_by(parameter) |>
      summarise(
        alpha = inv_gamma_mm(value)[1],
        beta = inv_gamma_mm(value)[2]
      )

    ticks_stats <- last_fx |>
      group_by(lifeStage, time) |>
      summarise(
        mu = mean(value),
        tau = 1 / var(value)
      )

    obs_index <- which(dates_tick == fx_start_date)

    obs <- data_tick |>
      filter(Date %in% dates_tick[obs_index])

    y <- matrix(NA, 4, horizon)
    y[1, 1] <- obs |> pull(Larvae)
    y[3, 1] <- obs |> pull(Nymphs)
    y[4, 1] <- obs |> pull(Adults)

    fx_mu <- ticks_stats |>
      select(-tau) |>
      arrange(time) |>
      pivot_wider(
        names_from = time,
        values_from = mu
      )

    fx_tau <- ticks_stats |>
      select(-mu) |>
      arrange(time) |>
      pivot_wider(
        names_from = time,
        values_from = tau
      )

    IC <- create_ic(fx_mu, fx_tau, fx_start_date)

    fx_sequence <- seq.Date(fx_start_date, by = 1, length.out = horizon)
    n_days <- length(fx_sequence)
  }

  # get the weather data for the forecast period
  nmme_cens <- FALSE
  if (use_nmme) {
    nmme_cens <- TRUE
    muf_missing <- FALSE
    nmme2get <- floor_date(fx_start_date, unit = "month")
    nmme_name <- gsub("-", "", nmme2get)

    # check if nmme file exists (some file don't)
    if (!nmme_name %in% nmme_files) {
      nmme2get <- floor_date(fx_start_date - 30, unit = "month")
      nmme_name <- gsub("-", "", nmme2get)
    }

    df_nmme <-
      read_csv(file.path(
        # dir.top,
        # dir_data,
        dir_nmme,
        nmme_name,
        "CARY",
        paste0("NOAANMME_", nmme_name, ".csv")
      )) |>
      suppressMessages()

    horizon <- as.numeric(max(unique(df_nmme$time)) - fx_start_date)
    fx_sequence <- seq.Date(fx_start_date, by = 1, length.out = horizon)

    y <- matrix(NA, 4, horizon)
    y[1, 1] <- obs |> pull(Larvae)
    y[3, 1] <- obs |> pull(Nymphs)
    y[4, 1] <- obs |> pull(Adults)

    target_nmme <- df_nmme |>
      select(-tmin) |>
      filter(time %in% fx_sequence) |>
      mutate(
        rhmin = rhmin * 100,
        rhmax = rhmax * 100
      )

    nmme_means <- target_nmme |>
      pivot_longer(
        cols = c(tmax, precipitation, rhmin, rhmax),
        names_to = "variable"
      ) |>
      group_by(variable, ensemble) |>
      summarise(mu = mean(value)) |>
      pivot_wider(
        names_from = variable,
        values_from = mu
      ) |>
      mutate(
        bias.tmax = hist_means$means["MAX_TEMP"] - tmax,
        bias.rhmax = hist_means$means["MAX_RH"] - rhmax,
        bias.rhmin = hist_means$means["MIN_RH"] - rhmin
      ) |>
      ungroup()

    precip_year <- target_nmme |>
      select(time, ensemble, precipitation) |>
      group_by(ensemble) |>
      summarise(tot.prec = sum(precipitation)) |>
      mutate(
        yearly.bias = mu_precip_year - tot.prec,
        day.add = yearly.bias / 365
      ) |>
      ungroup()

    ens_vec <- target_nmme$ensemble |> unique()
    nmme_correct <- tibble()
    for (i in ens_vec) {
      psub <- target_nmme |>
        filter(ensemble == i) |>
        mutate(
          precipitation = precipitation + precip_year$day.add[i],
          tmax = tmax + nmme_means$bias.tmax[i],
          rhmax = rhmax + nmme_means$bias.rhmax[i],
          rhmin = rhmin + nmme_means$bias.rhmin[i]
        )
      nmme_correct <- bind_rows(nmme_correct, psub)
    }

    # scale to historical period
    nmme_scaled <- nmme_correct |>
      filter(time %in% fx_sequence) |>
      select(time, ensemble, tmax, rhmax, rhmin, precipitation) |>
      mutate(
        tmax = (tmax - hist_means$means["MAX_TEMP"]) /
          hist_means$sds["MAX_TEMP"],
        rhmax = (rhmax - hist_means$means["MAX_RH"]) / hist_means$sds["MAX_RH"],
        rhmin = (rhmin - hist_means$means["MIN_RH"]) / hist_means$sds["MIN_RH"],
        precipitation = (precipitation - hist_means$means["TOT_PREC"]) /
          hist_means$sds["TOT_PREC"]
      )

    interval_x <- matrix(NA, 4, 2)
    interval_x[1, ] <- c(
      hist_means$scale.min["MAX_TEMP"],
      hist_means$scale.max["MAX_TEMP"]
    )
    interval_x[2, ] <- c(
      hist_means$scale.min["MAX_RH"],
      hist_means$scale.max["MAX_RH"]
    )
    interval_x[3, ] <- c(
      hist_means$scale.min["MIN_RH"],
      hist_means$scale.max["MIN_RH"]
    )
    interval_x[4, ] <- c(
      hist_means$scale.min["TOT_PREC"],
      hist_means$scale.max["TOT_PREC"]
    )

    y_ind <- y_censored <- array(NA, dim = c(10, 4, horizon))
    for (i in seq_along(fx_sequence)) {
      X <- nmme_scaled |>
        filter(time == fx_sequence[i]) |>
        select(tmax, rhmax, rhmin, precipitation) |>
        as.matrix()

      x_ind <- x_censored <- matrix(NA, nrow(X), ncol(X))
      for (j in seq_len(ncol(X))) {
        for (n in seq_len(nrow(X))) {
          if (is.na(X[n, j])) {
            X[n, j] <- rnorm(1)
          }

          if (j != 4) {
            if (X[n, j] < interval_x[j, 1]) {
              x_ind[n, j] <- 0
              x_censored[n, j] <- 0
            } else if (X[n, j] > interval_x[j, 2]) {
              x_ind[n, j] <- 0
              x_censored[n, j] <- 1
            } else {
              x_ind[n, j] <- 1
              x_censored[n, j] <- X[n, j]
            }
          } else {
            x_ind[n, j] <- as.numeric(X[n, j] > 0)
            x_censored[n, j] <- as.numeric(ifelse(
              X[n, j] > interval_x[j, 2],
              0,
              X[n, j]
            ))
          }
        }
      }

      y_ind[,, i] <- x_ind
      y_censored[,, i] <- x_censored
    }

    years_needed <- df_nmme |>
      select(time) |>
      mutate(year = year(time)) |>
      pull(year) |>
      unique() |>
      min()

    # get the cgd for the day before the forecast, apply to cgd calculation from nmme
    cary_cgdd <- weather_corrected |>
      filter(Date == fx_start_date - 1) |>
      pull(cgdd)

    cgdd <- nmme_correct |>
      filter(time %in% fx_sequence) |>
      mutate(year = year(time)) |>
      group_by(ensemble, year) |>
      mutate(
        gd = if_else(tmax > 10, tmax - 10, 0),
        cgdd = cumsum(gd)
      ) |>
      ungroup() |>
      mutate(cgdd = if_else(year == years_needed, cgdd + cary_cgdd, cgdd)) |>
      select(time, ensemble, cgdd) |>
      pivot_wider(names_from = time, values_from = cgdd) |>
      select(-ensemble)

    data$cgdd.mu <- apply(cgdd, 2, mean)
    data$cgdd.tau <- 1 / apply(cgdd, 2, var)
    data$y_ind <- y_ind
    data$y_censored <- y_censored
    data$mu_0 <- rep(0, ncol(X))
    data$lambda_0 <- solve(diag(100, ncol(X)))
    data$nu_0 <- ncol(X) + 1
    data$wts <- rep(1, nrow(X)) * nrow(X)
    data$Sigma_0 <- solve(diag(100, ncol(X)))
    data$interval <- interval_x

    constants$J <- ncol(X)
    constants$K <- nrow(X)
  } else {
    weather_hind <- weather_corrected |>
      ungroup() |>
      filter(Date %in% fx_sequence)

    cgdd <- weather_hind |>
      pull(cgdd)

    muf <- weather_hind |>
      select(MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) |>
      as.matrix()

    # need a prior on any missing weather from Cary
    muf_missing <- FALSE
    if (any(is.na(muf))) {
      muf_missing <- TRUE
      n_missing <- length(which(is.na(muf)))
      t_missing <- var_missing <- numeric(n_missing + 1) # need to add one to make a vector for single missing points
      xx <- 1
      for (i in seq_len(nrow(muf))) {
        for (j in seq_len(ncol(muf))) {
          if (is.na(muf[i, j])) {
            t_missing[xx] <- i
            var_missing[xx] <- j
            xx <- xx + 1
          }
        }
      }
      constants$n_missing <- n_missing
      constants$t_missing <- t_missing
      constants$var_missing <- var_missing
    }

    data$gdd <- cgdd
    data$muf <- muf
  }

  # finalize data
  data$y <- y
  data$IC <- IC
  data$pr_phi_l <- phi_l
  data$pr_phi_n <- phi_n
  data$pr_phi_a <- phi_a
  data$pr_theta_l2n <- theta_l2n
  data$pr_theta_n2a <- theta_n2a
  data$repro_mu <- repro_mu
  data$pr_beta <- pr_beta
  data$pr_sig <- pr_sig |>
    select(-parameter) |>
    as.matrix()

  if (use_mice) {
    mice_sub <- mna |>
      filter(time %in% fx_sequence) |>
      pull(mna)
    data$mice <- mice_sub
  }

  if (year(fx_start_date) == 2021) {
    # | length(mice_sub) < length(fx_sequence)){
    horizon <- min(length(cgdd), length(mice_sub))
    data$y <- y[, 1:horizon]
  }

  # finalize constants
  constants$n_beta <- n_beta
  constants$horizon <- horizon
  constants$ns <- 4

  # build inits
  inits <- function() {
    list(
      phi.l.mu = rnorm(1, phi_l[1], 1 / sqrt(phi_l[2])),
      phi.n.mu = rnorm(1, phi_n[1], 1 / sqrt(phi_n[2])),
      phi.a.mu = rnorm(1, phi_a[1], 1 / sqrt(phi_a[2])),
      theta.ln = rnorm(1, theta_l2n[1], 1 / sqrt(theta_l2n[2])),
      theta.na = rnorm(1, theta_n2a[1], 1 / sqrt(theta_n2a[2])),
      beta = rnorm(n_beta, pr_beta[, 1], 1 / sqrt(pr_beta[, 2])),
      sig = rinvgamma(4, pr_sig$alpha, pr_sig$beta),
      px = matrix(rpois(4 * horizon, 5), 4, horizon),
      x = matrix(rpois(4 * horizon, 5), 4, horizon)
    )
  }

  cl <- makeCluster(n_slots)
  out_nchains <- run_hindcast_nimble(
    cl = cl,
    code = model_code,
    data = data,
    constants = constants,
    inits = inits,
    n_iter = n_iter,
    thin = 5,
    nmme_cens = nmme_cens,
    mu_f_missing = muf_missing,
    use_mice = use_mice
  )
  stopCluster(cl)

  dat_hindcast <- do.call(rbind, out_nchains)

  nodes <- colnames(out_nchains[[1]])
  which_tick <- grep("x", nodes)
  which_params <- which(nodes[-grep("x", nodes)] %in% nodes)
  gelman_keep <- numeric(length(nodes))
  for (ff in seq_along(nodes)) {
    mcmc_check <- list()
    col <- nodes[ff]
    for (c in seq_along(out_nchains)) {
      mcmc_check[[c]] <- coda::mcmc(out_nchains[[c]][, col])
    }
    gelman_keep[ff] <- try(coda::gelman.diag(mcmc_check, transform = TRUE)$psrf[
      1
    ])
  }

  if (any(gelman_keep > 1.2)) {
    message("WARNING: Convergence not reached!")
    bad_nodes <- which(gelman_keep > 1.2)
    bad_params <- tibble(node = nodes[bad_nodes], psrf = gelman_keep[bad_nodes])
    print(head(as.matrix(bad_params)))
  } else {
    message("Convergence = TRUE")
  }

  draws <- sample.int(nrow(dat_hindcast), nmc, replace = TRUE)
  dat_draws <- dat_hindcast[draws, ]

  # save hindcast
  save_params <- dat_draws[, which_params] |>
    as_tibble() |>
    mutate(
      site = as.character(ticks_job),
      paramsFrom = as.character(params_job)
    )

  if (use_nmme) {
    which_muf <- grep("muf[", colnames(save_params), fixed = TRUE)
    which_pf <- grep("pf[", colnames(save_params), fixed = TRUE)
    save_muf <- save_params |> select(all_of(c(which_muf)), site, paramsFrom)
    save_pf <- save_params |> select(all_of(c(which_pf)), site, paramsFrom)
    save_params <- save_params |> select(-all_of(c(which_pf, which_muf)))
  }

  dates_tb <- tibble(
    fx_index = seq_along(fx_sequence),
    time = fx_sequence
  )

  ticks_long <- dat_draws[, which_tick] |>
    as_tibble() |>
    pivot_longer(cols = everything()) |>
    mutate(
      lifeStage = if_else(grepl("x[1", name, fixed = TRUE), "larvae", "x"),
      lifeStage = if_else(
        grepl("x[2", name, fixed = TRUE),
        "dormant",
        lifeStage
      ),
      lifeStage = if_else(
        grepl("x[3", name, fixed = TRUE),
        "nymphs",
        lifeStage
      ),
      lifeStage = if_else(
        grepl("x[4", name, fixed = TRUE),
        "adults",
        lifeStage
      ),
      fx_index = as.numeric(str_extract(name, "\\d*(?=\\])")),
      row = seq_len(n())
    )

  save_ticks <- left_join(ticks_long, dates_tb, by = "fx_index") |>
    mutate(
      site = as.character(ticks_job),
      paramsFrom = as.character(params_job)
    )

  save_ticks |>
    group_by(time, lifeStage) |>
    summarise(
      low = quantile(value, 0.025),
      med = quantile(value, 0.5),
      high = quantile(value, 0.975)
    ) |>
    # filter(lifeStage == "nymphs") |>
    ggplot() +
    aes(x = time) +
    geom_ribbon(
      aes(
        ymin = low,
        ymax = high
      ),
      fill = "lightblue"
    ) +
    geom_line(aes(y = med)) +
    facet_wrap(~lifeStage, scales = "free") +
    theme_bw()

  site_params <- paste0("ticksFrom_", ticks_job, "_paramsFrom_", params_job)
  exp <- if_else(use_nmme, paste0(experiment, "_nmme"), experiment)
  dir_write <- file.path(
    dir_save,
    site_params,
    model_job,
    exp,
    as.character(fx_start_date)
  )
  if (!dir.exists(dir_write)) {
    dir.create(dir_write, showWarnings = FALSE, recursive = TRUE)
  }

  write_csv(save_params, file = file.path(dir_write, "parameterSamples.csv"))
  write_csv(save_ticks, file = file.path(dir_write, "tickSamples.csv"))

  if (use_nmme) {
    write_csv(save_muf, file = file.path(dir_write, "mufSamples.csv"))
    write_csv(save_pf, file = file.path(dir_write, "pfSamples.csv"))
  }

  message("\n")
}
