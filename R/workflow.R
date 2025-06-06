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
source("Functions/run_hindcast_nimble.R")

n.slots <- Sys.getenv("NSLOTS") %>% as.numeric()
use.nmme <- TRUE
production <- TRUE
n.iter <- 20000
Nmc <- 2000
horizon <- 365
use.mice <- TRUE
restart <- FALSE

dir.data <- "Data"
dir.tick <- "FinalOut/A_Correct/DormantNymphs/DormantNymph_met_and_mice_nimble_All_Global"
dir.mice <- "FinalOut"
dir.nmme <- "C:/Users/John.Foster/Downloads/NMME"
dir.analysis <- file.path("Analysis/dormantStages/nimbleModels")
dir.save <- file.path(dir.mice, "Chapter2")

# field sites
sites <- c("Green", "Henry", "Tea")

# all possible combinations for hindcast experiments
all.jobs <- expand.grid(
  paramsFrom = sites,
  ticksFrom = sites,
  remove = c(TRUE, FALSE)
)

# all.jobs <- all.jobs %>%
#   slice(-c(10,14,18))

job.num <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(job.num)) job.num <- 9
remove <- all.jobs$remove[job.num]
paramsJob <- all.jobs$paramsFrom[job.num]
ticksJob <- all.jobs$ticksFrom[job.num]
miceJob <- ticksJob
modelJob <- if_else(use.mice, "WithWeatherAndMiceGlobal", "Weather")
uaFrom <- if_else(use.mice, "mice_ic_parameter_process", "ic_parameter_process")
calibration.mice <- paste(miceJob, "Control")
calibration.ticks <- paste(ticksJob, "Control")
experiment <- "base"
if (use.mice) experiment <- paste(experiment, "mna", sep = "_")
if (remove) experiment <- paste(experiment, "remove", sep = "_")
if (use.nmme) experiment <- paste(experiment, "nmme", sep = "_")

# =========================================== #
#       1 - state data intake
# =========================================== #

df.tick <- read_csv(file.path(dir.data, "Ticks2006_2021.csv"))
df.mice <- read_csv(file.path(dir.data, "Mice2006_2021.csv")) # currently tick data only through 2020

data.tick <- df.tick %>% filter(Grid == calibration.ticks)
data.mice <- df.mice %>% filter(Grid == calibration.mice)

dates.tick <- data.tick$Date
dates.mice <- data.mice$Date

start.date <- if_else(use.mice, min(dates.mice), min(dates.tick) - 31) # start one month before first tick observation
end.date <- max(dates.tick) # end at last tick observation
hindcast.period <- seq.Date(start.date, end.date, by = 1)

# hindcast.period represents every single day in the time series
# for a majority of the hindcasts, we only need to run the simulation
# after each observation because the weather is treated as known, then
# run hindcasts for every day once we reach the GEFS period (2020+)
# also only need to run the hindcasts during the sampling season (not winter)

# =========================================== #
#       1 - weather data intake
# =========================================== #

# we will need the means/SDs from the calibration period
source("Functions/scale_met_forecast.R")
hist.means <- scale_met_forecast()

# get the observed data from Cary and center to historical mean
cary.met <- read_csv(file.path(dir.data, "Cary_Met_Data_Daily.csv")) %>% suppressWarnings()
met.pre.2021 <- cary.met %>%
  slice(-c((nrow(cary.met) - 364):nrow(cary.met))) %>%
  mutate(Date = mdy(DATE)) %>%
  filter(Date >= "1995-01-01")
met.2021 <- cary.met %>%
  slice((nrow(cary.met) - 364):nrow(cary.met)) %>%
  mutate(Date = ymd(DATE))

weather.corrected <- bind_rows(met.pre.2021, met.2021) %>%
  select(Date, MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) %>%
  arrange(Date) %>%
  mutate(year = year(Date)) %>%
  group_by(year) %>%
  mutate(
    gd = if_else(MAX_TEMP > 10, MAX_TEMP - 10, 0),
    cgdd = cumsum(gd),
    MAX_TEMP = (MAX_TEMP - hist.means$means["MAX_TEMP"]) / hist.means$sds["MAX_TEMP"],
    MAX_RH = (MAX_RH - hist.means$means["MAX_RH"]) / hist.means$sds["MAX_RH"],
    MIN_RH = (MIN_RH - hist.means$means["MIN_RH"]) / hist.means$sds["MIN_RH"],
    TOT_PREC = (TOT_PREC - hist.means$means["TOT_PREC"]) / hist.means$sds["TOT_PREC"]
  )

mu.precip.year <- cary.met %>%
  mutate(Date = mdy(DATE)) %>%
  filter(
    Date <= "2005-12-31",
    Date >= "1995-01-01",
    TOT_PREC > min(TOT_PREC, na.rm = TRUE)
  ) %>%
  mutate(year = year(Date)) %>%
  group_by(year) %>%
  summarise(sum.precip = sum(TOT_PREC)) %>%
  pull(sum.precip) %>%
  mean()

# get NOAA NMME forecast dates
nmme.files <- list.files(file.path(dir.nmme))
dates.nmme <- ymd(nmme.files)

# which observations fall within the larval phenological window and are not peak nymph for that year?
if (remove) {
  data.tick <- data.tick %>%
    mutate(Larvae = NA)
}

if (use.mice) {
  source("Functions/mna_hindcast.R")
  mna <- mna_hindcast(calibration.mice) %>%
    suppressMessages()
}

# =========================================== #
#       get initial conditions
# =========================================== #

df.latent <- read_csv("Analysis/dormantNymphTimeSeries.csv")
data.latent <- df.latent %>%
  mutate(model = gsub("DormantNymph", "", model)) %>%
  filter(
    model == modelJob,
    type == "latent",
    statistic == "conf_50",
    ticksFrom == as.character(ticksJob),
    paramsFrom == as.character(paramsJob),
    ua == uaFrom,
    month(DATE) == month(start.date)
  ) %>%
  group_by(lifeStage) %>%
  summarise(
    mu = mean(value),
    prec = 1 / var(value)
  ) %>%
  pivot_wider(names_from = lifeStage, values_from = c(mu, prec))

IC <- tibble(
  mu = c(
    pull(data.latent, mu_larvae),
    pull(data.latent, mu_dormant),
    pull(data.latent, mu_nymphs),
    pull(data.latent, mu_adults)
  ),
  prec = c(
    pull(data.latent, prec_larvae),
    pull(data.latent, prec_dormant),
    pull(data.latent, prec_nymphs),
    pull(data.latent, prec_adults)
  )
) %>%
  as.matrix()

# =========================================== #
#       get informative priors
# =========================================== #


df.params <- read_csv("Analysis/dormantNymphParams.csv")
params.stats <- df.params %>%
  filter(
    model == modelJob,
    site == as.character(paramsJob)
  ) %>%
  select(parameter, value) %>%
  group_by(parameter) %>%
  summarise(
    mu = mean(value),
    tau = 1 / var(value)
  )

get_prior <- function(name) {
  pr <- numeric(2)
  xx <- params.stats %>%
    filter(parameter == name)
  pr[1] <- xx %>% pull(mu)
  pr[2] <- xx %>% pull(tau)
  pr
}

phi.l <- get_prior("phi.l.mu")
phi.n <- get_prior("phi.n.mu")
phi.a <- get_prior("phi.a.mu")
theta.l2n <- get_prior("theta.ln")
theta.n2a <- get_prior("theta.na")
repro <- get_prior("repro.mu")
repro.mu <- repro[1]

n.beta <- params.stats %>%
  filter(grepl("beta", parameter)) %>%
  nrow()

pr.beta <- matrix(NA, n.beta, 2)
for (i in seq_len(n.beta)) {
  pr.beta[i, ] <- get_prior(paste0("beta[", i, "]"))
}

# function to approximate moment the inverse gamma
inv_gamma_mm <- function(x) {
  mu <- mean(x)
  v <- var(x)
  alpha <- (mu^2 / v) + 2
  beta <- mu * ((mu^2 / v) + 1)
  return(c("alpha" = alpha, "beta" = beta))
}

# get invgamma parameters
pr.sig <- df.params %>%
  filter(
    model == modelJob,
    site == paramsJob,
    grepl("sig", parameter)
  ) %>%
  select(parameter, value) %>%
  group_by(parameter) %>%
  summarise(
    alpha = inv_gamma_mm(value)[1],
    beta = inv_gamma_mm(value)[2]
  )

# =========================================== #
#       run hindcast
# =========================================== #

# need to do the analysis from last observation to current observation
# need two y's

obs.only.tick <- dates.tick[dates.tick <= first(dates.nmme)]
obs.only.mice <- dates.mice[dates.mice <= first(dates.nmme)]
nmme.seq <- seq.Date(first(dates.nmme), last(dates.nmme), by = 1)

hindcast.seq.tick <- sort(c(start.date, dates.tick, end.date))
hindcast.seq.mice <- sort(c(start.date, dates.mice, end.date))

if (use.nmme) {
  hindcast.seq.tick <- hindcast.seq.tick[hindcast.seq.tick >= first(dates.nmme)]
  hindcast.seq.mice <- hindcast.seq.mice[hindcast.seq.mice >= first(dates.nmme)]
}

if (!production) {
  hindcast.seq.tick <- hindcast.seq.tick[1:10]
}

if (use.mice) {
  hindcast.seq.tick <- hindcast.seq.tick[hindcast.seq.tick < last(mna$time)]
}


#

# iterate ======================================================================================
# hindcast.seq.tick <- hindcast.seq.tick[hindcast.seq.tick > "2019-09-05"]
for (t in seq_along(hindcast.seq.tick)) {
  # for(t in 1:4){
  fx.start.date <- hindcast.seq.tick[t]
  message("---------------------------------------------------")
  mm <- paste0(fx.start.date, " (", round(t / length(hindcast.seq.tick) * 100, 2), "%)")
  message("\t", mm, "\n")

  # initialize nimble lists
  constants <- data <- list()

  if (t == 1 & !use.nmme) {
    fx.sequence <- seq.Date(fx.start.date, by = 1, length.out = horizon)
    # n.days <- horizon
    y <- matrix(NA, 4, horizon)
  } else {
    # read last forecast parameters and states
    if (t == 1) {
      site.params <- paste0("ticksFrom_", ticksJob, "_paramsFrom_", paramsJob)
      exp.read <- gsub("_nmme", "", experiment) # use last forecast from related experiment
      dir.read <- file.path(dir.save, site.params, modelJob, exp.read)
      fx.files <- list.files(dir.read) %>% ymd()
      fx.read <- fx.files[max(which(fx.files <= fx.start.date))]

      last.fx <- read_csv(file.path(dir.read, as.character(fx.read), "tickSamples.csv")) %>%
        suppressMessages()
      last.params <- read_csv(file.path(dir.read, as.character(fx.read), "parameterSamples.csv")) %>%
        suppressMessages()
    } else { # can just read from the last fx (t-1)
      dir.read <- file.path(dir.save, site.params, modelJob, exp, as.character(hindcast.seq.tick[t - 1]))
      last.params <- read_csv(file.path(dir.read, "parameterSamples.csv")) %>%
        suppressMessages()
      last.fx <- read_csv(file.path(dir.read, "tickSamples.csv")) %>%
        suppressMessages()
    }

    # get parameter posterior summary
    params.stats <- last.params %>%
      pivot_longer(cols = -c(site, paramsFrom), names_to = "parameter") %>%
      select(parameter, value) %>%
      group_by(parameter) %>%
      summarise(
        mu = mean(value),
        tau = 1 / var(value)
      )

    phi.l <- get_prior("phi.l.mu")
    phi.n <- get_prior("phi.n.mu")
    phi.a <- get_prior("phi.a.mu")
    theta.l2n <- get_prior("theta.ln")
    theta.n2a <- get_prior("theta.na")
    # repro <- get_prior("repro.mu")

    pr.beta <- matrix(NA, n.beta, 2)
    for (i in seq_len(n.beta)) {
      pr.beta[i, ] <- get_prior(paste0("beta[", i, "]"))
    }

    # get invgamma parameters
    pr.sig <- last.params %>%
      pivot_longer(cols = -c(site, paramsFrom), names_to = "parameter") %>%
      filter(grepl("sig", parameter)) %>%
      select(parameter, value) %>%
      group_by(parameter) %>%
      summarise(
        alpha = inv_gamma_mm(value)[1],
        beta = inv_gamma_mm(value)[2]
      )

    tick.stats <- last.fx %>%
      group_by(lifeStage, time) %>%
      summarise(
        mu = mean(value),
        tau = 1 / var(value)
      )

    obs.index <- which(dates.tick == fx.start.date)

    obs <- data.tick %>%
      filter(Date %in% dates.tick[obs.index])

    y <- matrix(NA, 4, horizon)
    y[1, 1] <- obs %>% pull(Larvae)
    y[3, 1] <- obs %>% pull(Nymphs)
    y[4, 1] <- obs %>% pull(Adults)

    fx.mu <- tick.stats %>%
      select(-tau) %>%
      arrange(time) %>%
      pivot_wider(
        names_from = time,
        values_from = mu
      )

    fx.tau <- tick.stats %>%
      select(-mu) %>%
      arrange(time) %>%
      pivot_wider(
        names_from = time,
        values_from = tau
      )

    IC <- matrix(NA, 4, 2)
    IC[1, 1] <- fx.mu %>%
      filter(lifeStage == "larvae") %>%
      pull(as.character(fx.start.date))
    IC[1, 2] <- fx.tau %>%
      filter(lifeStage == "larvae") %>%
      pull(as.character(fx.start.date))
    IC[2, 1] <- fx.mu %>%
      filter(lifeStage == "dormant") %>%
      pull(as.character(fx.start.date))
    IC[2, 2] <- fx.tau %>%
      filter(lifeStage == "dormant") %>%
      pull(as.character(fx.start.date))
    IC[3, 1] <- fx.mu %>%
      filter(lifeStage == "nymphs") %>%
      pull(as.character(fx.start.date))
    IC[3, 2] <- fx.tau %>%
      filter(lifeStage == "nymphs") %>%
      pull(as.character(fx.start.date))
    IC[4, 1] <- fx.mu %>%
      filter(lifeStage == "adults") %>%
      pull(as.character(fx.start.date))
    IC[4, 2] <- fx.tau %>%
      filter(lifeStage == "adults") %>%
      pull(as.character(fx.start.date))

    fx.sequence <- seq.Date(fx.start.date, by = 1, length.out = horizon)
    n.days <- length(fx.sequence)
  }

  nmme.cens <- FALSE
  if (use.nmme) {
    nmme.cens <- TRUE
    muf.missing <- FALSE
    nmme2get <- floor_date(fx.start.date, unit = "month")
    nmme.name <- gsub("-", "", nmme2get)

    # check if nmme file exists (some file don't)
    if (!nmme.name %in% nmme.files) {
      nmme2get <- floor_date(fx.start.date - 30, unit = "month")
      nmme.name <- gsub("-", "", nmme2get)
    }

    df.nmme <-
      read_csv(file.path(
        # dir.top,
        # dir.data,
        dir.nmme,
        nmme.name,
        "CARY",
        paste0("NOAANMME_", nmme.name, ".csv")
      )) %>%
      suppressMessages()

    horizon <- as.numeric(max(unique(df.nmme$time)) - fx.start.date)
    fx.sequence <- seq.Date(fx.start.date, by = 1, length.out = horizon)

    y <- matrix(NA, 4, horizon)
    y[1, 1] <- obs %>% pull(Larvae)
    y[3, 1] <- obs %>% pull(Nymphs)
    y[4, 1] <- obs %>% pull(Adults)

    target.nmme <- df.nmme %>%
      select(-tmin) %>%
      filter(time %in% fx.sequence) %>%
      mutate(
        rhmin = rhmin * 100,
        rhmax = rhmax * 100
      )

    nmme.means <- target.nmme %>%
      pivot_longer(
        cols = c(tmax, precipitation, rhmin, rhmax),
        names_to = "variable"
      ) %>%
      group_by(variable, ensemble) %>%
      summarise(mu = mean(value)) %>%
      pivot_wider(
        names_from = variable,
        values_from = mu
      ) %>%
      mutate(
        bias.tmax = hist.means$means["MAX_TEMP"] - tmax,
        bias.rhmax = hist.means$means["MAX_RH"] - rhmax,
        bias.rhmin = hist.means$means["MIN_RH"] - rhmin
      ) %>%
      ungroup()

    precip.year <- target.nmme %>%
      select(time, ensemble, precipitation) %>%
      group_by(ensemble) %>%
      summarise(tot.prec = sum(precipitation)) %>%
      mutate(
        yearly.bias = mu.precip.year - tot.prec,
        day.add = yearly.bias / 365
      ) %>%
      ungroup()

    ens.vec <- target.nmme$ensemble %>% unique()
    nmme.correct <- tibble()
    for (i in ens.vec) {
      psub <- target.nmme %>%
        filter(ensemble == i) %>%
        mutate(
          precipitation = precipitation + precip.year$day.add[i],
          tmax = tmax + nmme.means$bias.tmax[i],
          rhmax = rhmax + nmme.means$bias.rhmax[i],
          rhmin = rhmin + nmme.means$bias.rhmin[i]
        )
      nmme.correct <- bind_rows(nmme.correct, psub)
    }

    # nmme.correct %>%
    #   pivot_longer(cols = c(tmax, precipitation, rhmin, rhmax),
    #                names_to = "variable") %>%
    #   mutate(ensemble = as.character(ensemble)) %>%
    #   # filter(variable == "rhmax") %>%
    #   ggplot() +
    #   aes(x = time, y = value, color = ensemble) +
    #   geom_line() +
    #   facet_wrap(~ variable) +
    #   theme_bw()

    # scale to historical period
    nmme.scaled <- nmme.correct %>%
      filter(time %in% fx.sequence) %>%
      select(time, ensemble, tmax, rhmax, rhmin, precipitation) %>%
      mutate(
        tmax = (tmax - hist.means$means["MAX_TEMP"]) / hist.means$sds["MAX_TEMP"],
        rhmax = (rhmax - hist.means$means["MAX_RH"]) / hist.means$sds["MAX_RH"],
        rhmin = (rhmin - hist.means$means["MIN_RH"]) / hist.means$sds["MIN_RH"],
        precipitation = (precipitation - hist.means$means["TOT_PREC"]) / hist.means$sds["TOT_PREC"]
      )

    intervalX <- matrix(NA, 4, 2)
    intervalX[1, ] <- c(hist.means$scale.min["MAX_TEMP"], hist.means$scale.max["MAX_TEMP"])
    intervalX[2, ] <- c(hist.means$scale.min["MAX_RH"], hist.means$scale.max["MAX_RH"])
    intervalX[3, ] <- c(hist.means$scale.min["MIN_RH"], hist.means$scale.max["MIN_RH"])
    intervalX[4, ] <- c(hist.means$scale.min["TOT_PREC"], hist.means$scale.max["TOT_PREC"])

    y.ind <- y.censored <- array(NA, dim = c(10, 4, horizon))
    for (i in seq_along(fx.sequence)) {
      X <- nmme.scaled %>%
        filter(time == fx.sequence[i]) %>%
        select(tmax, rhmax, rhmin, precipitation) %>%
        as.matrix()

      x.ind <- x.censored <- matrix(NA, nrow(X), ncol(X))
      for (j in 1:ncol(X)) {
        for (n in 1:nrow(X)) {
          if (is.na(X[n, j])) {
            X[n, j] <- rnorm(1)
          }

          if (j != 4) {
            if (X[n, j] < intervalX[j, 1]) {
              x.ind[n, j] <- 0
              x.censored[n, j] <- 0
            } else if (X[n, j] > intervalX[j, 2]) {
              x.ind[n, j] <- 0
              x.censored[n, j] <- 1
            } else {
              x.ind[n, j] <- 1
              x.censored[n, j] <- X[n, j]
            }
          } else {
            x.ind[n, j] <- as.numeric(X[n, j] > 0)
            x.censored[n, j] <- as.numeric(ifelse(X[n, j] > intervalX[j, 2], 0, X[n, j]))
          }
        }
      }

      y.ind[, , i] <- x.ind
      y.censored[, , i] <- x.censored
    }

    years.needed <- df.nmme %>%
      select(time) %>%
      mutate(year = year(time)) %>%
      pull(year) %>%
      unique() %>%
      min()

    # get the cgd for the day before the forecast, apply to cgd calculation from nmme
    cary.cgdd <- weather.corrected %>%
      filter(Date == fx.start.date - 1) %>%
      pull(cgdd)

    cgdd <- nmme.correct %>%
      filter(time %in% fx.sequence) %>%
      mutate(year = year(time)) %>%
      group_by(ensemble, year) %>%
      mutate(
        gd = if_else(tmax > 10, tmax - 10, 0),
        cgdd = cumsum(gd)
      ) %>%
      ungroup() %>%
      mutate(cgdd = if_else(year == years.needed, cgdd + cary.cgdd, cgdd)) %>%
      select(time, ensemble, cgdd) %>%
      pivot_wider(names_from = time, values_from = cgdd) %>%
      select(-ensemble)

    data$cgdd.mu <- apply(cgdd, 2, mean)
    data$cgdd.tau <- 1 / apply(cgdd, 2, var)
    data$y.ind <- y.ind
    data$y.censored <- y.censored
    data$mu_0 <- rep(0, ncol(X))
    data$lambda_0 <- solve(diag(100, ncol(X)))
    data$nu_0 <- ncol(X) + 1
    data$wts <- rep(1, nrow(X)) * nrow(X)
    data$Sigma_0 <- solve(diag(100, ncol(X)))
    data$interval <- intervalX

    constants$J <- ncol(X)
    constants$K <- nrow(X)
  } else {
    weather.hind <- weather.corrected %>%
      ungroup() %>%
      filter(Date %in% fx.sequence)

    cgdd <- weather.hind %>%
      pull(cgdd)

    muf <- weather.hind %>%
      select(MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) %>%
      as.matrix()

    # need a prior on any missing weather from Cary
    muf.missing <- FALSE
    if (any(is.na(muf))) {
      muf.missing <- TRUE
      n.missing <- length(which(is.na(muf)))
      t.missing <- var.missing <- numeric(n.missing + 1) # need to add one to make a vector for single missing points
      xx <- 1
      for (i in 1:nrow(muf)) {
        for (j in 1:ncol(muf)) {
          if (is.na(muf[i, j])) {
            t.missing[xx] <- i
            var.missing[xx] <- j
            xx <- xx + 1
          }
        }
      }
      constants$n.missing <- n.missing
      constants$t.missing <- t.missing
      constants$var.missing <- var.missing
    }

    data$gdd <- cgdd
    data$muf <- muf
  }

  # finalize data
  data$y <- y
  data$IC <- IC
  data$pr.phi.l <- phi.l
  data$pr.phi.n <- phi.n
  data$pr.phi.a <- phi.a
  data$pr.theta.l2n <- theta.l2n
  data$pr.theta.n2a <- theta.n2a
  # data$pr.repro <- repro
  data$repro.mu <- repro.mu
  data$pr.beta <- pr.beta
  data$pr.sig <- pr.sig %>%
    select(-parameter) %>%
    as.matrix()

  if (use.mice) {
    mice.sub <- mna %>%
      filter(time %in% fx.sequence) %>%
      pull(mna)
    data$mice <- mice.sub
  }

  if (year(fx.start.date) == 2021) { # | length(mice.sub) < length(fx.sequence)){
    horizon <- min(length(cgdd), length(mice.sub))
    data$y <- y[, 1:horizon]
  }

  # finalize constants
  constants$n.beta <- n.beta
  constants$horizon <- horizon
  constants$ns <- 4

  # build inits
  inits <- function() {
    list(
      phi.l.mu = rnorm(1, phi.l[1], 1 / sqrt(phi.l[2])),
      phi.n.mu = rnorm(1, phi.n[1], 1 / sqrt(phi.n[2])),
      phi.a.mu = rnorm(1, phi.a[1], 1 / sqrt(phi.a[2])),
      theta.ln = rnorm(1, theta.l2n[1], 1 / sqrt(theta.l2n[2])),
      theta.na = rnorm(1, theta.n2a[1], 1 / sqrt(theta.n2a[2])),
      repro.mu = rnorm(1, repro[1], 1 / sqrt(repro[2])),
      beta = rnorm(n.beta, pr.beta[, 1], 1 / sqrt(pr.beta[, 2])),
      sig = rinvgamma(4, pr.sig$alpha, pr.sig$beta),
      px = matrix(rpois(4 * horizon, 5), 4, horizon),
      x = matrix(rpois(4 * horizon, 5), 4, horizon),
      y = matrix(rpois(4 * horizon, 5), 4, horizon)
    )
  }

  # source("Functions/run_hindcast_nimble.R")
  source("Scripts/DA_ticks.R")
  cl <- makeCluster(n.slots)
  out.nchains <- run_hindcast_nimble(
    cl = cl,
    model = model.code,
    data = data,
    constants = constants,
    inits = inits,
    n.iter = n.iter,
    thin = 5,
    nmme.cens = nmme.cens,
    muf.missing = muf.missing,
    use.mice = use.mice
  )
  stopCluster(cl)

  dat.hindcast <- do.call(rbind, out.nchains)

  p.monitor <- c(
    "phi.l.mu",
    "phi.n.mu",
    "phi.a.mu",
    "theta.ln",
    "theta.na",
    # "repro.mu",
    "sig[",
    "beta["
  )

  nodes <- colnames(out.nchains[[1]])
  which.tick <- grep("x", nodes)
  which.params <- which(nodes[-grep("x", nodes)] %in% nodes)
  gelman.keep <- numeric(length(nodes))
  for (ff in seq_along(nodes)) {
    mcmc.check <- list()
    col <- nodes[ff]
    for (c in 1:length(out.nchains)) {
      mcmc.check[[c]] <- coda::mcmc(out.nchains[[c]][, col])
    }
    gelman.keep[ff] <- try(coda::gelman.diag(mcmc.check, transform = TRUE)$psrf[1])
  }

  if (any(gelman.keep > 1.2)) {
    message("WARNING: Convergence not reached!")
    bad.nodes <- which(gelman.keep > 1.2)
    bad.params <- tibble(node = nodes[bad.nodes], psrf = gelman.keep[bad.nodes])
    print(head(as.matrix(bad.params)))
  } else {
    message("Convergence = TRUE")
  }

  draws <- sample.int(nrow(dat.hindcast), Nmc, replace = TRUE)
  dat.draws <- dat.hindcast[draws, ]

  # save hindcast
  save.params <- dat.draws[, which.params] %>%
    as_tibble() %>%
    mutate(
      site = as.character(ticksJob),
      paramsFrom = as.character(paramsJob)
    )

  if (use.nmme) {
    which.muf <- grep("muf[", colnames(save.params), fixed = TRUE)
    which.pf <- grep("pf[", colnames(save.params), fixed = TRUE)
    save.muf <- save.params %>% select(all_of(c(which.muf)), site, paramsFrom)
    save.pf <- save.params %>% select(all_of(c(which.pf)), site, paramsFrom)
    save.params <- save.params %>% select(-all_of(c(which.pf, which.muf)))
  }

  dates.tb <- tibble(
    fx.index = 1:length(fx.sequence),
    time = fx.sequence
  )

  ticks.long <- dat.draws[, which.tick] %>%
    as_tibble() %>%
    pivot_longer(cols = everything()) %>%
    mutate(
      lifeStage = if_else(grepl("x[1", name, fixed = TRUE), "larvae", "x"),
      lifeStage = if_else(grepl("x[2", name, fixed = TRUE), "dormant", lifeStage),
      lifeStage = if_else(grepl("x[3", name, fixed = TRUE), "nymphs", lifeStage),
      lifeStage = if_else(grepl("x[4", name, fixed = TRUE), "adults", lifeStage),
      fx.index = as.numeric(str_extract(name, "\\d*(?=\\])")),
      row = 1:n()
    )

  save.ticks <- left_join(ticks.long, dates.tb, by = "fx.index") %>%
    mutate(
      site = as.character(ticksJob),
      paramsFrom = as.character(paramsJob)
    )

  save.ticks %>%
    group_by(time, lifeStage) %>%
    summarise(
      low = quantile(value, 0.025),
      med = quantile(value, 0.5),
      high = quantile(value, 0.975)
    ) %>%
    # filter(lifeStage == "nymphs") %>%
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

  site.params <- paste0("ticksFrom_", ticksJob, "_paramsFrom_", paramsJob)
  exp <- if_else(use.nmme, paste0(experiment, "_nmme"), experiment)
  dir.write <- file.path(dir.save, site.params, modelJob, exp, as.character(fx.start.date))
  if (!dir.exists(dir.write)) dir.create(dir.write, showWarnings = FALSE, recursive = TRUE)

  write_csv(save.params, file = file.path(dir.write, "parameterSamples.csv"))
  write_csv(save.ticks, file = file.path(dir.write, "tickSamples.csv"))

  if (use.nmme) {
    write_csv(save.muf, file = file.path(dir.write, "mufSamples.csv"))
    write_csv(save.pf, file = file.path(dir.write, "pfSamples.csv"))
  }


  message("\n")
}
