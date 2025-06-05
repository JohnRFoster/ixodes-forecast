#=====================================================#
# This script creates the null model for the hindcast
# (chapter 2)
#
# The model is historical means for each week at each site
#=====================================================#

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(lubridate)
library(mgcv)

dir_top <- getwd()
dir_data <- "Data"
model <- experiment <- "null"

# set the job number to run each Grid
job_num <- 1

# first load the target data set
df_tick_calibration <- read_csv(file.path(dir_top, dir_data, "tick_cleaned"))
df_tick_hindcast <- read_csv(file.path(dir_top, dir_data, "Ticks2006_2021.csv"))
ls_vec <- c("Larvae", "Nymphs", "Adults")
grids <- c("Green Control", "Henry Control", "Tea Control")
grid_hind <- grids[job_num]
site_x <- gsub(" Control", "", grid_hind)

# filter the hindcast data to the site of interest
tick_cal <- df_tick_calibration |>
  select(Grid, DATE, n_larvae, n_nymphs, n_adults) |>
  rename(Date = DATE, Larvae = n_larvae, Nymphs = n_nymphs, Adults = n_adults)

all_ticks <- bind_rows(tick_cal, df_tick_hindcast) |>
  filter(Grid %in% grid_hind) |>
  mutate(year = year(Date), week = week(Date)) |>
  pivot_longer(
    cols = c(Larvae, Nymphs, Adults),
    names_to = "lifeStage",
    values_to = "count"
  ) |>
  arrange(Grid, Date)

hindcast_weeks <- all_ticks |>
  filter(Grid %in% grid_hind) |>
  mutate(week = week(Date)) |>
  pull(week) |>
  unique() |>
  sort()

# Forecast is weekly mean and sd by site
gam_forecast <- function(df, ls, target_date) {
  # filter out all "future observations"
  hist <- df |>
    filter(Date <= target_date[1], .data$lifeStage == ls) |>
    mutate(doy = yday(Date))

  b <- gam(count ~ s(doy), family = "poisson", data = hist, method = "REML")

  new_dat <- tibble(doy = yday(target_date))
  fit <- predict.gam(b, new_dat, type = "response", se.fit = TRUE)
  fx <- as.numeric(fit$fit)
  se_fit <- as.numeric(fit$se.fit)

  tibble(
    fx = fx,
    se_fit = se_fit
  ) |>
    mutate(
      lifeStage = ls,
      ymin = pmax(0, fx - (1.96 * se_fit)),
      ymax = fx + (1.96 * se_fit),
      var = fx,
      start.date = target_date[1],
      time = target_date,
      model = "Null",
      experiment = "Null",
      site = site_x,
      paramsFrom = site_x
    ) |>
    select(-se_fit)
}

hindcast_seq <- df_tick_hindcast |>
  filter(Grid == grid_hind) |>
  pull(Date) |>
  sort()

all_scores <- all.samples <- all_quants <- tibble()
pb <- txtProgressBar(min = 1, max = length(hindcast_seq), style = 1)
for (i in seq_along(hindcast_seq)) {
  f_date <- hindcast_seq[i]
  end_date <- f_date + 365
  fx_seq <- seq.Date(f_date, end_date, by = 1)
  fx_idx <- which(hindcast_seq %in% fx_seq)
  fx_date <- hindcast_seq[fx_idx]

  # get the forecasts we want
  fx_larvae <- gam_forecast(all_ticks, "Larvae", fx_seq)
  fx_nymphs <- gam_forecast(all_ticks, "Nymphs", fx_seq)
  fx_adults <- gam_forecast(all_ticks, "Adults", fx_seq)

  df_quants <- bind_rows(
    fx_larvae,
    fx_nymphs,
    fx_adults
  )

  tick_obs <- all_ticks |>
    filter(Date %in% fx_date) |>
    rename(time = Date, site = Grid) |>
    mutate(site = stringr::str_replace(site, " Control", ""))

  quants_with_data <- left_join(
    df_quants,
    tick_obs,
    by = c("site", "time", "lifeStage")
  ) |>
    mutate(grp = "95%") |>
    select(-year, -week)
  all_quants <- bind_rows(all_quants, quants_with_data)

  score_pois <- function(ls, fx_df) {
    y <- tick_obs |>
      filter(lifeStage == ls)

    left_join(y, fx_df, by = c("site", "time", "lifeStage")) |>
      mutate(score = scoringRules::crps_pois(count, fx)) |>
      mutate(
        lifeStage = ls,
        start.date = f_date,
        time = fx_date,
        model = model,
        experiment = experiment,
        site = site_x,
        paramsFrom = site_x,
        metric = "crps"
      )
  }

  pois_df <- bind_rows(
    score_pois("Larvae", fx_larvae),
    score_pois("Nymphs", fx_nymphs),
    score_pois("Adults", fx_adults)
  )

  all_scores <- bind_rows(all_scores, pois_df)

  setTxtProgressBar(pb, i)
}
close(pb)

site_params <- paste0("ticksFrom_", site_x, "_paramsFrom_", site_x)
dir_write <- file.path("Out", site_params, model, experiment)
if (!dir.exists(dir_write)) {
  dir.create(dir_write, showWarnings = FALSE, recursive = TRUE)
}
write_csv(all_scores, file = file.path(dir_write, "forecastScores.csv"))
write_csv(all_quants, file = file.path(dir_write, "forecastQuants.csv"))
