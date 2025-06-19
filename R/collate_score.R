# Collate, score, and calculate quantiles for process model forecasts

library(tidyverse)
library(lubridate)

dir_data <- "data"
dir_out <- "out"

ls_vec <- c("Larvae", "Nymphs", "Adults")

## calculate scores -----------------------------------------------------------------------------------
calc_scores <- function(df_fx) {
	all.scores <- tibble()
	for (s in seq_along(ls_vec)) {
		df_filter <- df_fx |>
			filter(lifeStage == ls_vec[s]) |>
			arrange(time)

		y <- df_filter |>
			select(time, data) |>
			distinct() |>
			pull(data)

		dat <- df_filter |>
			select(fx, time) |>
			group_by(time) |>
			mutate(mcmc = seq_len(n())) |>
			pivot_wider(names_from = time, values_from = fx) |>
			select(-mcmc) |>
			as.matrix() |>
			t()

		crps <- scoringutils::crps_sample(y, dat)
		logs <- scoringutils::logs_sample(y, dat)
		mean.sq.error <- scoringutils::se_mean_sample(y, dat)
		interval <- scoringutils::interval_score(
			y,
			apply(dat, 1, quantile, 0.025),
			apply(dat, 1, quantile, 0.975),
			95
		)

		df_score <- df_filter |>
			select(-fx, -data) |>
			distinct() |>
			mutate(
				crps = crps,
				logs = logs,
				mean.sq.error = mean.sq.error,
				interval = interval
			) |>
			pivot_longer(
				cols = c(crps, logs, mean.sq.error, interval),
				names_to = "metric",
				values_to = "score"
			)

		all.scores <- bind_rows(all.scores, df_score)
	}
	all.scores
}

# get the sequence if files we want
all.fx.files <- list.files(dir_out, recursive = TRUE)
all.fx.files <- all.fx.files[-grep("analysis", all.fx.files)]
state_files <- grep("tickSamples", all.fx.files, value = TRUE)

df_tick <- read_csv(file.path(dir_data, "Ticks2006_2021.csv"))
data_tick <- df_tick |>
	filter(Grid %in% paste(c("Green", "Henry", "Tea"), "Control")) |>
	mutate(Larvae = Larvae + 1) |>
	pivot_longer(
		cols = c(Larvae, Nymphs, Adults),
		names_to = "lifeStage",
		values_to = "data"
	) |>
	mutate(ticksFrom = gsub(" Control", "", Grid)) |>
	select(-Grid) |>
	rename(time = Date)

group <- str_extract(state_files, "(?<=\\/)[[:graph:]]*(?=/2)")
model <- str_extract(group, "[[:alpha:]]*(?=\\/)")
experiment <- str_extract(group, "(?<=\\/)[[:graph:]]*")
experiment <- gsub("nmme_nmme", "nmme", experiment)

fx_quants_data_long <- tibble()
fx_quants_all_long <- tibble()
fx_samples_long <- tibble()
scores_long <- tibble()

pb <- txtProgressBar(min = 1, max = length(state_files), style = 1)

message("----------------- state files  ------------------------")

for (i in seq_along(state_files)) {
	df <- read_csv(file.path(dir_out, state_files[i])) |> suppressMessages()
	start.date <- min(df$time) # nolint
	group <- str_extract(state_files[i], "(?<=\\/)[[:graph:]]*(?=/2)")
	model <- str_extract(group, "[[:alpha:]]*(?=\\/)")
	# if(is.na(model)) next
	experiment <- str_extract(group, "(?<=\\/)[[:graph:]]*")
	experiment <- gsub("nmme_nmme", "nmme", experiment)
	p_from <- str_extract(state_files[i], "(?<=paramsFrom_)[[:alpha:]]*(?=\\/)")
	t_from <- str_extract(state_files[i], "(?<=ticksFrom_)[[:alpha:]]*(?=\\_)")

	df_mutate <- df |>
		filter(lifeStage != "dormant") |>
		mutate(
			horizon = as.numeric(time - start.date),
			ticksFrom = t_from,
			paramsFrom = p_from,
			experiment = experiment,
			model = model,
			start.date = start.date
		) |>
		mutate(
			lifeStage = if_else(lifeStage == "larvae", "Larvae", lifeStage),
			lifeStage = if_else(lifeStage == "nymphs", "Nymphs", lifeStage),
			lifeStage = if_else(lifeStage == "adults", "Adults", lifeStage)
		) |>
		select(-site)

	if (model != "null") {
		df_mutate <- df_mutate |>
			rename(fx = value) |>
			select(-row, -name, -fx.index)
	} else {
		df_mutate <- df_mutate |>
			mutate(horizon = NA)
	}

	df_quants <- df_mutate |>
		group_by(time, lifeStage, start.date) |>
		summarise(
			q0.025 = quantile(fx, 0.025),
			q0.25 = quantile(fx, 0.125),
			q0.5 = quantile(fx, 0.5),
			q0.75 = quantile(fx, 0.875),
			q0.975 = quantile(fx, 0.975),
			var = var(fx)
		) |>
		mutate(
			horizon = as.numeric(time - start.date),
			ticksFrom = t_from,
			paramsFrom = p_from,
			experiment = experiment,
			model = model,
			start.date = start.date
		) |>
		suppressMessages()

	quants_with_data <- left_join(
		df_quants,
		data_tick,
		by = c("time", "lifeStage", "ticksFrom")
	) |>
		pivot_longer(
			cols = c(q0.025, q0.25),
			values_to = "ymin",
			names_to = "lowQ"
		) |>
		pivot_longer(
			cols = c(q0.75, q0.975),
			values_to = "ymax",
			names_to = "highQ"
		) |>
		filter(
			lowQ == "q0.025" & highQ == "q0.975" | lowQ == "q0.25" & highQ == "q0.75"
		) |>
		mutate(grp = if_else(highQ == "q0.975", "95%", "75%")) |>
		select(-lowQ, -highQ)

	# plot quants
	for (s in 1:3) {
		plot_single_forecast(
			quants_with_data,
			s,
			file.path(dir.plot, model, experiment)
		)
	}

	# fx_quants_all_long <- bind_rows(fx_quants_all_long, quants_with_data)

	quants_data_only <- quants_with_data |>
		filter(!is.na(data))
	fx_quants_data_long <- bind_rows(fx_quants_data_long, quants_data_only)

	samples_data_only <- left_join(
		df_mutate,
		data_tick,
		by = c("time", "lifeStage", "ticksFrom")
	) |>
		filter(!is.na(data))
	fx_samples_long <- bind_rows(fx_samples_long, samples_data_only)

	# calculate scores
	scores <- calc_scores(samples_data_only)
	scores_long <- bind_rows(scores_long, scores)

	setTxtProgressBar(pb, i)
}
close(pb)

dir_write <- "analysis"
if (!dir.exists(dir_write)) {
	dir.create(dir_write, recursive = TRUE, showWarnings = FALSE)
}

write_csv(
	fx_quants_data_long,
	file = file.path(dir_write, "allForecastQuantiles.csv")
)
# write_csv(fx_quants_all_long, file = file.path(dir_write, "forecastQuantilesAllTime.csv"))
write_csv(
	fx_samples_long,
	file = file.path(dir_write, "allForecastSamples.csv")
)
write_csv(scores_long, file = file.path(dir_write, "allForecastScores.csv"))


# parameters ------------------------------------------------------------------------------------------------

message("----------------- parameter files  ------------------------")

param_files <- grep("parameterSamples", all.fx.files, value = TRUE)
pb <- txtProgressBar(min = 1, max = length(param_files), style = 1)

all.params <- all.quants <- tibble()
for (i in seq_along(param_files)) {
	df <- read_csv(file.path(dir_out, param_files[i])) |> suppressMessages()
	start.date <- str_extract(param_files[i], "\\d{4}-\\d{2}-\\d{2}") # nolint
	group <- str_extract(param_files[i], "(?<=\\/)[[:graph:]]*(?=/2)")
	model <- str_extract(group, "[[:alpha:]]*(?=\\/)")
	experiment <- str_extract(group, "(?<=\\/)[[:graph:]]*")
	experiment <- gsub("nmme_nmme", "nmme", experiment)
	p_from <- str_extract(param_files[i], "(?<=paramsFrom_)[[:alpha:]]*(?=\\/)")
	t_from <- str_extract(param_files[i], "(?<=ticksFrom_)[[:alpha:]]*(?=\\_)")

	df_wide <- df |>
		mutate(
			ticksFrom = t_from,
			paramsFrom = p_from,
			experiment = experiment,
			model = model,
			start.date = start.date
		) |>
		select(-site)

	all.params <- bind_rows(all.params, df_wide)

	params_quantiles <- df_wide |>
		pivot_longer(cols = where(is.numeric), names_to = "parameter") |>
		group_by(parameter) |>
		summarise(
			q0.025 = quantile(value, 0.025),
			q0.05 = quantile(value, 0.05),
			q0.5 = quantile(value, 0.5),
			q0.95 = quantile(value, 0.95),
			q0.975 = quantile(value, 0.975)
		) |>
		mutate(
			ticksFrom = t_from,
			paramsFrom = p_from,
			experiment = experiment,
			model = model,
			start.date = start.date
		)

	all.quants <- bind_rows(all.quants, params_quantiles)

	setTxtProgressBar(pb, i)
}
close(pb)


write_csv(all.params, file = file.path(dir_write, "allParameterSamples.csv"))
write_csv(all.quants, file = file.path(dir_write, "allParameterQuants.csv"))
