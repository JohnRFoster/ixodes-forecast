# script to generate figures for the manuscript

library(tidyverse)
library(lubridate)
library(MetBrewer)
library(ggpubr)
library(mgcv)
library(ggeffects)

source("R/functions_figures.R")

dir_data <- "data"
dir_out <- "out"
dir_plot <- "plots"

site_vec <- c("Green", "Henry", "Tea")
exp.vec <- c("base", "base_remove", "base_nmme_nmme", "base_remove_nmme_nmme")
model_vec <- c("Weather")

fx_scores <- read_csv(file.path(dir_data, "allForecastScores.csv"))

crps_process <- fx_scores |>
	filter(metric == "crps", model != "Null", model != "null") |>
	mutate(
		driver = if_else(grepl("nmme", experiment), "NMME", "CARY"),
		remove = if_else(grepl("remove", experiment), "No larvae", "Larvae"),
		mice = if_else(grepl("mna", experiment), "Mice", "No mice"),
		model = if_else(model == "null", "Null", model),
		model = if_else(
			model == "WithWeatherAndMiceGlobal",
			"Weather and Mice",
			model
		)
	)

# get the null forecast scores
dir_null <- "Null"
site_x <- c("Green", "Henry", "Tea")
crps_null <- tibble()
for (i in 1:3) {
	site_params <- paste0("ticksFrom_", site_x[i], "_paramsFrom_", site_x[i])
	dir_read <- file.path(dir_out, dir_null, site_params, "null", "null")
	null_score <- read_csv(file.path(dir_read, "forecastScores.csv"))
	crps_null <- bind_rows(crps_null, null_score)
}

crps_null <- crps_null |>
	mutate(
		driver = "Null",
		remove = "Null",
		mice = "Null",
		model = "Null",
		experiment = "Null",
		horizon = as.numeric(time - start.date)
	) |>
	rename(ticksFrom = site) |>
	select(all_of(colnames(crps_process)))

crps <- bind_rows(crps_process, crps_null)

cary_met <- read_csv("Data/Cary_Met_Data_Daily.csv")
met_pre_2021 <- cary_met |>
	slice(-c((nrow(cary_met) - 364):nrow(cary_met))) |>
	mutate(DATE = mdy(DATE)) |>
	filter(DATE >= "1995-01-01")
met_2021 <- cary_met |>
	slice((nrow(cary_met) - 364):nrow(cary_met)) |>
	mutate(DATE = ymd(DATE))

met <- bind_rows(met_pre_2021, met_2021) |>
	select(DATE, MAX_TEMP) |>
	mutate(
		MAX_TEMP = if_else(is.na(MAX_TEMP), mean(MAX_TEMP, na.rm = TRUE), MAX_TEMP),
		year = year(DATE)
	) |>
	mutate(gdd = pmax(0, MAX_TEMP - 10)) |>
	group_by(year) |>
	mutate(cgdd = cumsum(gdd)) |>
	ungroup() |>
	select(DATE, cgdd) |>
	rename(time = DATE)

crps_cgdd <- left_join(crps, met)

crps_period_nymph <- crps_cgdd |>
	filter(lifeStage == "Nymphs") |>
	mutate(period = if_else(cgdd >= 400 & cgdd <= 2500, "Questing", "Dormant"))

crps_period_adult <- crps_cgdd |>
	filter(lifeStage == "Adults") |>
	mutate(period = if_else(cgdd <= 1000 | cgdd >= 2500, "Questing", "Dormant"))

crps_period_larvae <- crps_cgdd |>
	filter(lifeStage == "Larvae") |>
	mutate(period = if_else(cgdd >= 1400 & cgdd <= 2500, "Questing", "Dormant"))

crps_period <- bind_rows(
	crps_period_nymph,
	crps_period_adult,
	crps_period_larvae
)

# forecast skill wrt horizon ----

crps_filter <- crps_period

crps_min <- crps_filter |>
	filter(horizon == 0) |>
	rename(score0 = score) |>
	select(
		lifeStage,
		ticksFrom,
		paramsFrom,
		start.date,
		score0,
		driver,
		remove,
		mice
	)

crps_not_min <- crps_filter |>
	select(
		lifeStage,
		horizon,
		ticksFrom,
		paramsFrom,
		start.date,
		score,
		driver,
		remove,
		mice
	)

crps_join <- left_join(crps_not_min, crps_min) |>
	filter(!is.na(score0))

crps_delta <- crps_join |>
	mutate(target = score - score0) |>
	left_join(crps_filter) |>
	mutate(year_start = year(start.date), year_end = year(time))

crps_skill <- crps_filter |>
	rename(target = score) |>
	mutate(
		year_start = year(start.date),
		year_end = year(time),
		frame = if_else(year_start == year_end, "Subannual", "Interannual")
	)

pred_data_nymphs <- make_prediction_data("Nymphs", crps_skill, crps_delta)

vline <- crps_skill |>
	filter(frame == "Interannual") |>
	pull(horizon) |>
	min()

crps_ylim <- c(-12, 20)
delta_ylim <- c(-2, 12)

plot_df <- make_plot_df(pred_data_nymphs)

## nymph lead time ----
g1_q <- wrap_plots(
	plot_df,
	pred_data_nymphs,
	"Questing",
	vline,
	crps_ylim,
	delta_ylim
)
g1_d <- wrap_plots(
	plot_df,
	pred_data_nymphs,
	"Dormant",
	vline,
	crps_ylim,
	delta_ylim
)


df_tick <- read_csv(file.path(dir_data, "Ticks2006_2021.csv"))
tick_grid <- df_tick |>
	rename(time = Date) |>
	mutate(site = gsub(" Control", "", Grid), year = year(time)) |>
	filter(site %in% site_vec) |>
	select(-Grid) |>
	pivot_longer(
		cols = c(Larvae, Nymphs, Adults),
		names_to = "lifeStage",
		values_to = "data"
	) |>
	group_by(year, site, lifeStage) |>
	mutate(peak = if_else(data == max(data), 1, 0)) |>
	ungroup() |>
	rename(ticksFrom = site) |>
	left_join(met)

g3 <- plot_ticks_on_drag_cloths(tick_grid, crps_skill, "Nymphs")

ghq <-
	ggarrange(
		g1_q,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
ghq
ggsave(
	"leadTimeGamCary_nymphs_questing.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)

ghd <-
	ggarrange(
		g1_d,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
ghd
ggsave(
	"leadTimeGamCary_nymphs_dormant.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)

## nymph forecast limit ----
plot_df |>
	mutate(frame = if_else(x < 175, "Subannual", "Interannual")) |>
	group_by(x, driver, mice, remove, period, frame) |>
	summarise(mu = mean(delta)) |>
	ungroup() |>
	filter(mu >= 0) |>
	group_by(period, driver, mice, remove, frame) |>
	filter(x == min(x)) |>
	ungroup() |>
	select(-mu) |>
	distinct() |>
	pivot_wider(values_from = x, names_from = driver) |>
	filter(period == "Questing") |>
	arrange(desc(frame), remove, mice) |>
	select(frame, remove, mice, `FX weather`, `Obs. weather`)

## larvae lead time ----
pred_data_larvae <- make_prediction_data("Larvae", crps_skill, crps_delta)
plot_df <- make_plot_df(pred_data_larvae)

crps_ylim <- c(-400, 50)
delta_ylim <- c(-275, 500)

# larvae lead time gam plots
g1_q <- wrap_plots(
	plot_df,
	pred_data_larvae,
	"Questing",
	vline,
	crps_ylim,
	delta_ylim
)
g1_d <- wrap_plots(
	plot_df,
	pred_data_larvae,
	"Dormant",
	vline,
	crps_ylim,
	delta_ylim
)

g3 <- plot_ticks_on_drag_cloths(tick_grid, crps_skill, "Larvae")

ghq <-
	ggarrange(
		g1_q,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
ghq
ggsave(
	"leadTimeGamCary_larvae_questing.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)

ghd <-
	ggarrange(
		g1_d,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
ghd
ggsave(
	"leadTimeGamCary_larvae_dormant.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)

## larvae forecast limit ----
plot_df |>
	group_by(x, driver, mice, remove, period) |>
	summarise(mu = mean(delta)) |>
	ungroup() |>
	filter(mu >= 0) |>
	group_by(period, driver, mice, remove) |>
	filter(x == min(x)) |>
	ungroup() |>
	select(-mu) |>
	distinct() |>
	pivot_wider(values_from = x, names_from = driver) |>
	arrange(period, remove, mice) |>
	select(period, remove, mice, `FX weather`, `Obs. weather`)


## adult lead time ----
pred_data_adults <- make_prediction_data("Adults", crps_skill, crps_delta)
plot_df <- make_plot_df(pred_data_adults)

crps_ylim <- c(-1.5, 6)
delta_ylim <- c(-1, 6)

# adult lead time gam plots
g1_q <- wrap_plots(
	plot_df,
	pred_data_adults,
	"Questing",
	vline,
	crps_ylim,
	delta_ylim
)
g1_d <- wrap_plots(
	plot_df,
	pred_data_adults,
	"Dormant",
	vline,
	crps_ylim,
	delta_ylim
)

g3 <- plot_ticks_on_drag_cloths(tick_grid, crps_skill, "Adults")

ghq <-
	ggarrange(
		g1_q,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
ghq
ggsave(
	"leadTimeGamCary_adults_questing.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)

ghd <-
	ggarrange(
		g1_d,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
ghd
ggsave(
	"leadTimeGamCary_adults_dormant.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)

## adult forecast limit ----
plot_df |>
	group_by(x, driver, mice, remove, period) |>
	summarise(mu = mean(delta)) |>
	ungroup() |>
	filter(mu >= 0) |>
	group_by(period, driver, mice, remove) |>
	filter(x == min(x)) |>
	ungroup() |>
	select(-mu) |>
	distinct() |>
	pivot_wider(values_from = x, names_from = driver) |>
	arrange(period, remove, mice) |>
	select(period, remove, mice, `FX weather`, `Obs. weather`)

# Forecast improvements ----

params <- read_csv(file.path("analysis", "allParameterQuants.csv"))

params_mutate <- params |>
	parameter_names()

pf <- "Henry"

site_cols <- c(
	"Green" = "#66c2a5",
	"Henry" = "#fc8d62",
	"Tea" = "#8da0cb"
)

my_theme <- theme(
	axis.text = element_text(size = 8),
	axis.title = element_text(size = 8),
	strip.text = element_text(size = 6),
	legend.text = element_text(size = 8),
	legend.title = element_text(size = 8)
)

g1 <- params_mutate |>
	filter(
		lifeStage == "Nymphs",
		paramsFrom == pf,
		mice == "Mice",
		remove == "Larvae",
		parameter.name == "Weather Effect",
		weather == "CARY"
	) |>
	ggplot() +
	aes(x = start.date, fill = site) +
	geom_ribbon(aes(ymin = q0.025, ymax = q0.975), alpha = 0.5) +
	geom_hline(aes(yintercept = 0), linetype = "dotted") +
	scale_fill_manual(values = site_cols) +
	scale_alpha_manual(values = 0.5) +
	facet_wrap(~driver, scales = "free") +
	labs(
		x = "Forecast start date",
		y = "Effect",
		linetype = "Model",
		fill = "Site"
	) +
	theme_bw() +
	theme(
		legend.position = "right",
		axis.text.x = element_text(angle = 90, vjust = 0.5)
	) +
	my_theme

m <- "Mice"
r <- "Larvae"
p <- "Questing"
d <- "CARY"
ls <- "Nymphs"

df_process <- crps_period |>
	filter(
		#paramsFrom == pf,
		mice == m,
		remove == r,
		period == p,
		driver == d
	)

start_dates_plot <- unique(df_process$start.date)

df_null <- crps_period |>
	filter(
		driver == "Null",
		period == p,
		start.date %in% start_dates_plot
	)

df_null |>
	mutate(year = year(start.date)) |>
	group_by(ticksFrom, lifeStage) |>
	summarise(m = lm(score ~ year)$coefficients[2]) |>
	ungroup() |>
	mutate(m = round(m, 2)) |>
	pivot_wider(names_from = ticksFrom, values_from = m)

# average rate of CRPS change per year
crps_period_filter <- crps_period |>
	filter(
		mice == m,
		remove == r,
		period == p,
		driver == d
	) |>
	mutate(year = year(start.date))

crps_period_filter |>
	group_by(lifeStage, ticksFrom) |>
	summarise(m = lm(score ~ year)$coefficients[2]) |>
	ungroup() |>
	mutate(m = round(m, 2)) |>
	pivot_wider(names_from = ticksFrom, values_from = m)

get_predicted <- function(tf) {
	tmp1 <- crps_period_filter |> filter(ticksFrom == tf, lifeStage == "Nymphs")
	fit <- lm(score ~ year, data = tmp1)
	p1 <- predict(fit, newdata = data.frame(year = 2006))
	p2 <- predict(fit, newdata = data.frame(year = 2020))
	round((p1 - p2) / p1 * 100)
}

# linear improvement in skill first year vs last year
get_predicted("Green")
get_predicted("Henry")
get_predicted("Tea")

smooth <- "lm"

g2 <- crps_period |>
	filter(
		mice == m,
		remove == r,
		period == p,
		driver == d
	) |>
	ggplot() +
	aes(x = start.date, y = score) +
	geom_smooth(
		aes(color = ticksFrom, linetype = "Process"),
		method = smooth,
		se = FALSE
	) +
	geom_smooth(
		data = df_null,
		aes(color = ticksFrom, linetype = "Null"),
		method = smooth,
		se = FALSE
	) +
	scale_color_manual(values = site_cols) +
	scale_linetype_manual(values = c("Process" = "solid", "Null" = "dotted")) +
	labs(
		x = "Forecast start date",
		y = "CRPS",
		color = "Site",
		linetype = "Model"
	) +
	theme_bw() +
	theme(
		legend.position = "right",
		axis.text.x = element_text(angle = 90, vjust = 0.5)
	) +
	my_theme

g2

ggarrange(g1, g2, nrow = 1, legend = "right", labels = "AUTO")

ggsave(
	"skillClimateRelationship.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 8,
	units = "cm"
)

## phenology ----

draw <- seq(1, nrow(crps_period), length.out = 5)
labels <- month(sort(month(crps_period$time))[draw], label = TRUE)
breaks <- sort(yday(crps_period$time))[draw]

xmin <- crps_period |>
	filter(period == "Questing") |>
	mutate(doy = yday(time), year = year(time)) |>
	group_by(year) |>
	filter(doy == min(doy)) |>
	pull(doy) |>
	mean()
xmax <- crps_period |>
	filter(period == "Questing") |>
	mutate(doy = yday(time), year = year(time)) |>
	group_by(year) |>
	filter(doy == max(doy)) |>
	pull(doy) |>
	mean()

phenology_process <- crps_period |>
	filter(driver != "Null") |>
	mutate(doy = yday(time))
phenology_null <- crps_period |>
	filter(driver == "Null") |>
	mutate(doy = yday(time), instance = "Day-of-year")

g1 <- plot_phenology_scores("CARY", "Nymphs", phenology_process, phenology_null)
g2 <- plot_phenology_scores("NMME", "Nymphs", phenology_process, phenology_null)

ggarrange(
	g1,
	g2,
	align = "hv",
	common.legend = TRUE,
	legend = "bottom",
	labels = "AUTO"
)

ggsave(
	"phenologyScores.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 8,
	units = "cm"
)

## time series ----

dir_top <- "C:/Users/John.Foster/Downloads/nefi"

all.fx.files <- list.files(dir_top, recursive = TRUE)
state_files <- grep("tickSamples", all.fx.files, value = TRUE)
state_files <- grep(
	"ticksFrom_Green_paramsFrom_Green",
	state_files,
	value = TRUE
)

group <- str_extract(state_files, "(?<=\\/)[[:graph:]]*(?=/2)")
model <- str_extract(group, "[[:alpha:]]*(?=\\/)")
experiment <- str_extract(group, "(?<=\\/)[[:graph:]]*")
experiment <- gsub("nmme_nmme", "nmme", experiment)

quants <- tibble()
for (i in seq_along(state_files)) {
	df <- read_csv(file.path(dir_top, state_files[i])) |> suppressMessages()
	start_date <- min(df$time)
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
			horizon = as.numeric(time - start_date),
			ticksFrom = t_from,
			paramsFrom = p_from,
			experiment = experiment,
			model = model,
			start_date = start_date
		) |>
		mutate(
			lifeStage = if_else(lifeStage == "larvae", "Larvae", lifeStage),
			lifeStage = if_else(lifeStage == "nymphs", "Nymphs", lifeStage),
			lifeStage = if_else(lifeStage == "adults", "Adults", lifeStage)
		) |>
		select(-site) |>
		group_by(
			lifeStage,
			time,
			paramsFrom,
			ticksFrom,
			experiment,
			model,
			start_date
		) |>
		summarise(
			fx = median(value),
			ymin = quantile(value, 0.025),
			ymax = quantile(value, 0.975)
		)

	quants <- bind_rows(df_mutate, quants)
}

fname <- file.path(
	dir_out,
	dir_null,
	"ticksFrom_Green_paramsFrom_Green/null/null/forecastQuants.csv"
)
null_fx3 <- read_csv(fname)

start_dates_vec <- sort(unique(quants$start_date), decreasing = FALSE)
start_dates_vec <- start_dates_vec[grep("2018", start_dates_vec)]
g <- list()
for (i in seq_along(start_dates_vec)) {
	start_date_filter <- start_dates_vec[i]

	dfq <- quants |>
		filter(start_date == start_date_filter, lifeStage == "Nymphs") |>
		ungroup() |>
		mutate(
			driver = if_else(
				grepl("nmme", experiment),
				"Forecasted weather",
				"Observed weather"
			),
			remove = if_else(grepl("remove", experiment), "No larvae", "Larvae"),
			mice = if_else(grepl("mna", experiment), "Mice", "No mice")
		) |>
		select(ticksFrom, fx, ymin, ymax, time, lifeStage, remove, mice, driver)

	null <- null_fx3 |>
		filter(
			start.date == start_date_filter,
			site == dfq$ticksFrom[1],
			lifeStage == "Nymphs"
		) |>
		select(fx, ymin, ymax, time, count, lifeStage) |>
		mutate(remove = "Null", mice = "Null", driver = "Null")

	obs <- null |>
		filter(lifeStage == "Nymphs") |>
		select(time, count, lifeStage)

	g[[i]] <- bind_rows(null, dfq) |>
		left_join(obs) |>
		filter(
			lifeStage == "Nymphs",
			remove == "No larvae" | remove == "Null",
			mice == "Mice" | mice == "Null"
		) |>
		mutate(
			driver = factor(
				driver,
				levels = c("Null", "Observed weather", "Forecasted weather")
			)
		) |>
		ggplot() +
		aes(
			x = time,
			y = fx,
			ymin = ymin,
			ymax = ymax,
			fill = driver,
			linetype = driver
		) +
		geom_ribbon(alpha = 0.25) +
		scale_fill_brewer(type = "qual", palette = "Dark2") +
		geom_line() +
		geom_point(aes(y = count, shape = "Nymphs on drag cloths")) +
		coord_cartesian(ylim = c(0, 35)) +
		labs(
			title = paste0("Forecast issued on: ", start_date),
			linetype = "Model type",
			fill = "Model type",
			x = "Date",
			y = "Ticks/450 sq. m",
			shape = ""
		) +
		theme_bw()
}

length(g)
ggarrange(
	g[[1]],
	g[[2]],
	g[[3]],
	g[[4]],
	g[[5]],
	nrow = 3,
	ncol = 2,
	align = "hv",
	common.legend = TRUE,
	legend = "bottom"
)

ggsave(
	"timeSeries2018.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 24,
	units = "cm"
)

## transferability ----

w <- 0.5

transfer_cols <- c(
	"Within-site" = "#f1a340",
	"Across-site" = "#998ec3"
)

crps_instance_p <- crps_period |>
	filter(driver != "Null") |>
	mutate(
		frame = if_else(year(time) == year(start_date), "Subannual", "Interannual"),
		driver = if_else(
			driver == "CARY",
			"Observed weather",
			"Forecasted weather"
		),
		instance = paste(remove, mice, sep = " "),
		transfer = if_else(ticksFrom == paramsFrom, "Within-site", "Across-site")
	)

n1 <- crps_period |>
	filter(driver == "Null") |>
	mutate(
		instance = "Day-of-year",
		driver = "Observed weather",
		transfer = "Within-site",
		frame = if_else(year(time) == year(start_date), "Subannual", "Interannual")
	)
n2 <- crps_period |>
	filter(driver == "Null") |>
	mutate(
		instance = "Day-of-year",
		driver = "Forecasted weather",
		transfer = "Within-site",
		frame = if_else(year(time) == year(start_date), "Subannual", "Interannual")
	)

n3 <- bind_rows(n1, n2)

transfer_cis <- bind_rows(crps_instance_p, n3) |>
	group_by(period, frame, driver, instance, transfer, lifeStage) |>
	summarise(
		med = DescTools::MedianCI(score)["median"],
		lwr.ci = DescTools::MedianCI(score)["lwr.ci"],
		upr.ci = DescTools::MedianCI(score)["upr.ci"]
	) |>
	filter(period == "Questing") |>
	mutate(
		instance = factor(
			instance,
			levels = c(
				"No larvae No mice",
				"Larvae No mice",
				"No larvae Mice",
				"Larvae Mice",
				"Day-of-year"
			)
		)
	) |>
	ungroup()

transfer_plot(transfer_cis, "Nymphs")
ggsave(
	"transferability_nymphs.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 10,
	units = "cm"
)

transfer_plot(transfer_cis, "Larvae")
ggsave(
	"transferability_larvae.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 10,
	units = "cm"
)

transfer_plot(transfer_cis, "Adults")
ggsave(
	"transferability_adults.jpeg",
	dpi = "retina",
	path = dir_plot,
	width = 18,
	height = 10,
	units = "cm"
)
