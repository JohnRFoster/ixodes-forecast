library(tidyverse)
library(lubridate)
library(MetBrewer)
library(ggpubr)
library(mgcv)
library(ggeffects)

source("Scripts/hindcastManuscriptFigureFunctions.R")

# dir_top <- "/projectnb/dietzelab/fosterj"
dir_data <- "Data"
dir_analysis <- "Analysis"
dir_plot <- "Plots"
files <- list.files(dir_analysis, recursive = TRUE)

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


site_x <- c("Green", "Henry", "Tea")
crps_null <- tibble()
for (i in 1:3) {
	site_params <- paste0("ticksFrom_", site_x[i], "_paramsFrom_", site_x[i])
	dir_read <- file.path(dir_data, site_params, "null", "null")
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
		horizon = as.numeric(time - start_date)
	) |>
	rename(ticksFrom = site) |>
	select(all_of(colnames(crps_process)))

crps <- bind_rows(crps_process, crps_null)

cary_met <- read_csv("Data/Cary_Met_Data_Daily.csv")
met_pre_2021 <- cary_met %>%
	slice(-c((nrow(cary_met) - 364):nrow(cary_met))) %>%
	mutate(DATE = mdy(DATE)) %>%
	filter(DATE >= "1995-01-01")
met_2021 <- cary_met %>%
	slice((nrow(cary_met) - 364):nrow(cary_met)) %>%
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

# horizon ----

crps_filter <- crps_period

crps_min <- crps_filter |>
	filter(horizon == 0) |>
	rename(score0 = score) |>
	select(
		lifeStage,
		ticksFrom,
		paramsFrom,
		start_date,
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
		start_date,
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
	mutate(year_start = year(start_date), year_end = year(time))

crps_skill <- crps_filter |>
	rename(target = score) |>
	mutate(
		year_start = year(start_date),
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
g <- list()
g[[1]] <- plot_df |>
	plot_gam("score", "Obs. weather", vline) +
	coord_cartesian(ylim = crps_ylim)
g[[2]] <- plot_df |>
	plot_gam("score", "FX weather", vline) +
	coord_cartesian(ylim = crps_ylim)
g[[3]] <- pred_data_nymphs |>
	filter(instance != "Day-of-year") |>
	rename(delta = value) |>
	plot_gam("delta", "Obs. weather", vline) +
	coord_cartesian(ylim = delta_ylim) +
	geom_hline(yintercept = 0, linetype = "dashed")
g[[4]] <- pred_data_nymphs |>
	filter(instance != "Day-of-year") |>
	rename(delta = value) |>
	plot_gam("delta", "FX weather", vline) +
	coord_cartesian(ylim = delta_ylim) +
	geom_hline(yintercept = 0, linetype = "dashed")

g1 <- ggarrange(
	plotlist = g,
	nrow = 2,
	ncol = 2,
	align = "hv",
	legend = "bottom",
	labels = "AUTO",
	common.legend = TRUE
)

g1

df_tick <- read_csv(file.path(dir_data, "Ticks2006_2021.csv"))
tick_grid <- df_tick %>%
	rename(time = Date) %>%
	mutate(site = gsub(" Control", "", Grid), year = year(time)) %>%
	filter(site %in% site_vec) %>%
	select(-Grid) %>%
	pivot_longer(
		cols = c(Larvae, Nymphs, Adults),
		names_to = "lifeStage",
		values_to = "data"
	) %>%
	group_by(year, site, lifeStage) %>%
	mutate(peak = if_else(data == max(data), 1, 0)) %>%
	ungroup() |>
	rename(ticksFrom = site) |>
	left_join(met)

g3 <- plot_ticks_on_drag_cloths(tick_grid, crps_skill, "Nymphs")

gh <-
	ggarrange(
		g1,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
gh
ggsave(
	"leadTimeGamCary_nymphs.jpeg",
	dpi = "retina",
	path = "Plots/hindcast_manuscript",
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

# nymph lead time gam plots
g <- list()
g[[1]] <- plot_df |>
	plot_gam("score", "Obs. weather", vline) +
	coord_cartesian(ylim = crps_ylim)
g[[2]] <- plot_df |>
	plot_gam("score", "FX weather", vline) +
	coord_cartesian(ylim = crps_ylim)
g[[3]] <- pred_data_larvae |>
	filter(instance != "Day-of-year") |>
	rename(delta = value) |>
	plot_gam("delta", "Obs. weather", vline) +
	coord_cartesian(ylim = delta_ylim) +
	geom_hline(yintercept = 0, linetype = "dashed")
g[[4]] <- pred_data_larvae |>
	filter(instance != "Day-of-year") |>
	rename(delta = value) |>
	plot_gam("delta", "FX weather", vline) +
	coord_cartesian(ylim = delta_ylim) +
	geom_hline(yintercept = 0, linetype = "dashed")

g1 <- ggarrange(
	plotlist = g,
	nrow = 2,
	ncol = 2,
	align = "hv",
	legend = "bottom",
	labels = "AUTO",
	common.legend = TRUE
)

g1

g3 <- plot_ticks_on_drag_cloths(tick_grid, crps_skill, "Larvae")

gh <-
	ggarrange(
		g1,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
gh
ggsave(
	"leadTimeGamCary_larvae.jpeg",
	dpi = "retina",
	path = "Plots/hindcast_manuscript",
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


# adult gam plots
pred_data_adults <- make_prediction_data("Adults", crps_skill, crps_delta)
plot_df <- make_plot_df(pred_data_adults)

crps_ylim <- c(-1.5, 6)
delta_ylim <- c(-1, 6)

# nymph lead time gam plots
g <- list()
g[[1]] <- plot_df |>
	plot_gam("score", "Obs. weather", vline) +
	coord_cartesian(ylim = crps_ylim)
g[[2]] <- plot_df |>
	plot_gam("score", "FX weather", vline) +
	coord_cartesian(ylim = crps_ylim)
g[[3]] <- pred_data_adults |>
	filter(instance != "Day-of-year") |>
	rename(delta = value) |>
	plot_gam("delta", "Obs. weather", vline) +
	coord_cartesian(ylim = delta_ylim) +
	geom_hline(yintercept = 0, linetype = "dashed")
g[[4]] <- pred_data_adults |>
	filter(instance != "Day-of-year") |>
	rename(delta = value) |>
	plot_gam("delta", "FX weather", vline) +
	coord_cartesian(ylim = delta_ylim) +
	geom_hline(yintercept = 0, linetype = "dashed")

g1 <- ggarrange(
	plotlist = g,
	nrow = 2,
	ncol = 2,
	align = "hv",
	legend = "bottom",
	labels = "AUTO",
	common.legend = TRUE
)

g1

g3 <- plot_ticks_on_drag_cloths(tick_grid, crps_skill, "Adults")

gh <-
	ggarrange(
		g1,
		g3,
		nrow = 2,
		heights = c(2, 1),
		labels = c("", "E"),
		legend = "bottom"
	)
gh
ggsave(
	"leadTimeGamCary_adults.jpeg",
	dpi = "retina",
	path = "Plots/hindcast_manuscript",
	width = 18,
	height = 18,
	units = "cm"
)


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

#### forecast limit

# Forecast improvements ----

params <- read_csv("data/allParameterQuants.csv")

params_mutate <- params %>%
	mutate(
		weather = if_else(grepl("nmme", experiment), "NMME", "CARY"),
		remove = if_else(grepl("remove", experiment), "No larvae", "Larvae"),
		mice = if_else(grepl("mna", experiment), "Mice", "No Mice"),
		within = if_else(paramsFrom == ticksFrom, "Within", "Across"),
		parameter.name = if_else(grepl("beta", parameter), "Weather Effect", "x"),
		parameter.name = if_else(
			grepl("phi.a.mu", parameter),
			"Adult Survival Rate",
			parameter.name
		),
		parameter.name = if_else(
			grepl("phi.l.mu", parameter),
			"Larvae Survival Rate",
			parameter.name
		),
		parameter.name = if_else(
			grepl("phi.n.mu", parameter),
			"Nymph Survival Rate",
			parameter.name
		),
		parameter.name = if_else(
			parameter == "sig[1]",
			"Larvae process variance",
			parameter.name
		),
		parameter.name = if_else(
			parameter == "sig[2]",
			"Dormant process variance",
			parameter.name
		),
		parameter.name = if_else(
			parameter == "sig[3]",
			"Nymph process variance",
			parameter.name
		),
		parameter.name = if_else(
			parameter == "sig[4]",
			"Adult process variance",
			parameter.name
		),
		parameter.name = if_else(
			grepl("theta.ln", parameter),
			"Larvae-Nymph Transition Rate",
			parameter.name
		),
		parameter.name = if_else(
			grepl("theta.na", parameter),
			"Nymph-Adult Transition Rate",
			parameter.name
		),
		parameter.name = if_else(
			parameter == "beta[13]",
			"Larvae-Nymph Effect",
			parameter.name
		),
		parameter.name = if_else(
			parameter == "beta[14]",
			"Nymph-Adult Effect",
			parameter.name
		),
		driver = if_else(
			parameter %in% c("beta[1]", "beta[5]", "beta[9]"),
			"Max Temp",
			"x"
		),
		driver = if_else(
			parameter %in% c("beta[2]", "beta[6]", "beta[10]"),
			"Max RH",
			driver
		),
		driver = if_else(
			parameter %in% c("beta[3]", "beta[7]", "beta[11]"),
			"Min RH",
			driver
		),
		driver = if_else(
			parameter %in% c("beta[4]", "beta[8]", "beta[12]"),
			"Precip",
			driver
		),
		driver = if_else(parameter %in% c("beta[13]", "beta[14]"), "Mice", driver),
		lifeStage = if_else(
			parameter %in% c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),
			"Larvae",
			"transition"
		),
		lifeStage = if_else(
			parameter %in% c("beta[5]", "beta[6]", "beta[7]", "beta[8]"),
			"Nymphs",
			lifeStage
		),
		lifeStage = if_else(
			parameter %in% c("beta[9]", "beta[10]", "beta[11]", "beta[12]"),
			"Adults",
			lifeStage
		),
		lifeStage = if_else(grepl("phi.l.mu", parameter), "Larvae", lifeStage),
		lifeStage = if_else(grepl("phi.n.mu", parameter), "Nymphs", lifeStage),
		lifeStage = if_else(grepl("phi.a.mu", parameter), "Adults", lifeStage)
	) %>%
	rename(site = ticksFrom)

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
	aes(x = start_date, fill = site) +
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

df_process <- crps_period |>
	filter(
		paramsFrom == pf,
		mice == m,
		remove == r,
		period == p,
		driver == d
	)

start_dates_plot <- unique(df_process$start_date)

df_null <- crps_period |>
	filter(driver == "Null", period == p, start_date %in% start_dates_plot)

df_null |>
	mutate(year = year(start_date)) |>
	group_by(ticksFrom) |>
	summarise(m = lm(score ~ year)$coefficients[2]) |>
	ungroup() |>
	mutate(m = round(m, 2)) |>
	pivot_wider(names_from = ticksFrom, values_from = m)

# average rate of CRPS change per year
crps_period |>
	filter(
		mice == m,
		remove == r,
		period == p,
		driver == d
	) |>
	mutate(year = year(start_date)) |>
	group_by(ticksFrom) |>
	summarise(m = lm(score ~ year)$coefficients[2]) |>
	ungroup() |>
	mutate(m = round(m, 2)) |>
	pivot_wider(names_from = ticksFrom, values_from = m)

smooth <- "lm"

g2 <- crps_period |>
	filter(
		mice == m,
		remove == r,
		period == p,
		driver == d
	) |>
	ggplot() +
	aes(x = start_date, y = score) +
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
	path = "Plots/hindcast/manuscript",
	width = 18,
	height = 8,
	units = "cm"
)


## interannual variability

iav_scores <- crps_period |>
	mutate(year = year(time)) |>
	filter(period == "Questing", driver == "CARY" | driver == "Null") |>
	mutate(instance = paste(driver, remove, mice, sep = " ")) |>
	mutate(
		instance = stringr::str_replace(instance, "CARY ", ""),
		instance = if_else(grepl("Null", instance), "Day-of-year", instance)
	) |>
	group_by(year, instance) |>
	summarise(
		mu = median(score)
	) |>
	ungroup()

g3 <- iav_scores |>
	ggplot() +
	aes(x = year, y = mu, color = instance) +
	geom_line() +
	# geom_smooth(method = "lm", se = FALSE) +
	geom_point() +
	labs(x = "Year", y = "Median CRPS", color = "Sampling strategy") +
	scale_color_manual(values = instance_cols) +
	theme_bw() +
	theme(axis.title.x = element_blank())
g3

iav <- tick_grid |>
	left_join(met) |>
	filter(lifeStage == "Nymphs", period == "Questing") |>
	distinct() |>
	group_by(year) |>
	summarise(mu = mean(data)) |>
	ungroup() |>
	mutate(iav = sd(mu))

g4 <- iav |>
	ggplot() +
	aes(x = year, y = mu) +
	geom_point() +
	geom_line() +
	geom_hline(aes(yintercept = iav, linetype = "Interannual\nvariability")) +
	scale_linetype_manual(values = "dashed") +
	labs(x = "Year", y = "Mean ticks / 450 km2", linetype = "") +
	theme_bw() +
	theme(strip.text = element_blank(), strip.background = element_blank())

g4

ggarrange(
	g3,
	g4,
	nrow = 2,
	ncol = 1,
	legend = "right",
	labels = "AUTO",
	align = "v"
)

ggsave(
	"interannualVariability.jpeg",
	dpi = "retina",
	path = "Plots/hindcast/manuscript",
	width = 18,
	height = 12,
	units = "cm"
)

iav_scores_null <- iav_scores |>
	filter(instance == "Null") |>
	rename(null_score = mu) |>
	select(-instance)
iav_scores_proc <- iav_scores |>
	filter(instance != "Null") |>
	rename(proc_score = mu)

score_delta <- left_join(iav_scores_proc, iav_scores_null) |>
	mutate(delta = null_score - proc_score)

iav_with_scores <- left_join(score_delta, iav)

iav_with_scores |> filter(mu > iav) |> pull(year) |> unique() |> length()

iav_with_scores |>
	filter(mu > iav, delta < 0)

iav_with_scores |>
	filter(mu <= iav, delta < 0)

left_join(score_delta, iav) |>
	filter(year == 2008)

## phenology ----

draw <- seq(1, nrow(crps_period), length.out = 5)
labels <- month(crps_period$time) %>%
	sort() %>%
	.[draw] %>%
	month(., label = TRUE)
breaks <- yday(crps_period$time) %>% sort() %>% .[draw]

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

phenology_process <- crps_period %>%
	filter(driver != "Null") |>
	mutate(doy = yday(time))
phenology_null <- crps_period %>%
	filter(driver == "Null") |>
	mutate(doy = yday(time), instance = "Day-of-year")

rect <- tibble(
	lifeStage = c("Larvae", "Nymphs", "Adults"),
	xmin = c(yday("2022-07-01"), yday("2022-05-01"), yday("2022-08-01")),
	xmax = c(yday("2022-10-01"), yday("2022-08-01"), yday("2022-05-01"))
)

plot_phenology_scores <- function(d) {
	gg <- phenology_process |>
		filter(driver == d) |>
		mutate(
			driver = if_else(
				driver == "CARY",
				"Observed weather",
				"Forecasted weather"
			)
		) |>
		mutate(instance = paste(remove, mice, sep = " ")) |>
		mutate(
			instance = stringr::str_replace(instance, "CARY ", ""),
			instance = if_else(grepl("Null", instance), "Day-of-year", instance)
		) |>
		bind_rows(phenology_null) |>
		ggplot() +
		aes(x = doy) +
		geom_smooth(
			aes(y = score, color = instance),
			se = FALSE,
			method = "loess"
		) +
		scale_color_manual(values = instance_cols) +
		scale_y_continuous(limits = c(0, NA)) +
		scale_x_continuous(labels = labels, breaks = breaks) +
		labs(
			x = "Day of Year",
			y = "CRPS",
			color = "Data in model",
			title = if_else(d == "CARY", "Observed weather", "Forecasted weather")
		) +
		theme_bw() +
		annotate(
			geom = "rect",
			xmin = xmin,
			xmax = xmax,
			ymin = 0,
			ymax = Inf,
			alpha = 0.25
		)
	gg
}

g1 <- plot_phenology_scores("CARY")
g2 <- plot_phenology_scores("NMME")

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
	path = "Plots/hindcast/manuscript",
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
	df <- read_csv(file.path(dir_top, state_files[i])) %>% suppressMessages()
	start_date <- min(df$time)
	group <- str_extract(state_files[i], "(?<=\\/)[[:graph:]]*(?=/2)")
	model <- str_extract(group, "[[:alpha:]]*(?=\\/)")
	# if(is.na(model)) next
	experiment <- str_extract(group, "(?<=\\/)[[:graph:]]*")
	experiment <- gsub("nmme_nmme", "nmme", experiment)
	p_from <- str_extract(state_files[i], "(?<=paramsFrom_)[[:alpha:]]*(?=\\/)")
	t_from <- str_extract(state_files[i], "(?<=ticksFrom_)[[:alpha:]]*(?=\\_)")

	df_mutate <- df %>%
		filter(lifeStage != "dormant") %>%
		mutate(
			horizon = as.numeric(time - start_date),
			ticksFrom = t_from,
			paramsFrom = p_from,
			experiment = experiment,
			model = model,
			start_date = start_date
		) %>%
		mutate(
			lifeStage = if_else(lifeStage == "larvae", "Larvae", lifeStage),
			lifeStage = if_else(lifeStage == "nymphs", "Nymphs", lifeStage),
			lifeStage = if_else(lifeStage == "adults", "Adults", lifeStage)
		) %>%
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

# null_fx1 <- read_csv("Data/ticksFrom_Tea_paramsFrom_Tea/null/null/forecastQuants.csv")
# null_fx2 <- read_csv("Data/ticksFrom_Henry_paramsFrom_Henry/null/null/forecastQuants.csv")
null_fx3 <- read_csv(
	"Data/ticksFrom_Green_paramsFrom_Green/null/null/forecastQuants.csv"
)

# null_fx_bind <- bind_rows(null_fx1, null_fx2, null_fx3)

start_dates <- sort(unique(quants$start_date), decreasing = FALSE)
start_dates <- start_dates[grep("2018", start_dates)]
g <- list()
for (i in seq_along(start_dates)) {
	start_date <- start_dates[i]
	# "2018-09-05"

	dfq <- quants |>
		filter(start_date == start_date) |>
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
		filter(start_date == start_date, site == dfq$ticksFrom[1]) |>
		select(fx, ymin, ymax, time, count, lifeStage) |>
		mutate(remove = "Null", mice = "Null", driver = "Null")

	obs <- null |>
		select(time, count, lifeStage)

	g[[i]] <- bind_rows(null, dfq) |>
		left_join(obs) |>
		filter(
			lifeStage == "Nymphs",
			# remove == "Larvae" | remove == "Null",
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
	path = "Plots/hindcast/manuscript",
	width = 18,
	height = 24,
	units = "cm"
)

## transferability

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

transfer_plot <- function(ls) {
	transfer_cis |>
		filter(lifeStage == ls) |>
		ggplot() +
		aes(xmin = lwr.ci, x = med, xmax = upr.ci, y = instance, color = transfer) +
		geom_point(position = position_dodge(width = w)) +
		geom_linerange(position = position_dodge(width = w)) +
		facet_grid(frame ~ driver) +
		scale_color_manual(values = transfer_cols) +
		labs(x = "CRPS", y = "Data in model", color = "Parameter-site match") +
		theme_bw()
}


transfer_plot("Nymphs")
ggsave(
	"transferability_nymphs.jpeg",
	dpi = "retina",
	path = "Plots/hindcast/manuscript",
	width = 18,
	height = 10,
	units = "cm"
)

transfer_plot("Larvae")
ggsave(
	"transferability_larvae.jpeg",
	dpi = "retina",
	path = "Plots/hindcast/manuscript",
	width = 18,
	height = 10,
	units = "cm"
)

transfer_plot("Adults")
ggsave(
	"transferability_adults.jpeg",
	dpi = "retina",
	path = "Plots/hindcast/manuscript",
	width = 18,
	height = 10,
	units = "cm"
)


np <- crps_period |>
	filter(lifeStage == "Nymphs") |>
	mutate(
		frame = if_else(year(start_date) == year(time), "Subannual", "Interannual")
	) |>
	select(score, driver, remove, mice, period, horizon, frame) |>
	mutate(instance = paste(driver, remove, mice, sep = " "))

np_quants <- np |>
	filter(!grepl("NMME", instance)) |>

	group_by(period, instance, frame) |>
	summarise(
		`5%` = quantile(score, 0.05),
		`25%` = quantile(score, 0.25),
		`50%` = quantile(score, 0.5),
		`75%` = quantile(score, 0.75),
		`95%` = quantile(score, 0.95)
	) |>
	ungroup()

np_quants |>
	ggplot() +
	aes(y = instance, x = `50%`, xmin = `5%`, xmax = `95%`) +
	geom_linerange() +
	geom_linerange(aes(xmin = `25%`, xmax = `75%`), linewidth = 2) +
	geom_point(size = 3.75, color = "black") +
	geom_point(size = 3, fill = "white", color = "white") +
	facet_grid(frame ~ period) +
	labs(x = "CRPS", y = "Data streams in model") +
	theme_bw() +
	theme(
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12),
		strip.text = element_text(size = 14)
	)

ggsave("allCRPS.jpeg", dpi = "print", path = "Plots/hindcast/manuscript")


np |>
	filter(!grepl("NMME", instance)) |>
	mutate(
		instance = stringr::str_replace(instance, "CARY ", ""),
		instance = stringr::str_replace(instance, "Null Null ", "")
	) |>
	ggplot() +
	aes(x = score, y = instance) +
	geom_boxplot(outliers = FALSE) +
	facet_grid(frame ~ period) +
	theme_bw()


# data removal ----

remove_quants <- np |>
	filter(remove != "Null") |>
	group_by(remove, frame, period, driver) |>
	summarise(
		`5%` = quantile(score, 0.05),
		`25%` = quantile(score, 0.25),
		`50%` = quantile(score, 0.5),
		`75%` = quantile(score, 0.75),
		`95%` = quantile(score, 0.95)
	)

w <- 0.5
remove_quants |>
	ggplot() +
	aes(y = frame, x = `50%`, xmin = `5%`, xmax = `95%`, color = remove) +
	geom_linerange(position = position_dodge(width = w)) +
	geom_linerange(
		aes(xmin = `25%`, xmax = `75%`),
		linewidth = 2,
		position = position_dodge(width = w)
	) +
	geom_point(size = 3.75, position = position_dodge(width = w)) +
	geom_point(size = 3, fill = "white", position = position_dodge(width = w)) +
	facet_grid(driver ~ period) +
	labs(x = "CRPS", y = "Data streams in model") +
	theme_bw() +
	theme(
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12),
		strip.text = element_text(size = 14)
	)

mice_quants <- np |>
	filter(remove != "Null") |>
	group_by(mice, frame, period, driver) |>
	summarise(
		`5%` = quantile(score, 0.05),
		`25%` = quantile(score, 0.25),
		`50%` = quantile(score, 0.5),
		`75%` = quantile(score, 0.75),
		`95%` = quantile(score, 0.95)
	)

mice_quants |>
	ggplot() +
	aes(y = frame, x = `50%`, xmin = `5%`, xmax = `95%`, color = mice) +
	geom_linerange(position = position_dodge(width = w)) +
	geom_linerange(
		aes(xmin = `25%`, xmax = `75%`),
		linewidth = 2,
		position = position_dodge(width = w)
	) +
	geom_point(size = 3.75, position = position_dodge(width = w)) +
	geom_point(size = 3, fill = "white", position = position_dodge(width = w)) +
	facet_grid(driver ~ period) +
	labs(x = "CRPS", y = "Data streams in model") +
	theme_bw() +
	theme(
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12),
		strip.text = element_text(size = 14)
	)


save_gg <- function(dest, gg, path) {
	if (!dir.exists(path)) {
		dir.create(path, showWarnings = FALSE, recursive = TRUE)
	}
	ggsave(
		filename = dest,
		plot = gg,
		device = "jpeg",
		path = path,
		width = 7,
		height = 5,
		units = "in",
		dpi = "retina",
		bg = "white"
	)
	# dev.off()
}


fx_scores <- read_csv(file.path(dir_data, "allForecastScores.csv"))

crps_process <- fx_scores |>
	filter(metric == "crps", model != "Null") |>
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


site_x <- c("Green", "Henry", "Tea")
crps_null <- tibble()
for (i in 1:3) {
	site_params <- paste0("ticksFrom_", site_x[i], "_paramsFrom_", site_x[i])
	dir_read <- file.path(dir_data, site_params, "null", "null")
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
		horizon = as.numeric(time - start_date)
	) |>
	rename(ticksFrom = site)

crps <- bind_rows(crps_process, crps_null)

cary_met <- read_csv("Data/Cary_Met_Data_Daily.csv")
met_pre_2021 <- cary_met %>%
	slice(-c((nrow(cary_met) - 364):nrow(cary_met))) %>%
	mutate(DATE = mdy(DATE)) %>%
	filter(DATE >= "1995-01-01")
met_2021 <- cary_met %>%
	slice((nrow(cary_met) - 364):nrow(cary_met)) %>%
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

crps_period <- crps_cgdd |>
	filter(lifeStage == "Nymphs") |>
	mutate(period = if_else(cgdd >= 400 & cgdd <= 2500, "Questing", "Dormant"))


tick_grid |>
	left_join(met) |>
	filter(lifeStage == "Nymphs", cgdd >= 400 & cgdd <= 2500) |>
	group_by(ticksFrom, year) |>
	summarise(data = sum(data)) |>
	mutate(
		ym1 = c(0, data[-length(data)]),
		index = if_else(data > ym1, data / ym1, ym1 / data)
	) |>
	filter(!is.infinite(index)) |>
	group_by(ticksFrom) |>
	summarise(m = mean(index))

ggplot() +
	aes(x = year, y = index, color = ticksFrom) +
	geom_point() +
	geom_line()


w <- 3


iav <- iav |> rename(iav = cov)
crps_period <- crps_period |> mutate(year = year(time))

left_join(crps_period, iav) |>
	mutate(
		frame = if_else(year(start_date) == year(time), "Short", "Intermediate")
	) |>
	filter(period == "Questing") |>
	filter(driver != "Null") |>
	ggplot() +
	aes(x = iav, y = score, color = mice) +
	geom_smooth(method = "gam", se = FALSE) +
	facet_wrap(~frame) +
	theme_bw()


crps_id <- crps_period |>
	mutate(
		fx_id = paste(
			time,
			paramsFrom,
			ticksFrom,
			start_date,
			remove,
			mice,
			sep = "-"
		)
	)

nmme_ids <- crps_id |>
	filter(driver == "NMME") |>
	pull(fx_id) |>
	unique()

cary_ids <- crps_id |>
	filter(driver == "CARY", fx_id %in% nmme_ids) |>
	pull(fx_id) |>
	unique()

compare_ids <- nmme_ids[!nmme_ids %in% setdiff(nmme_ids, cary_ids)]

crps_compare <- crps_id |>
	filter(fx_id %in% compare_ids)

table(crps_compare$driver)


glimpse(crps_compare)

w <- 0.75
crps_compare |>
	filter(period == "Questing") |>
	mutate(
		frame = if_else(year(start_date) == year(time), "Short", "Intermediate")
	) |>
	group_by(driver, remove, mice, frame) |>
	summarise(
		med = median(score),
		low = quantile(score, 0.05),
		q1 = quantile(score, 0.25),
		q3 = quantile(score, 0.75),
		high = quantile(score, 0.95)
	) |>
	ggplot() +
	aes(x = driver, y = med, color = mice, linetype = remove) +
	geom_point(position = position_dodge(width = w), size = 4) +
	geom_linerange(
		aes(ymin = low, ymax = high),
		position = position_dodge(width = w),
		linewidth = 1.5
	) +
	facet_wrap(~frame) +
	theme_bw()

crps_compare |>
	filter(period == "Questing", remove == "Larvae") |>
	mutate(
		frame = if_else(year(start_date) == year(time), "Short", "Intermediate")
	) |>
	ggplot() +
	aes(x = frame, y = score, fill = driver, color = mice) +
	geom_boxplot(outliers = TRUE) +
	# facet_wrap(~ mice) +
	scale_fill_manual(values = c("white", "gray")) +
	theme_bw()


fx_samps <- read_csv("/Users/John.Foster/Downloads/allForecastSamples.csv") |>
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
	) |>
	mutate(
		fx_id = paste(
			time,
			paramsFrom,
			ticksFrom,
			start_date,
			remove,
			mice,
			sep = "-"
		)
	)


fx_cov <- fx_samps |>
	filter(lifeStage == "Nymphs") |>
	group_by(
		time,
		paramsFrom,
		horizon,
		ticksFrom,
		experiment,
		model,
		start_date
	) |>
	summarise(cov = sd(fx) / mean(fx), n = n())

cov_compare <- fx_cov |>
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
	) |>
	mutate(
		fx_id = paste(
			time,
			paramsFrom,
			ticksFrom,
			start_date,
			remove,
			mice,
			sep = "-"
		)
	) |>
	filter(fx_id %in% compare_ids) |>
	left_join(met) |>
	mutate(
		frame = if_else(year(start_date) == year(time), "Short", "Intermediate")
	) |>
	mutate(
		period = if_else(cgdd >= 400 & cgdd <= 2500, "Questing", "Dormant"),
		driver = if_else(model == "Null", "Null", driver),
		remove = if_else(model == "Null", "Null", remove),
		mice = if_else(model == "Null", "Null", mice)
	)

g1 <- cov_compare |>
	filter(model != "Null", period == "Questing") |>
	mutate(frame = factor(frame, levels = c("Short", "Intermediate"))) |>
	ggplot() +
	aes(x = driver, y = cov, fill = mice, color = remove) +
	geom_boxplot(outliers = FALSE, coef = 1.58) +
	facet_grid(~frame) +
	scale_fill_manual(values = c("white", "gray")) +
	scale_color_brewer(type = "qual", palette = 6) +
	labs(
		color = "Larvae data",
		y = "COV",
		x = "Weather data source",
		fill = "Mice data"
	) +
	theme_bw()

g2 <- fx_scores |>
	filter(lifeStage == "Nymphs", metric == "crps") |>
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
	) |>
	mutate(
		fx_id = paste(
			time,
			paramsFrom,
			ticksFrom,
			start_date,
			remove,
			mice,
			sep = "-"
		)
	) |>
	left_join(met) |>
	mutate(
		frame = if_else(year(start_date) == year(time), "Short", "Intermediate")
	) |>
	mutate(
		period = if_else(cgdd >= 400 & cgdd <= 2500, "Questing", "Dormant"),
		driver = if_else(model == "Null", "Null", driver),
		remove = if_else(model == "Null", "Null", remove),
		mice = if_else(model == "Null", "Null", mice)
	) |>
	filter(fx_id %in% compare_ids, period == "Questing", model != "Null") |>
	mutate(frame = factor(frame, levels = c("Short", "Intermediate"))) |>
	ggplot() +
	aes(x = driver, y = score, fill = mice, color = remove) +
	geom_violin() +
	geom_boxplot(
		outliers = FALSE,
		width = 0.12,
		position = position_dodge(width = 0.9)
	) +
	facet_grid(period ~ frame) +
	scale_fill_manual(values = c("white", "gray")) +
	scale_color_brewer(type = "qual", palette = 6) +
	labs(
		color = "Larvae data",
		y = "CRPS",
		x = "Weather data source",
		fill = "Mice data"
	) +
	theme_bw()

g2


ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "right")
