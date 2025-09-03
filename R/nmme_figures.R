library(tidyverse)
library(lubridate)
library(MetBrewer)
library(ggpubr)

dir.top <- ""
dir.data <- "Data"
# dir.out <- "/projectnb/dietzelab/fosterj/FinalOut/Chapter2"
dir.analysis <- "Analysis"
dir.plot <- "plots"
dir.nmme <- "data/NMME"

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

# get the observed data from Cary and center to historical mean
cary.met <- read_csv(file.path(dir.data, "Cary_Met_Data_Daily.csv")) %>%
	suppressWarnings()
met.pre.2021 <- cary.met %>%
	slice(-c((nrow(cary.met) - 364):nrow(cary.met))) %>%
	mutate(Date = mdy(DATE)) %>%
	filter(Date >= "1995-01-01")
met.2021 <- cary.met %>%
	slice((nrow(cary.met) - 364):nrow(cary.met)) %>%
	mutate(Date = ymd(DATE))

source("R/function_scale_met_forecast.R")
hist.means <- scale_met_forecast()

hdf <- bind_rows(met.pre.2021, met.2021) %>%
	select(Date, MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) %>%
	arrange(Date) %>%
	select(Date, MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) %>%
	rename(
		`Max Temperature` = MAX_TEMP,
		`Max relative humidity` = MAX_RH,
		`Min relative humidity` = MIN_RH,
		`Daily precipitation` = TOT_PREC
	) %>%
	pivot_longer(
		cols = c(
			`Max Temperature`,
			`Max relative humidity`,
			`Min relative humidity`,
			`Daily precipitation`
		),
		names_to = "variable",
		values_to = "observed"
	) |>
	mutate(source = "Observed")

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


dir.out <- "C:/Users/John.Foster/Downloads/NMME"

files <- list.files(dir.out, recursive = TRUE)
# weather.files <- grep("CARY", files, value = TRUE)
# weather.files <- grep("pfSamples.csv", files, value = TRUE)

df.tick <- read_csv(file.path(dir.data, "Ticks2006_2021.csv"))
hindcast_seq <- df.tick |>
	filter(
		Grid %in% c("Green Control", "Henry Control", "Tea Control"),
		Date >= ymd("2018-01-01")
	) |>
	pull(Date) |>
	sort() |>
	unique()

pb = txtProgressBar(min = 1, max = length(hindcast_seq), style = 1)
all_df <- tibble()
for (i in seq_along(hindcast_seq)) {
	start.date <- hindcast_seq[i]

	nmme2get <- floor_date(ymd(start.date), unit = "month")
	nmme.name <- gsub("-", "", nmme2get)

	fff <- file.path(
		dir.nmme,
		nmme.name,
		"CARY",
		paste0("NOAANMME_", nmme.name, ".csv")
	)

	# check if nmme file exists (some file don't)
	if (!file.exists(fff)) {
		nmme2get <- floor_date(ymd(start.date) - 30, unit = "month")
		nmme.name <- gsub("-", "", nmme2get)
		fff <- file.path(
			dir.nmme,
			nmme.name,
			"CARY",
			paste0("NOAANMME_", nmme.name, ".csv")
		)
	}

	df.nmme <-
		read_csv(fff) %>%
		suppressMessages()

	target.nmme <- df.nmme %>%
		filter(time >= start.date) |>
		select(-tmin) %>%
		mutate(rhmin = rhmin * 100, rhmax = rhmax * 100)

	nmme.means <- target.nmme %>%
		pivot_longer(
			cols = c(tmax, precipitation, rhmin, rhmax),
			names_to = "variable"
		) %>%
		group_by(variable, ensemble) %>%
		summarise(mu = mean(value)) %>%
		pivot_wider(names_from = variable, values_from = mu) %>%
		mutate(
			bias.tmax = hist.means$means['MAX_TEMP'] - tmax,
			bias.rhmax = hist.means$means['MAX_RH'] - rhmax,
			bias.rhmin = hist.means$means['MIN_RH'] - rhmin
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

	# scale to historical period
	nmme.scaled <- nmme.correct %>%
		select(time, ensemble, tmax, rhmax, rhmin, precipitation) %>%
		mutate(
			# tmax = (tmax - hist.means$means['MAX_TEMP']) / hist.means$sds['MAX_TEMP'],
			# rhmax = (rhmax - hist.means$means['MAX_RH']) / hist.means$sds['MAX_RH'],
			# rhmin = (rhmin - hist.means$means['MIN_RH']) / hist.means$sds['MIN_RH'],
			# precipitation = (precipitation - hist.means$means['TOT_PREC']) / hist.means$sds['TOT_PREC']
		)

	df.nmme2 <- nmme.scaled |>
		rename(
			Date = time,
			`Max Temperature` = tmax,
			`Daily precipitation` = precipitation,
			`Max relative humidity` = rhmax,
			`Min relative humidity` = rhmin
		) |>
		arrange(ensemble, Date) |>
		group_by(ensemble) |>
		mutate("Cumulative precipitation" = cumsum(`Daily precipitation`)) |>
		ungroup() |>
		mutate(source = "Forecasted", ensemble = paste0("nmme_", ensemble)) |>
		pivot_longer(
			cols = -c(source, Date, ensemble),
			values_to = "observed",
			names_to = "variable"
		)

	nmme_dates <- df.nmme2 |>
		select(Date, source) |>
		distinct()

	good_dates <- hdf |>
		select(Date) |>
		distinct() |>
		left_join(nmme_dates) |>
		filter(!is.na(source)) |>
		pull(Date)

	if (length(good_dates) == 0) {
		next
	}

	total_precip <- hdf |>
		filter(
			Date >= start.date,
			Date %in% good_dates
		) |>
		filter(variable == "Daily precipitation") |>
		mutate(observed = cumsum(observed), variable = "Cumulative precipitation")

	plot_df <- hdf |>
		filter(Date %in% good_dates) |>
		mutate(ensemble = "cary_1") |>
		bind_rows(total_precip) |>
		bind_rows(df.nmme2) |>
		mutate(
			source = factor(source, levels = c("Observed", "Forecasted")),
			start.date = start.date,
			start.date = start.date
		)

	ensemble_fac <- plot_df |>
		pull(ensemble) |>
		unique() |>
		sort(decreasing = TRUE)

	plot_var <- function(v) {
		plot_df |>
			mutate(ensemble = factor(ensemble, levels = ensemble_fac)) |>
			filter(variable == v) |>
			ggplot() +
			aes(
				x = Date,
				y = observed,
				group = ensemble,
				color = source,
				linewidth = source,
				linetype = source
			) +
			geom_line() +
			facet_wrap(~variable, scales = "free_y") +
			scale_color_manual(
				values = c("Observed" = "red", "Forecasted" = "black")
			) +
			scale_linetype_manual(
				values = c("Observed" = "solid", "Forecasted" = "dotted")
			) +
			scale_linewidth_manual(values = c(0.2, 0.2)) +
			labs(
				color = "Weather data source",
				linewidth = "Weather data source",
				linetype = "Weather data source",
				x = element_blank()
			) +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
	}

	# g <- list()
	# g[[1]] <- plot_var("Max Temperature") + labs(y = "Degrees C.")
	# g[[2]] <- plot_var("Max relative humidity") + labs(y = "Percent")
	# g[[3]] <- plot_var("Min relative humidity") + labs(y = "Percent")
	# g[[4]] <- plot_var("Daily precipitation") + labs(y = "Millimeters")
	# g[[5]] <- plot_var("Cumulative precipitation") + labs(y = "Millimeters")
	#
	# ggarrange(plotlist = g, nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom", labels = "AUTO")
	#
	# ggsave(paste0(start.date, "nmme_vs_observed.jpeg"), dpi = "retina",
	#        path = "Plots/hindcast/manuscript/nmmeVsObs", width = 18, height = 12, units = "cm")

	setTxtProgressBar(pb, i)

	all_df <- bind_rows(all_df, plot_df)
}
close(pb)


cary <- all_df |>
	filter(source == "Observed") |>
	select(-ensemble, -source) |>
	distinct() |>
	mutate(start.date = ymd(start.date), horizon = as.numeric(Date - start.date))

nmme <- all_df |>
	filter(source == "Forecasted") |>
	select(-source) |>
	mutate(
		start.date = ymd(start.date),
		horizon = as.numeric(Date - start.date)
	) |>
	rename(forecast = observed)

met_join <- left_join(nmme, cary) |>
	filter(!is.na(observed), !is.na(forecast))

met_fx_horizon <- met_join |>
	group_by(variable, horizon) |>
	summarise(Bias = mean(forecast - observed), MAE = abs(Bias))

my_plot <- function(y, v, ylab) {
	g <- met_fx_horizon |>
		filter(variable == v) |>
		ggplot() +
		aes(x = horizon, y = .data[[y]]) +
		geom_point(size = 0.1) +
		geom_smooth(linewidth = 0.5) +
		facet_wrap(~variable) +
		labs(x = "Lead time (days)", y = ylab) +
		theme_bw()

	if (y == "Bias") {
		g <- g + geom_hline(yintercept = 0, linetype = "dashed")
	}

	g
}

g <- list()
g[[1]] <- my_plot("Bias", "Daily precipitation", "Millimeters")
g[[2]] <- my_plot("Bias", "Max relative humidity", "Percent RH")
g[[3]] <- my_plot("Bias", "Min relative humidity", "Percent RH")
g[[4]] <- my_plot("Bias", "Max Temperature", "Celsius")
g[[5]] <- my_plot("MAE", "Daily precipitation", "Millimeters")
g[[6]] <- my_plot("MAE", "Max relative humidity", "Percent RH")
g[[7]] <- my_plot("MAE", "Min relative humidity", "Percent RH")
g[[8]] <- my_plot("MAE", "Max Temperature", "Celsius")

ggarrange(plotlist = g, nrow = 2, ncol = 4, align = "hv", labels = "AUTO")

ggsave(
	paste0("nmmeAllMetrics.jpeg"),
	dpi = "retina",
	path = "plots",
	width = 18,
	height = 12,
	units = "cm"
)

met_fx_horizon |>
	pivot_longer(cols = c(Bias, MAE), names_to = "metric") |>
	filter(variable != "Cumulative precipitation") |>
	ggplot() +
	aes(x = horizon, y = value) +
	geom_point(size = 0.1) +
	geom_smooth(linewidth = 0.5) +
	facet_grid(metric ~ variable, scales = "free_y") +
	labs(x = "Forecast lead time (days)", y = "Value") +
	theme_bw()

df <- muf %>%
	pivot_longer(cols = -c(site, paramsFrom), names_to = "parameter") %>%
	mutate(
		variable = if_else(grepl("1]", parameter), "Max Temperature", "x"),
		variable = if_else(
			grepl("2]", parameter),
			"Max relative humidity",
			variable
		),
		variable = if_else(
			grepl("3]", parameter),
			"Min relative humidity",
			variable
		),
		variable = if_else(grepl("4]", parameter), "Precipitation", variable),
		horizon = as.numeric(str_extract(parameter, "(?<=muf\\[)\\d*")) - 1,
		value = if_else(
			variable == "Max Temperature",
			value * hist.means$sds['MAX_TEMP'] + hist.means$means['MAX_TEMP'],
			value
		),
		value = if_else(
			variable == "Max relative humidity",
			value * hist.means$sds['MAX_RH'] + hist.means$means['MAX_RH'],
			value
		),
		value = if_else(
			variable == "Min relative humidity",
			value * hist.means$sds['MIN_RH'] + hist.means$means['MIN_RH'],
			value
		),
		value = if_else(
			variable == "Precipitation",
			value * hist.means$sds['TOT_PREC'] + hist.means$means['TOT_PREC'],
			value
		)
	) %>%
	group_by(variable, horizon) %>%
	summarise(
		low = quantile(value, 0.05),
		med = quantile(value, 0.5),
		high = quantile(value, 0.95)
	) %>%
	mutate(
		source = "Estimated",
		start.date = ymd(start.date),
		Date = start.date + horizon
	)


df.p <- left_join(df, hdf, by = c("Date", "variable"))

gg <- df.p %>%
	# filter(variable == "rhmax") %>%
	ggplot() +
	aes(x = Date) +
	geom_line(aes(y = observed, linetype = "Observed"), size = 0.25) +
	geom_ribbon(aes(ymin = low, ymax = high, fill = "Estimated"), alpha = 0.8) +
	scale_fill_manual(values = "#fc8d59") +
	labs(y = "Value", linetype = "", fill = "") +
	facet_wrap(~variable) +
	theme_pubr() +
	theme(
		legend.position = "bottom",
		axis.text.x = element_text(angle = 60, vjust = 0.5)
	)
