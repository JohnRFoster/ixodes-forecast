library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lubridate)

my_theme <- function() {
	theme(
		axis.title = element_text(size = 8),
		axis.text = element_text(size = 8),
		legend.text = element_text(size = 10),
		legend.title = element_text(size = 10),
		strip.text = element_text(size = 6)
	)
}

dir_data <- "data"
dir_out <- "out"
dir_plot <- "plots"

site_vec <- c("Green", "Henry", "Tea")

fx_quants <- read_csv("data/allForecastQuantiles.csv") |>
	filter(model != "null")
glimpse(fx_quants)

resid_cols <- c(
	"start.date",
	"time",
	"lifeStage",
	"count",
	"fx",
	"ymin",
	"ymax",
	"var",
	"in_95",
	"residual",
	"driver",
	"remove",
	"mice",
	"frame",
	"ticksFrom",
	"paramsFrom"
)

fx_resids <- fx_quants |>
	filter(grp == "95%") |>
	distinct() |>
	mutate(
		residual = data - q0.5,
		in_95 = if_else(data >= ymin & data <= ymax, 1, 0),
		driver = if_else(grepl("nmme", experiment), "NMME", "CARY"),
		remove = if_else(grepl("remove", experiment), "No larvae", "Larvae"),
		mice = if_else(grepl("mna", experiment), "Mice", "No mice"),
		model = if_else(
			model == "WithWeatherAndMiceGlobal",
			"Weather and Mice",
			model
		),
		frame = if_else(horizon < 175, "Subannual", "Interannual")
	) |>
	rename(fx = q0.5, count = data) |>
	select(all_of(resid_cols))

dir_null <- "Null"
site_x <- c("Green", "Henry", "Tea")
crps_null <- tibble()
for (i in 1:3) {
	site_params <- paste0("ticksFrom_", site_x[i], "_paramsFrom_", site_x[i])
	dir_read <- file.path("out", dir_null, site_params, "null", "null")
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
		residual = count - fx,
		in_95 = if_else(count >= ymin & count <= ymax, 1, 0),
		horizon = as.numeric(time - start.date),
		frame = if_else(horizon < 175, "Subannual", "Interannual")
	) |>
	rename(ticksFrom = site) |>
	select(all_of(resid_cols))


cary_met <- read_csv("Data/Cary_Met_Data_Daily.csv")
met_pre_2021 <- cary_met |>
	slice(-c((nrow(cary_met) - 364):nrow(cary_met))) |>
	mutate(DATE = mdy(DATE)) |>
	filter(DATE >= "1995-01-01")
met_2021 <- cary_met |>
	slice((nrow(cary_met) - 364):nrow(cary_met)) |>
	mutate(DATE = ymd(DATE))

met <- bind_rows(met_pre_2021, met_2021) |>
	arrange(DATE) |>
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

ticks_clean <- read_csv("data/tick_cleaned")

df_tick <- read_csv(file.path(dir_data, "Ticks2006_2021.csv"))
tick_grid <- df_tick |>
	rename(time = Date) |>
	mutate(site = gsub(" Control", "", Grid), year = year(time)) |>
	filter(site %in% c("Green", "Henry", "Tea")) |>
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

ticks_phase <- ticks_clean |>
	filter(grepl("Control", Grid)) |>
	select(Grid, DATE, starts_with("n_")) |>
	pivot_longer(
		cols = starts_with("n_"),
		names_to = "lifeStage",
		values_to = "count"
	) |>
	rename(time = DATE, ticksFrom = Grid, data = count) |>
	mutate(year = year(time)) |>
	mutate(
		lifeStage = case_when(
			lifeStage == "n_larvae" ~ "Larvae",
			lifeStage == "n_nymphs" ~ "Nymphs",
			lifeStage == "n_adults" ~ "Adults"
		)
	) |>
	bind_rows(tick_grid) |>
	left_join(met) |>
	filter(!is.na(cgdd)) |>
	mutate(
		phase = case_when(
			lifeStage == "Nymphs" ~ if_else(
				cgdd >= 400 & cgdd <= 2500,
				"Questing",
				"Dormant"
			),
			lifeStage == "Adults" ~ if_else(
				cgdd <= 1000 | cgdd >= 2500,
				"Questing",
				"Dormant"
			),
			lifeStage == "Larvae" ~ if_else(
				cgdd >= 1400 & cgdd <= 2500,
				"Questing",
				"Dormant"
			)
		)
	)

nymph_q_dates <- ticks_phase |>
	filter(phase == "Questing", lifeStage == "Nymphs") |>
	pull(time) |>
	unique()

all_resids <- bind_rows(fx_resids, crps_null) |>
	mutate(
		mu = fx,
		y_pred = NA,
		ymin_q = NA,
		ymax_q = NA,
		bayes_p = NA
	)

all_resids <- left_join(all_resids, ticks_phase) |>
	select(-count) |>
	distinct()

sim_y <- function(mu, v, y_obs, n = 1000) {
	set.seed(1205)

	lambda <- pmax(0, rnorm(n, mu, sqrt(v)))
	y_pred <- rpois(n, lambda)

	bayes_p <- if_else(
		y_obs == 0,
		mean(y_pred == 0),
		mean(y_pred >= y_obs)
	)

	list(
		bayes_p = bayes_p,
		y_median = median(y_pred),
		y_low = quantile(y_pred, 0.025),
		y_high = quantile(y_pred, 0.975)
	)
}

pb <- txtProgressBar(min = 1, max = nrow(all_resids), style = 3)
for (i in seq_len(nrow(all_resids))) {
	m_i <- all_resids$mu[i]
	v_i <- all_resids$var[i]
	y_obs_i <- all_resids$data[i]
	y_dist_i <- sim_y(m_i, v_i, y_obs_i)
	all_resids$y_pred[i] <- y_dist_i$y_median
	all_resids$ymin_q[i] <- y_dist_i$y_low
	all_resids$ymax_q[i] <- y_dist_i$y_high
	all_resids$bayes_p[i] <- y_dist_i$bayes_p
	setTxtProgressBar(pb, i)
}
close(pb)

all_resids <- all_resids |>
	mutate(
		b = y_pred - data,
		in_95_q = if_else(data >= ymin_q & data <= ymax_q, 1, 0),
		low_high = if_else(data > ymax_q, "fx too low", "fx too high"),
		low_high = if_else(in_95_q == 1, "in", low_high)
	)

table(all_resids$low_high)
# mu = fx                high: 3622 in: 38134 low: 3742

nmme_range <- all_resids |>
	filter(driver == "NMME") |>
	pull(start.date) |>
	range()

pm <- function(df, l, f, r, m, d) {
	df |>
		filter(
			lifeStage == l,
			frame == f,
			remove == r,
			mice == m,
			driver == d
		)
}

all_resids |>
	group_by(remove, mice, driver, lifeStage) |>
	reframe(
		n = n(),
		RMSE = round(sqrt(sum(b^2) / n), 1),
		bias = round(mean(b), 1),
		coverage = round(sum(in_95_q) / n * 100, 1)
	) |>
	arrange(lifeStage, remove, mice, driver)


instance <- all_resids |>
	mutate(
		remove = if_else(remove == "Larvae", "L", "NL"),
		mice = if_else(mice == "Mice", "M", "NM"),
		data_instance = paste(driver, remove, mice, sep = "-"),
		data_instance = if_else(grepl("Null", data_instance), "Null", data_instance)
	)

# need to put the above frequencies calculations into this function
hex_horizon <- function(df, bs) {
	tmp <- df |>
		mutate(
			bayes_p_bin = cut(
				bayes_p,
				breaks = seq(0, 1, by = 0.2),
				include.lowest = TRUE
			),
			horizon = as.numeric(time - start.date),
			horizon_bin = if_else(horizon < 175, "SA", "IA")
		)

	plyr::count(tmp$bayes_p_bin)
	plyr::count(tmp$horizon_bin)

	horizon_sample_size <- tmp |>
		group_by(data_instance, horizon_bin) |>
		reframe(n = n())

	horizon_freq <- tmp |>
		group_by(data_instance, horizon_bin, bayes_p_bin) |>
		reframe(nh = n()) |>
		left_join(horizon_sample_size) |>
		mutate(p = nh / n)

	horizon_freq |>
		ggplot() +
		aes(x = bayes_p_bin, y = horizon_bin, fill = p) +
		facet_wrap(~data_instance, scales = "free_y") +
		geom_tile() +
		labs(x = "Bayesian p-value", y = "Horizon", fill = "Proportion") +
		scale_fill_gradientn(
			colors = c("#0571b0", "white", "#fddbc7", "#ef8a62", "#b2182b"),
			limits = c(0, 0.8),
			breaks = c(0, 0.2, 0.4, 0.6, 0.8)
		) +
		theme_bw() +
		labs_pubr(base_size = bs) +
		theme(
			strip.text = element_text(size = 8),
			axis.text.x = element_text(angle = 90, hjust = 1)
		)
}

g <- list()
bs <- 8

g[[1]] <- instance |>
	hex_horizon(bs) +
	labs(title = "All forecasts") +
	my_theme()

g[[2]] <- instance |>
	filter(phase == "Questing") |>
	filter(lifeStage == "Nymphs") |>
	hex_horizon(bs) +
	labs(title = "Questing nymph forecasts") +
	my_theme()

hex_doy <- function(df, bs) {
	doy_seq <- c(1, 90, 181, 273, 365)
	doy_labels <- c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec")

	tmp <- df |>
		mutate(
			bayes_p_bin = cut(
				bayes_p,
				breaks = seq(0, 1, by = 0.2),
				include.lowest = TRUE
			),
			doy = yday(time),
			doy_bin = cut(doy, breaks = doy_seq, labels = doy_labels)
		)

	plyr::count(tmp$bayes_p_bin)
	plyr::count(tmp$doy_bin)

	doy_sample_size <- tmp |>
		group_by(data_instance, doy_bin) |>
		reframe(n = n())

	doy_freq <- tmp |>
		group_by(data_instance, doy_bin, bayes_p_bin) |>
		reframe(nh = n()) |>
		left_join(doy_sample_size) |>
		mutate(p = nh / n)

	doy_freq |>
		ggplot() +
		aes(x = bayes_p_bin, y = doy_bin, fill = p) +
		facet_wrap(~data_instance, scales = "free_y") +
		geom_tile() +
		labs(x = "Bayesian p-value", y = "Month", fill = "Proportion") +
		scale_fill_gradientn(
			colors = c("#0571b0", "white", "#fddbc7", "#ef8a62", "#b2182b"),
			limits = c(0, 0.8),
			breaks = c(0, 0.2, 0.4, 0.6, 0.8)
		) +
		theme_bw() +
		labs_pubr(base_size = bs) +
		theme(
			strip.text = element_text(size = 8),
			axis.text.x = element_text(angle = 90, hjust = 1)
		)
}

g[[3]] <- instance |>
	hex_doy(bs) +
	labs(title = "All forecasts") +
	my_theme()

g[[4]] <- instance |>
	filter(phase == "Questing") |>
	filter(lifeStage == "Nymphs") |>
	hex_doy(bs) +
	labs(title = "Questing nymph forecasts") +
	my_theme()

ggarrange(
	plotlist = g,
	nrow = 2,
	ncol = 2,
	labels = "AUTO",
	common.legend = TRUE,
	legend = "bottom",
	align = "hv"
)

ggsave(
	"figure_S1.tiff",
	dpi = 600,
	path = dir_plot,
	width = 18,
	height = 18,
	units = "cm"
)
