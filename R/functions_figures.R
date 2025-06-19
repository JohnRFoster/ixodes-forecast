### functions for plots

make_plot_df <- function(pred_data) {
	doy <- pred_data |>
		filter(y == "score", instance == "Day-of-year") |>
		rename(value_doy = value) |>
		select(value_doy, x, period, driver)

	proc <- pred_data |>
		filter(y == "score", instance != "Day-of-year")

	left_join(proc, doy) |>
		mutate(delta = value - value_doy)
}

# make gams for questing and dormant phases by lifstage
gam_predict <- function(df, ls, t, w, r, m) {
	df_gam1 <- df |>
		filter(
			lifeStage == ls,
			driver == w,
			remove == r,
			mice == m,
			period == "Questing"
		)

	df_gam2 <- df |>
		filter(
			lifeStage == ls,
			driver == w,
			remove == r,
			mice == m,
			period == "Dormant"
		)

	if (t == "score") {
		fit1 <- gam(
			target ~ s(horizon, k = 5),
			data = df_gam1,
			method = "REML",
			family = gaussian(link = "log")
		)
		fit2 <- gam(
			target ~ s(horizon, k = 5),
			data = df_gam2,
			method = "REML",
			family = gaussian(link = "log")
		)
	} else {
		fit1 <- gam(target ~ s(horizon, k = 5), data = df_gam1, method = "REML")
		fit2 <- gam(target ~ s(horizon, k = 5), data = df_gam2, method = "REML")
	}

	h <- 0:365
	mq <- predict.gam(fit1, tibble(horizon = h), type = "response")
	md <- predict.gam(fit2, tibble(horizon = h), type = "response")

	pq <- tibble(
		value = mq,
		lifeStage = ls,
		x = h,
		driver = w,
		remove = r,
		mice = m,
		period = "Questing"
	)

	pd <- tibble(
		value = md,
		lifeStage = ls,
		x = h,
		driver = w,
		remove = r,
		mice = m,
		period = "Dormant"
	)

	bind_rows(pq, pd) |>
		mutate(y = t)
}

make_prediction_data <- function(LS, df_skill, df_delta) {
	tmp <- bind_rows(
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "CARY",
			r = "Larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "CARY",
			r = "Larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "CARY",
			r = "No larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "CARY",
			r = "No larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "NMME",
			r = "Larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "NMME",
			r = "Larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "NMME",
			r = "No larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "NMME",
			r = "No larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_skill,
			ls = LS,
			t = "score",
			w = "Null",
			r = "Null",
			m = "Null"
		),

		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "CARY",
			r = "Larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "CARY",
			r = "Larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "CARY",
			r = "No larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "CARY",
			r = "No larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "NMME",
			r = "Larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "NMME",
			r = "Larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "NMME",
			r = "No larvae",
			m = "Mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "NMME",
			r = "No larvae",
			m = "No mice"
		),
		gam_predict(
			df = df_delta,
			ls = LS,
			t = "delta",
			w = "Null",
			r = "Null",
			m = "Null"
		)
	) |>
		mutate(
			driver = if_else(driver == "CARY", "Obs. weather", driver),
			driver = if_else(driver == "NMME", "FX weather", driver)
		)

	n1 <- tmp |>
		filter(driver == "Null") |>
		mutate(instance = "Day-of-year", driver = "Obs. weather")
	n2 <- tmp |>
		filter(driver == "Null") |>
		mutate(instance = "Day-of-year", driver = "FX weather")

	tmp <- tmp |>
		filter(driver != "Null") |>
		mutate(instance = paste(remove, mice, sep = " ")) |>
		bind_rows(n1) |>
		bind_rows(n2)
}

instance_cols <- c(
	"Day-of-year" = "black",
	"Larvae Mice" = "#ca0020",
	"Larvae No mice" = "#f4a582",
	"No larvae Mice" = "#92c5de",
	"No larvae No mice" = "#0571b0"
)

my_theme <- function() {
	theme(
		axis.title = element_text(size = 8),
		axis.text = element_text(size = 8),
		legend.text = element_text(size = 10),
		legend.title = element_text(size = 10),
		strip.text = element_text(size = 8)
	)
}

plot_gam <- function(df, target, d, vline) {
	g <- df |>
		filter(y == target, driver == d) |>
		ggplot() +
		aes(x = x, y = delta, color = instance) +
		geom_line(linewidth = 1) +
		geom_vline(xintercept = vline) +
		geom_hline(yintercept = 0, linetype = "dashed") +
		facet_grid(driver ~ period) +
		scale_color_manual(values = instance_cols) +
		labs(
			color = "Data in model",
			x = "Lead time (days)",
			y = if_else(target == "score", "Relative CRPS", "CRPS - CRPS0")
		) +
		theme_bw() +
		my_theme()
}

plot_ticks_on_drag_cloths <- function(df_tick, df_skill, ls) {
	tmp <- left_join(df_skill, df_tick) |>
		filter(driver == "CARY", lifeStage == ls)

	if (ls == "Nymphs") {
		tmp2 <- tmp |>
			mutate(
				period = if_else(cgdd >= 400 & cgdd <= 2500, "Questing", "Dormant")
			)
	} else if (ls == "Adults") {
		tmp2 <- tmp |>
			mutate(
				period = if_else(cgdd <= 1000 | cgdd >= 2500, "Questing", "Dormant")
			)
	} else if (ls == "Larvae") {
		tmp2 <- tmp |>
			mutate(
				period = if_else(cgdd >= 1400 & cgdd <= 2500, "Questing", "Dormant")
			)
	}

	tmp2 |>
		ggplot() +
		aes(x = horizon, y = data, color = period) +
		geom_point(size = 0.2) +
		geom_vline(xintercept = 175) +
		labs(
			x = "Lead time (days)",
			y = "Ticks / 450 sq. m",
			color = paste0(ls, " on\ndrag cloths")
		) +
		scale_color_brewer(type = "qual", palette = 2) +
		theme_bw() +
		my_theme()
}

get_limit <- function(df) {
	n3 <- df |>
		filter(instance == "Day-of-year", y == "score") |>
		rename(null_value = value) |>
		select(null_value, x, driver, period) |>
		distinct()

	limit <- df |>
		filter(instance != "Day-of-year", y == "score") |>
		select(-y) |>
		left_join(n3, by = join_by(x, driver, period)) |>
		filter(!is.na(null_value)) |>
		mutate(
			limit = null_value - value,
			frame = if_else(x < 175, "Subannual", "Interannual")
		) |>
		group_by(x, driver, frame, mice, remove, period) |>
		summarise(mu = mean(limit)) |>
		ungroup() |>
		suppressMessages()

	tmp <- limit |>
		filter(driver == d, frame == f, mice == m, remove == r, period == p) |>
		arrange(x)

	idx <- max(which(tmp$mu >= 0))
	tmp$x[idx]
}


parameter_names <- function(df) {
	df |>
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
			driver = if_else(
				parameter %in% c("beta[13]", "beta[14]"),
				"Mice",
				driver
			),
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
}

plot_phenology_scores <- function(d, ls, df_process, df_null) {
	rect <- tibble(
		lifeStage = c("Larvae", "Nymphs", "Adults"),
		xmin = c(yday("2022-07-01"), yday("2022-05-01"), yday("2022-08-01")),
		xmax = c(yday("2022-10-01"), yday("2022-08-01"), yday("2022-05-01"))
	)

	df_null <- df_null |>
		filter(lifeStage == ls)

	gg <- df_process |>
		filter(driver == d, lifeStage == ls) |>
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
		bind_rows(df_null) |>
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

transfer_plot <- function(df, ls) {
	df |>
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
