# miscellaneous functions used throughout the workkflow

create_ic <- function(df_mu, df_tau, start_date) {
	IC <- matrix(NA, 4, 2)
	IC[1, 1] <- df_mu |>
		filter(lifeStage == "larvae") |>
		pull(as.character(fx_start_date))
	IC[1, 2] <- df_tau |>
		filter(lifeStage == "larvae") |>
		pull(as.character(fx_start_date))
	IC[2, 1] <- df_mu |>
		filter(lifeStage == "dormant") |>
		pull(as.character(fx_start_date))
	IC[2, 2] <- df_tau |>
		filter(lifeStage == "dormant") |>
		pull(as.character(fx_start_date))
	IC[3, 1] <- df_mu |>
		filter(lifeStage == "nymphs") |>
		pull(as.character(fx_start_date))
	IC[3, 2] <- df_tau |>
		filter(lifeStage == "nymphs") |>
		pull(as.character(fx_start_date))
	IC[4, 1] <- df_mu |>
		filter(lifeStage == "adults") |>
		pull(as.character(fx_start_date))
	IC[4, 2] <- df_tau |>
		filter(lifeStage == "adults") |>
		pull(as.character(fx_start_date))

	IC
}

# function to approximate moment the inverse gamma
inv_gamma_mm <- function(x) {
	mu <- mean(x)
	v <- var(x)
	alpha <- (mu^2 / v) + 2
	beta <- mu * ((mu^2 / v) + 1)
	c("alpha" = alpha, "beta" = beta)
}

get_prior <- function(name) {
	pr <- numeric(2)
	xx <- params_stats |>
		filter(parameter == name)
	pr[1] <- xx |> pull(mu)
	pr[2] <- xx |> pull(tau)
	pr
}
