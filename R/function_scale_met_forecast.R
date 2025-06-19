#' Scale meteorological data for forecasting
#' the scaling is based on the calibration period (1995-2005)

library(plantecophys)
library(dplyr)

scale_met_forecast <- function() {
	met <- read.csv("Data/Cary_Met_Data_Daily.csv")
	met$DATE <- lubridate::mdy(met$DATE)
	met <- met %>%
		filter(DATE >= "1995-05-02") %>%
		filter(DATE <= "2005-08-17") %>%
		select(c("MAX_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC"))

	met_means <- apply(met, 2, mean, na.rm = TRUE)
	met_sd <- apply(met, 2, sd, na.rm = TRUE)

	met_scale <- apply(met, 2, scale)
	scale.max <- c(
		10000,
		(100 - met_means["MAX_RH"]) / met_sd["MAX_RH"],
		(100 - met_means["MIN_RH"]) / met_sd["MIN_RH"],
		10000
	)

	scale.min <- c(
		-10000,
		(0 - met_means["MAX_RH"]) / met_sd["MAX_RH"],
		(0 - met_means["MIN_RH"]) / met_sd["MIN_RH"],
		(0 - met_means["TOT_PREC"]) / met_sd["TOT_PREC"]
	)

	names(scale.max) <- c("MAX_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")
	names(scale.min) <- c("MAX_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")

	list(
		means = met_means,
		sds = met_sd,
		scale.max = scale.max,
		scale.min = scale.min
	)
}
