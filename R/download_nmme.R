# get NMME data from NOAA

library(ncdf4)
library(utils)
library(tidyverse)

# variables to grab
vars <- c(
	"tasmin", # daily min temp
	"tasmax", # daily max temp
	"hus", # specific humidity
	"pr" # precipitation rate
)


url <- "https://www.ncei.noaa.gov/data/north-american-multi-model-ensemble/access/ccsm4/"
nmme_base <- readLines(url)

potential_start <- potential_end <- potential_var <- potential_year <- vector()
for (i in seq_along(nmme_base)) {
	if (str_detect(nmme_base[i], "\\d{4}")) {
		year_grab <- str_extract(nmme_base[i], "\\d{4}")
		url_year <- paste0(url, year_grab, "/")
		Sys.sleep(0.5)
		nmme_page <- readLines(url_year)
		for (j in seq_along(nmme_page)) {
			if (any(str_detect(nmme_page[j], vars))) {
				var_match <- vars[str_detect(nmme_page[j], vars)]
				month_range <- str_extract(nmme_page[j], "\\d{8}-\\d{8}")
				start_month <- str_extract(month_range, "^\\d{8}")
				end_month <- str_extract(month_range, "\\d{8}$")
				potential_start <- c(potential_start, start_month)
				potential_end <- c(potential_end, end_month)
				potential_var <- c(potential_var, var_match)
				potential_year <- c(potential_year, year_grab)
			}
		}
	}
}

nmme_df <- tibble(
	potential_start = potential_start,
	potential_end = potential_end,
	potential_var = potential_var,
	potential_year = potential_year
) %>%
	distinct()

url_file_exist <- function(url) {
	HTTP_STATUS_OK <- 200 # nolint
	hd <- httr::HEAD(url)
	status <- hd$all_headers[[1]]$status
	exists <- status == HTTP_STATUS_OK
	list(exists = exists, status = status)
}

get_nmme <- function(
	nmme_var,
	nmme_year,
	nmme_start_month,
	nmme_end_month,
	dest_dir = "data/NMME"
) {
	for (i in 1:10) {
		# 10 files per var x month

		# construct file name
		file_name <- paste0(
			nmme_var,
			"_day_ccsm4_",
			nmme_start_month,
			"_r",
			i,
			"i1p1_",
			nmme_start_month,
			"-",
			nmme_end_month,
			".nc"
		)

		# the actual file from noaa
		url_2_get <- paste0(
			url,
			nmme_year,
			"/",
			file_name
		)

		file_avail <- url_file_exist(url_2_get)
		if (file_avail$exists) {
			save_dir <- file.path(dest_dir, "raw", nmme_start_month)
			if (!dir.exists(save_dir)) {
				dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
			}

			destfile <- file.path(save_dir, file_name)
			if (file.exists(destfile)) {
				message(file_name, " already exists")
			} else {
				message(
					"Downloading ",
					nmme_var,
					" for ",
					nmme_start_month,
					" ensemble ",
					i
				)

				download.file(
					url = url_2_get,
					destfile = destfile
				)
				Sys.sleep(2)
			}
		} else {
			message("File not available ", url_2_get)
			message("Status ", file_avail$status)
		}
	}
}

for (m in seq_along(potential_start)) {
	get_nmme(
		nmme_var = nmme_df$potential_var[m],
		nmme_year = nmme_df$potential_year[m],
		nmme_start_month = nmme_df$potential_start[m],
		nmme_end_month = nmme_df$potential_end[m]
	)
}
