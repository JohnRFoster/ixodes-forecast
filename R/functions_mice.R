# functions that deal with the mouse mark-racpture data

library(tidyverse)
library(lubridate)

known_states <- function(ch) {
	state <- ch
	for (i in seq_len(dim(ch)[1])) {
		n1 <- min(which(ch[i, ] != 0))
		n2 <- max(which(ch[i, ] != 0))
		if (n2 > n1) {
			state[i, n1:n2] <- 1
		}
	}
	state
}

ch_cary <- function(grid, path = "") {
	file <- "data/Cary_mouse.csv"

	dat <- read.csv(
		paste(file, sep = ""),
		na.strings = c(c("", " ", "      "), "NA")
	) # read in data

	smam <- dat[, c(
		"Grid",
		"Full.Date.1",
		"Full.Date.2",
		"Day.1",
		"Day.2",
		"Tag..",
		"Fate"
	)]

	alive <- c(1, 2) # codes for alive individuals
	smam <- subset(smam, Fate == 1 | Fate == 2)

	smam[4:5] <- as.integer(!is.na(smam[4:5])) # converts trap histories to 1 or 0

	smam$Full.Date.1 <- as.Date(smam$Full.Date.1) # convert factors to dates
	smam$Full.Date.2 <- as.Date(smam$Full.Date.2) # convert factors to dates

	m <- subset(smam, Grid == grid)
	m <- m[, c("Tag..", "Day.1", "Day.2", "Full.Date.1", "Full.Date.2")]
	m <- subset(m, !is.na(Tag..))

	m$Tag.. <- as.character(m$Tag..)
	m_ls <- split(m, m$Full.Date.1) # split on first capture date for sampling occasion
	for (i in seq_along(m_ls)) {
		m_ls[[i]] <- m_ls[[i]][c(-4, -5)]
	}

	day_1 <- as.character(unique(m$Full.Date.1)) # 1st capture date of sampling occasion
	day_2 <- as.character(unique(m$Full.Date.2)) # 2nd capture date of sampling occasion
	days <- c(rbind(day_1, day_2)) # vector of unique trapping days (for colnames)

	ch_base <- merge(m_ls[[1]], m_ls[[2]], by = "Tag..", all = TRUE)
	l_m <- length(m_ls) * 1
	for (i in 3:l_m) {
		# loop through the rest
		g <- as.data.frame(m_ls[[i]])
		ch_base <- merge(ch_base, g, by = "Tag..", all = TRUE)
	}
	ch_base <- as.matrix(ch_base[, -1]) # convert all NAs to 0
	for (i in seq_len(nrow(ch_base))) {
		for (t in seq_len(ncol(ch_base))) {
			if (is.na(ch_base[i, t])) {
				ch_base[i, t] <- 0
			}
			as.numeric(ch_base[i, t])
		}
	}
	colnames(ch_base) <- days
	return(ch_base)
}


mna_nimble <- function(site_run, return_mean = FALSE) {
	if (!grepl("Control", site_run)) {
		site_run <- paste(site_run, "Control")
	}
	ch <- suppressWarnings(ch_cary(site_run))
	ks <- known_states(ch)
	mna <- apply(ks, 2, sum)

	mice_obs <- ymd(colnames(ch)) # unique sampling days: mice

	# every day in mouse sequence
	mice_seq <- seq.Date(mice_obs[1], mice_obs[length(mice_obs)], by = 1)

	mna_all_days <- rep(NA, length(mice_seq))
	mna_count <- 1
	for (i in seq_along(mice_seq)) {
		if (mice_seq[i] %in% mice_obs) {
			mna_all_days[i] <- mna[mna_count]
			mna_count <- mna_count + 1
		} else {
			mna_all_days[i] <- mna[mna_count]
		}
	}

	# tick data
	dat <- read.csv("data/tick_cleaned") # tick data
	tick <- dat[, c("Grid", "DATE", "n_larvae", "n_nymphs", "n_adults")]
	tick <- subset(tick, Grid == site_run)
	tick$DATE <- as.Date(tick$DATE)

	# match tick dates to mice dates
	start_tick <- which(tick$DATE[1] == mice_seq)
	end_tick <- which(tick$DATE[length(tick$DATE)] == mice_seq)

	# index mice estimates - convert sd to prec
	mna_for_nimble <- mna_all_days[start_tick:end_tick]

	# center and scale
	mna_scaled <- scale(mna_for_nimble)

	if (return_mean) {
		return(
			list(
				mna = mna_scaled,
				mean = mean(mna_for_nimble),
				sd = sd(mna_for_nimble)
			)
		)
	} else {
		return(mna_scaled)
	}
}


mna_hindcast <- function(site_run) {
	df_mice <- read_csv("data/Mice2006_2021.csv")

	ch <- df_mice |>
		filter(Grid == site_run) |>
		select(-Grid) |>
		distinct() |>
		arrange(Date) |>
		pivot_wider(
			id_cols = Tag,
			names_from = Date,
			values_from = capture,
			values_fill = 0
		) |>
		select(-Tag)

	ks <- known_states(ch)
	mna <- apply(ks, 2, sum)
	mice_obs <- ymd(colnames(ch)) # unique sampling days: mice

	# every day in mouse sequence
	mice_seq <- seq.Date(mice_obs[1], mice_obs[length(mice_obs)], by = 1)

	mna_all_days <- rep(NA, length(mice_seq))
	mna_count <- 1
	for (i in seq_along(mice_seq)) {
		if (mice_seq[i] %in% mice_obs) {
			mna_all_days[i] <- mna[mna_count]
			mna_count <- mna_count + 1
		} else {
			mna_all_days[i] <- mna[mna_count]
		}
	}

	# historical mna
	mna_hist <- mna_nimble(site_run, return_mean = TRUE)

	# center and scale
	mna_scaled <- (mna_all_days - mna_hist$mean) / mna_hist$sd

	tibble(
		mna = mna_scaled,
		time = mice_seq
	)
}
