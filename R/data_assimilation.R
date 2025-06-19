library(nimble)

model_code <- nimbleCode({
	### priors
	phi.l.mu ~ dnorm(pr_phi_l[1], tau = pr_phi_l[2]) # larvae survival
	phi.n.mu ~ dnorm(pr_phi_n[1], tau = pr_phi_n[2]) # nymph survival
	phi.a.mu ~ dnorm(pr_phi_a[1], tau = pr_phi_a[2]) # adult survival
	theta.ln ~ dnorm(pr_theta_l2n[1], tau = pr_theta_l2n[2]) # larvae -> dormant nymph daily transition
	theta.na ~ dnorm(pr_theta_n2a[1], tau = pr_theta_n2a[2]) # larvae -> questing nymph daily transition

	for (j in 1:n_beta) {
		beta[j] ~ dnorm(pr_beta[j, 1], tau = pr_beta[j, 2])
	}

	for (i in 1:ns) {
		sig[i] ~ dinvgamma(pr_sig[i, 1], pr_sig[i, 2])
	}

	### precision priors
	OMEGA[1, 1] <- sig[1]
	OMEGA[1, 2] <- 0
	OMEGA[1, 3] <- 0
	OMEGA[1, 4] <- 0
	OMEGA[2, 1] <- 0
	OMEGA[2, 2] <- sig[2]
	OMEGA[2, 3] <- 0
	OMEGA[2, 4] <- 0
	OMEGA[3, 1] <- 0
	OMEGA[3, 2] <- 0
	OMEGA[3, 3] <- sig[3]
	OMEGA[3, 4] <- 0
	OMEGA[4, 1] <- 0
	OMEGA[4, 2] <- 0
	OMEGA[4, 3] <- 0
	OMEGA[4, 4] <- sig[4]

	# Cholesky decomposition
	omega_chol[1:ns, 1:ns] <- chol(OMEGA[1:ns, 1:ns])

	### first latent process
	for (i in 1:4) {
		x[i, 1] ~ dnorm(IC[i, 1], tau = IC[i, 2])
		px[i, 1] <- x[i, 1]
	}

	# if there are missing weather from Cary
	if (mu_f_missing) {
		for (i in 1:n_missing) {
			muf[t_missing[i], var_missing[i]] ~ dnorm(0, tau = 1)
		}
	}

	### define parameters, loop over every day in time series
	for (t in 1:horizon) {
		if (use_mice) {
			logit(l2n[t]) <- theta.ln + beta[13] * mice[t]
			logit(n2a[t]) <- theta.na + beta[14] * mice[t]
		} else {
			logit(l2n[t]) <- theta.ln
			logit(n2a[t]) <- theta.na
		}

		theta_n2a[t] <- if_else_nimble(
			(gdd[t] <= 1000) | (gdd[t] >= 2500),
			n2a[t],
			0
		)
		lambda[t] <- if_else_nimble(
			(gdd[t] >= 1400) & (gdd[t] <= 2500),
			repro_mu,
			0
		)
		l2n_quest[t] <- if_else_nimble((gdd[t] >= 400) & (gdd[t] <= 2500), 1, 0)

		if (nmme_cens) {
			for (i in 1:K) {
				y.censored[i, 1:J, t] ~
					dwtmnorm(mean = muf[t, 1:J], prec = pf[1:J, 1:J, t], wt = wts[i])

				y.ind[i, 2, t] ~
					dconstraint(
						interval[2, 1] <= y.censored[i, 2, t] &
							y.censored[i, 2, t] <= interval[2, 2]
					)
				y.ind[i, 3, t] ~
					dconstraint(
						interval[3, 1] <= y.censored[i, 3, t] &
							y.censored[i, 3, t] <= interval[3, 2]
					)
				y.ind[i, 4, t] ~ dconstraint(interval[4, 1] <= y.censored[i, 4, t])
			}

			muf[t, 1:J] ~ dmnorm(mean = mu_0[1:J], prec = Sigma_0[1:J, 1:J])
			pf[1:J, 1:J, t] ~ dwish(S = lambda_0[1:J, 1:J], df = nu_0)
			gdd[t] ~ dnorm(cgdd.mu[t], tau = cgdd.tau[t])
		}

		# daily survival model
		logit(phi_l[t]) <- phi.l.mu +
			beta[1] * muf[t, 1] +
			beta[2] * muf[t, 2] +
			beta[3] * muf[t, 3] +
			beta[4] * muf[t, 4]

		logit(phi_n[t]) <- phi.n.mu +
			beta[5] * muf[t, 1] +
			beta[6] * muf[t, 2] +
			beta[7] * muf[t, 3] +
			beta[8] * muf[t, 4]

		logit(phi_a[t]) <- phi.a.mu +
			beta[9] * muf[t, 1] +
			beta[10] * muf[t, 2] +
			beta[11] * muf[t, 3] +
			beta[12] * muf[t, 4]

		A[1, 1, t] <- phi_l[t] * (1 - l2n[t])
		A[1, 2, t] <- 0
		A[1, 3, t] <- 0
		A[1, 4, t] <- lambda[t]
		A[2, 1, t] <- phi_l[t] * l2n[t]
		A[2, 2, t] <- 1 - l2n_quest[t]
		A[2, 3, t] <- 0
		A[2, 4, t] <- 0
		A[3, 1, t] <- 0
		A[3, 2, t] <- l2n_quest[t]
		A[3, 3, t] <- phi_n[t] * (1 - theta_n2a[t])
		A[3, 4, t] <- 0
		A[4, 1, t] <- 0
		A[4, 2, t] <- 0
		A[4, 3, t] <- phi_n[t] * theta_n2a[t]
		A[4, 4, t] <- phi_a[t]

		### Data Model ###
		y[1, t] ~ dpois(x[1, t])
		y[3, t] ~ dpois(x[3, t])
		y[4, t] ~ dpois(x[4, t])
	}

	for (t in 2:horizon) {
		# expected number questing
		ex[1:ns, t] <- A[1:ns, 1:ns, t - 1] %*% x[1:ns, t - 1]

		# process error
		px[1:ns, t] ~
			dmnorm(
				mean = ex[1:ns, t],
				cholesky = omega_chol[1:ns, 1:ns],
				prec_param = 0
			)

		x[1, t] <- px[1, t]
		x[2, t] <- max(px[2, t], 0)
		x[3, t] <- px[3, t]
		x[4, t] <- px[4, t]
	}
})
