run_hindcast_nimble <- function(
  cl,
  code,
  data,
  constants,
  inits,
  n_iter,
  thin,
  nmme_cens,
  mu_f_missing,
  use_mice
) {
  # Run hindcast using nimble
  # cl: cluster object
  # code: nimble model object
  # data: list of data objects
  # constants: list of constants
  # inits: function to generate initial values
  # n_iter: number of iterations for MCMC
  # thin: thinning interval for MCMC
  # nmme_cens: logical, whether to include censored weather nodes
  # mu_f_missing: logical, whether to handle missing weather data with micee
  # use_mice: logical, whether to use micee for missing data imputation

  require(parallel)
  require(nimble)
  require(coda)

  source("R/functions_nimble.R")

  n_cores <- length(cl) # number of cores used

  export_vec <- c(
    "code",
    "constants",
    "data",
    "n_iter",
    "thin",
    "n_cores",
    "if_else_nimble",
    "dwtmnorm",
    "rwtmnorm",
    "nmme_cens",
    "mu_f_missing",
    "use_mice"
  )

  clusterExport(cl, export_vec, envir = environment())

  # export inits to clusters
  for (j in seq_len(n_cores)) {
    set.seed(j)
    init <- inits()
    clusterExport(cl[j], "init", envir = environment())
  }

  message("Running mcmc...")
  # sample on each cluster
  out <- clusterEvalQ(cl, {
    library(nimble)
    library(coda)

    nimbleOptions(MCMCusePosteriorPredictiveSampler = FALSE)

    registerDistributions(list(
      dwtmnorm = list(
        BUGSdist = "dwtmnorm(mean, prec, wt)",
        types = c(
          "value = double(1)",
          "mean = double(1)",
          "prec = double(2)",
          "wt = double(0)"
        )
      )
    ))

    model <- nimbleModel(
      code = code,
      constants = constants,
      data = data,
      inits = init
    )

    c_model <- compileNimble(model)

    monitor <- c(
      "phi.l.mu",
      "phi.n.mu",
      "phi.a.mu",
      "theta.ln",
      "theta.na",
      "sig",
      "beta",
      "x"
    )

    mcmc_conf <- configureMCMC(c_model, monitors = monitor, thin = thin)

    # add appropriate samples for y.censored nodes
    if (nmme_cens) {
      message("--- Add samplers for y.censored nodes ---")
      mcmc_conf$addMonitors(c("muf", "pf"))
      for (j in seq_along(data$muf)) {
        for (n in seq_len(nrow(data$y.ind))) {
          for (t in seq_len(constants$horizon)) {
            if (data$y.ind[n, j, t] == 0) {
              node <- paste0("y.censored[", n, ", ", j, ", ", t, "]")
              mcmc_conf$addMonitors(node)
              mcmc_conf$addSampler(target = node, type = "RW")
            }
          }
        }
      }
    } else {
      if (mu_f_missing) {
        message("--- Remove samplers for observed weather nodes ---")
        for (i in seq_len(nrow(data$muf))) {
          for (j in seq_len(ncol(data$muf))) {
            if (!is.na(data$muf[i, j])) {
              node <- paste0("muf[", i, ", ", j, "]")
              mcmc_conf$removeSamplers(node)
            }
          }
        }
      }
    }

    for (i in seq_len(nrow(data$y))) {
      for (j in seq_len(ncol(data$y))) {
        if (j > 1) {
          node <- paste0("px[", i, ", ", j, "]")
          mcmc_conf$addSampler(node, "RW")
        }
      }
    }

    mcmc_conf$printSamplers(byType = TRUE)

    r_mcmc <- buildMCMC(mcmc_conf)
    c_mcmc <- compileNimble(r_mcmc)
    c_mcmc$run(niter = n_iter, nburnin = 0.3 * n_iter)

    return(as.mcmc(as.matrix(c_mcmc$mvSamples)))
  })

  # x1 <- as.mcmc(as.matrix(c_mcmc$mvSamples))
  # x2 <- x1[, grep("x", colnames(x1))]
  # summary(x2[, 1:20])
  # which(as.vector(apply(x2, 2, min)) < 0)

  message("MCMC run complete.")
  as.mcmc.list(out)
}
