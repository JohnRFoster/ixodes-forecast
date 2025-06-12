library(nimble)

#' random weighted multivariate normal
#' @param n sample size
#' @param mean mean
#' @param prec precision
#' @param wt weight
#' @export
rwtmnorm <- nimbleFunction(
	run = function(
		n = integer(0),
		mean = double(1),
		prec = double(2),
		wt = double(0)
	) {
		returnType(double(1))
		if (n != 1) {
			nimPrint("rwtmnorm only allows n = 1; using n = 1.")
		}
		prob <-
			rmnorm_chol(n = 1, mean, chol(prec), prec_param = TRUE) * wt
		return(prob)
	}
)

#' weighted multivariate normal density
#' @param x random variable
#' @param mean mean
#' @param prec precision
#' @param wt weight
#' @param log log
#' @export
dwtmnorm <- nimbleFunction(
	run = function(
		x = double(1),
		mean = double(1),
		prec = double(2),
		wt = double(0),
		log = integer(0, default = 0)
	) {
		returnType(double(0))

		log_prob <-
			dmnorm_chol(
				x = x,
				mean = mean,
				cholesky = chol(prec),
				prec_param = TRUE,
				log = TRUE
			) *
			wt

		if (log) {
			return((log_prob))
		} else {
			return((exp(log_prob)))
		}
	}
)

if_else_nimble <- nimbleFunction(
	run = function(
		condition = integer(0),
		value_if = double(0),
		values_else = double(0)
	) {
		returnType(double(0))
		if (condition == TRUE) {
			return(value_if)
		} else {
			return(values_else)
		}
	}
)
