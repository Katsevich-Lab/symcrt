######################################################################
#
# Helpers functions for numerical simulations.
#
######################################################################

#' Generate dxd AR(1) covariance matrix with parameter rho
#'
#' @param rho The autocorrelation parameter
#' @param d The dimension of the matrix
#'
#' @return A dxd covariance matrix.
#' @export
generate_cov_ar1 = function(rho, d){
  stats::toeplitz(rho^(0:(d-1)))
}

#' A fast way to generate data from multivariate Gaussian distribution
#'
#' @param mean The mean vector of MVN.
#' @param covariance The covariance matrix of MVN.
#' @param num_samples The number of samples to generate
#'
#' @return A matrix with \code{sum_samples} rows and \code{d = ncol(covariance)}
#' columns, such that each row is drawn i.i.d. from the multivariate normal
#' specified by the arguments.
#'
#' @export
fast_generate_mvn <- function(mean, covariance, num_samples){
  A <- Matrix::chol(covariance)
  d <- ncol(covariance)
  indep_data <- t(matrix(stats::rnorm(num_samples*d), nrow = num_samples, ncol = d)) + mean
  trans_data <- t(indep_data) %*% A
  return(trans_data)
}
