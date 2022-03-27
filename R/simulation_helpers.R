######################################################################
#
# Helpers functions for numerical simulations.
#
######################################################################

# generate dxd AR(1) covariance matrix with parameter rho

#' Generate dxd AR(1) covariance matrix with parameter rho
#'
#' @param rho The autocorrelation parameter
#' @param d The dimension of the matrix
#'
#' @return A dxd covariance matrix.
#' @export
generate_cov_ar1 = function(rho, d){
  sig <- matrix(0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (j in 1:d) {
      sig[i,j] <- rho**(abs(i-j))
    }
  }
  sig
}

#' A fast way to generate data from multivariate Gaussian distribution
#'
#' @param mean The mean vector of MVN.
#' @param covariance The covariance matrix of MVN.
#' @param B The number of replicate
#' @param n The number of sample
#'
#' @return A data matrix generated via Cholesky
#' @export
fast_data_generate <- function(mean, covariance, B, n){
  A <- Matrix::chol(covariance)
  d <- ncol(covariance)
  indep_data <- t(matrix(rnorm(B*n*d), nrow = B*n, ncol = d)) + mean
  trans_data <- t(indep_data) %*% A
  return(trans_data)
}