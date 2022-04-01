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
generate_cov_ar1 <- function(rho, d) {
  stats::toeplitz(rho^(0:(d - 1)))
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
fast_generate_mvn <- function(mean, covariance, num_samples) {
  A <- Matrix::chol(covariance)
  d <- ncol(covariance)
  indep_data <- t(matrix(stats::rnorm(num_samples * d), nrow = num_samples, ncol = d)) + mean
  trans_data <- t(indep_data) %*% A
  return(trans_data)
}

#' Generate response data based on a GLM model
#'
#' @param design_matrix The design matrix (n x p)
#' @param coefficients Coefficient vector (vector of length p)
#' @param family One of \code{gaussian, binomial, poisson}
#'
#' @return A response vector of length n
#' @export
generate_glm_response_data <- function(design_matrix, coefficients, family){
  # obtain the natural parameter for each observation
  eta = design_matrix %*% coefficients
  n = length(eta)

  # generate mean function and samples based on the family
  switch(family,
         gaussian = {
           mu <- eta
           response <- stats::rnorm(n = n, mean = mu, sd = 1)
         },
         binomial = {
           mu <- exp(eta)/(1+exp(eta))
           response <- stats::rbinom(n = n, size = 1, prob = mu)
         },
         poisson = {
           mu <- exp(eta)
           response <- stats::rpois(n = n, lambda = mu)
         }
  )

  # return the response vector
  response
}

generate_method_list <- function(methods_df){
  num_methods <- nrow(methods_df)
  for(method_idx in 1:num_methods){
    test_type <- methods_df$test_type[method_idx]
    X_on_Z_reg <- methods_df$X_on_Z_reg[method_idx]
    Y_on_Z_reg <- methods_df$Y_on_Z_reg[method_idx]
    test_hyperparams <- methods_df$test_hyperparams[method_idx]



    switch(test_type,
           squared_residual = {
             (response - conditional_mean)^2
           },
           homoskedastic = {
             n <- length(response)
             rep(mean((response - conditional_mean)^2), n)
           },
           {
             stop("Invalid specification of conditional variance estimation method.")
           }
    )
  }

}
