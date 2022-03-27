######################################################################
#
# Conditional independence testing methods used for simulations for
# symcrt project.
#
######################################################################

#' The MX(2) F-test.
#'
#' \code{MX2_F_test_f} is a function carrying out the MX(2) F-test.
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#' @param X_on_Z_var The method to obtain the conditional variance of X|Z.
#'
#' @return A data frame with columns "parameter," "target," "value",
#' with two rows, one for the test statistic, and one for the p-value.
#'
#' @export
MX2_F_test_f <- function(data, X_on_Z_reg, Y_on_Z_reg, X_on_Z_var) {
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # fit conditional mean of X given Z
  if(X_on_Z_reg == "oracle"){
    E_X_given_Z <- data$cond_mean
    if(is.null(E_X_given_Z)) stop("Must specify cond_mean if X_on_Z_reg is oracle")
  } else{
    E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)
  }

  # fit conditional variance of X given Z
  if(X_on_Z_var == "oracle"){
    Var_X_given_Z <- data$cond_var
    if(is.null(Var_X_given_Z)) stop("Must specify cond_var if X_on_Z_var is oracle")
  } else{
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_var)
  }

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)

  # define the test statistic
  X_residuals <- X - E_X_given_Z
  Y_residuals <- Y - E_Y_given_Z
  S_hat <- sqrt(mean(Var_X_given_Z*Y_residuals^2))
  test_statistic <- 1/(sqrt(n)*S_hat)*sum(X_residuals*Y_residuals)

  # define the p-value (for now, define as two-sided)
  p_value <- 2*stats::pnorm(abs(test_statistic), lower.tail = FALSE)

  # output the results
  data.frame(parameter = c("test_statistic", "p_value"),
             target = "conditional_independence",
             value = c(test_statistic, p_value))
}

#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM_f} is a function carrying out the dCRT.
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#'
#' @return A data frame with columns "parameter," "target," "value",
#' with two rows, one for the test statistic, and one for the p-value.
#'
#' @export
GCM_f <- function(data, X_on_Z_reg, Y_on_Z_reg) {
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # fit conditional mean of X given Z
  E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)

  # define the test statistic
  X_residuals <- X - E_X_given_Z
  Y_residuals <- Y - E_Y_given_Z
  R = X_residuals*Y_residuals
  test_statistic <- sqrt(n)*mean(R)/sqrt(mean(R^2)-(mean(R))^2)

  # define the p-value (for now, define as two-sided)
  p_value <- 2*stats::pnorm(abs(test_statistic), lower.tail = FALSE)

  # output the results
  data.frame(parameter = c("test_statistic", "p_value"),
             target = "conditional_independence",
             value = c(test_statistic, p_value))
}


#' The distilled CRT test
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#' @param X_on_Z_var The method to obtain the conditional variance of X|Z.
#' @param resample_dist The resampling distribution for CRT
#' @param no_resample The number of resampling
#'
#' @return A data frame with columns "parameter," "target," "value" with p-value.
#' @export
dCRT_f <- function(data,
                   X_on_Z_reg,
                   Y_on_Z_reg,
                   X_on_Z_var,
                   resample_dist,
                   no_resample = 2000) {
  # TODO: Write a test for this function as follows: Give it a Gaussian resampling
  # distribution, and then compare the output to the MX(2) F-test

  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # fit conditional mean of X given Z
  if(X_on_Z_reg == "oracle"){
    E_X_given_Z <- data$cond_mean
    if(is.null(E_X_given_Z)) stop("Must specify cond_mean if X_on_Z_reg is oracle")
  } else{
    E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)
  }

  # fit conditional variance of X given Z
  if(X_on_Z_var == "oracle"){
    Var_X_given_Z <- data$cond_var
    if(is.null(Var_X_given_Z)) stop("Must specify cond_var if X_on_Z_var is oracle")
  } else{
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_var)
  }

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)

  # define the test statistic
  X_residuals <- X - E_X_given_Z
  Y_residuals <- Y - E_Y_given_Z
  S_hat <- sqrt(mean(Var_X_given_Z*Y_residuals^2))
  test_statistic <- 1/(sqrt(n)*S_hat)*sum(X_residuals*Y_residuals)

  # resample matrix from the specified distribution
  resample_matrix <- resample_dCRT(E_X_given_Z, Var_X_given_Z, no_resample, resample_dist)

  # compute the residuals and variance vector for each resample
  resample_X_residuals <- resample_matrix - E_X_given_Z
  residual_star <- apply(resample_X_residuals*Y_residuals, 2, sum)

  # compute the resample test statistic and quantile
  resample_test_statistic <- 1/(sqrt(n)*S_hat)*residual_star
  no_exceed <- length(which(abs(resample_test_statistic) >= abs(test_statistic)))

  # compute the p-value (two sided)
  p_value <- (no_exceed+1)/(no_resample+1)

  # output the results
  data.frame(parameter = c("p_value"),
             target = "conditional_independence",
             value = c(p_value))
}
