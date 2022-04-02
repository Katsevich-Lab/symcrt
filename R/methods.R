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
#' @param test_hyperparams Additional method hyperparameters (currently not used)
#'
#' @return A data frame with columns "parameter," "target," "value",
#' with two rows, one for the test statistic, and one for the p-value.
#'
#' @export
MX2_F_test <- function(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams) {
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # fit conditional mean of X given Z
  if (X_on_Z_reg$mean_method_type == "oracle") {
    E_X_given_Z <- data$E_X_given_Z_oracle
    if (is.null(E_X_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)
  }

  # fit conditional variance of X given Z
  if (X_on_Z_reg$var_method_type == "oracle") {
    Var_X_given_Z <- data$Var_X_given_Z_oracle
    if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  } else {
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_reg$var_method_type)
  }

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)

  # define the test statistic
  X_residuals <- X - E_X_given_Z
  Y_residuals <- Y - E_Y_given_Z
  S_hat <- sqrt(mean(Var_X_given_Z * Y_residuals^2))
  test_statistic <- 1 / (sqrt(n) * S_hat) * sum(X_residuals * Y_residuals)

  # define the p-value (for now, define as two-sided)
  p_value <- 2 * stats::pnorm(abs(test_statistic), lower.tail = FALSE)

  # output the results
  data.frame(
    parameter = c("test_statistic", "p_value"),
    target = "conditional_independence",
    value = c(test_statistic, p_value)
  )
}

#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM_f} is a function carrying out the dCRT.
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#' @param test_hyperparams Additional method hyperparameters (currently not used)

#' @return A data frame with columns "parameter," "target," "value",
#' with two rows, one for the test statistic, and one for the p-value.
#'
#' @export
GCM <- function(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams) {
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
  R <- X_residuals * Y_residuals
  test_statistic <- sqrt(n) * mean(R) / sqrt(mean(R^2) - (mean(R))^2)

  # define the p-value (for now, define as two-sided)
  p_value <- 2 * stats::pnorm(abs(test_statistic), lower.tail = FALSE)

  # output the results
  data.frame(
    parameter = c("test_statistic", "p_value"),
    target = "conditional_independence",
    value = c(test_statistic, p_value)
  )
}



#' A function that outputs the bias term in GCM test
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#' @param test_hyperparams Additional method hyperparameters (currently not used)
#'
#' @return A data frame with columns "parameter," "target," "value",
#' with three rows, are respectively cross_bias, X_Z_bias, Y_Z_bias
#' @export

GCM_debug <- function(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams) {
  # extract X, Y, Z, cond_mean_X_Z, cond_mean_Y_Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  if (is.null(data$cond_mean_X_Z)){
    stop("GCM_debug should know the oracle conditional mean of X given Z!")
  }else{
    cond_mean_X_Z <- data$cond_mean_X_Z
  }
  if (is.null(data$cond_mean_Y_Z)){
    stop("GCM_debug should know the oracle conditional mean of Y given Z!")
  }else{
    cond_mean_Y_Z <- data$cond_mean_Y_Z
  }
  n <- nrow(Z)

  # fit conditional mean of X given Z
  E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg) # hat{E(X|Z)}

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg) # hat{E(Y|Z)}

  # compute the oracle residuals
  X_oracle_residuals <- X - cond_mean_X_Z # X - E(X|Z)
  Y_oracle_residuals <- Y - cond_mean_Y_Z # Y - E(Y|Z)
  cond_mean_diff_X_Z <- cond_mean_X_Z - E_X_given_Z # E(X|Z) - hat{E(X|Z)}
  cond_mean_diff_Y_Z <- cond_mean_Y_Z - E_Y_given_Z # E(Y|Z) - hat{E(Y|Z)}

  # compute three bias terms
  cross_bias <- sqrt(n)*mean(cond_mean_diff_X_Z*cond_mean_diff_Y_Z)
  X_Z_bias <- sqrt(n)*mean(Y_oracle_residuals*cond_mean_diff_X_Z)
  Y_Z_bias <- sqrt(n)*mean(X_oracle_residuals*cond_mean_diff_Y_Z)

  # output the results
  data.frame(parameter = c("cross_bias", "X_Z_bias", "Y_Z_bias"),
             target = "bias",
             value = c(cross_bias, X_Z_bias, Y_Z_bias))
}




#' The distilled CRT test
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#' @param test_hyperparams Additional test hyperparameters

#' @return A data frame with columns "parameter," "target," "value" with p-value.
#' @export
dCRT <- function(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams) {
  # TODO: Write a test for this function as follows: Give it a Gaussian resampling
  # distribution, and then compare the output to the MX(2) F-test
  test_hyperparams <- set_default_test_hyperparams("dCRT", test_hyperparams)
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # fit conditional mean of X given Z
  if (X_on_Z_reg$mean_method_type == "oracle") {
    E_X_given_Z <- data$E_X_given_Z_oracle
    if (is.null(E_X_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)
  }

  # fit conditional variance of X given Z
  if (X_on_Z_reg$var_method_type == "oracle") {
    Var_X_given_Z <- data$Var_X_given_Z_oracle
    if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  } else {
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_reg$var_method_type)
  }

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)

  # define the test statistic
  X_residuals <- c(X - E_X_given_Z)
  Y_residuals <- c(Y - E_Y_given_Z)
  S_hat <- sqrt(mean(Var_X_given_Z * Y_residuals^2))
  test_statistic <- 1 / (sqrt(n) * S_hat) * sum(X_residuals * Y_residuals)

  # resample matrix from the specified distribution
  resample_matrix <- resample_dCRT(conditional_mean = E_X_given_Z,
                                   conditional_variance = Var_X_given_Z,
                                   no_resample = test_hyperparams$no_resample,
                                   resample_dist = test_hyperparams$resample_family)

  # compute the residuals and variance vector for each resample
  resample_X_residuals <- resample_matrix - E_X_given_Z
  residual_star <- apply(resample_X_residuals * Y_residuals, 2, sum)

  # compute the resample test statistic and quantile
  resample_test_statistic <- 1 / (sqrt(n) * S_hat) * residual_star
  no_exceed <- length(which(abs(resample_test_statistic) >= abs(test_statistic)))

  # compute the p-value (two sided)
  p_value <- (no_exceed + 1) / (test_hyperparams$no_resample + 1)

  # output the results
  data.frame(
    parameter = c("p_value"),
    target = "conditional_independence",
    value = c(p_value)
  )
}
