######################################################################
#
# Conditional independence testing methods used for simulations.
#
# NOTE: The design matrix argument is not actually used by these functions
#       and should be removed when Tim gets back to us.
#
######################################################################

##################################################################################
#
# The MX(2) F-test.
#
# Arguments:
#   - data: A named list with fields X, Y, Z
#   - X_on_Z_reg: The regression method to apply for X|Z
#   - Y_on_Z_reg: The regression method to apply for Y|Z
#   - X_on_Z_var: The method to obtain the conditional variance of X|Z.
#
# Output: A data frame with columns "parameter," "target," "value",
#         with two rows, one for the test statistic, and one for the p-value.
#
##################################################################################


#' MX2_F_test_f
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
#'
#' @examples
MX2_F_test_f <- function(data, X_on_Z_reg, Y_on_Z_reg, X_on_Z_var) {
  source("method_helpers.R")

  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # fit conditional mean of X given Z
  E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)

  # fit conditional variance of X given Z
  Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_var)

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

# Lasso regression & logistic-regression
## known variance for X|Z
# Note: y is the response vector of length 2n, the first column is the response of ols and the second is glm
# the output is a data frame with columns "parameter," "target," "value"
MX2_Gaussf <- function(y, design_matrix) {
  n <- length(y)
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_residual <- y[,d+1] - design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2])
  XZ_residual <- y[,d+2] - design_matrix %*% (glmnet::coef.glmnet(object = XZ_lasso, "lambda.min")[,1][-1])
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
  prod_res <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"),
             value = c(stats::pnorm(prod_res), prod_res))
}

# OLS regression & logistic-regression
## known variance for X|Z
# Note: y is the response vector of length 2n, the first column is the response of ols and the second is glm
# the output is a data frame with columns "parameter," "target," "value"
dCRT_Gaussf <- function(y, design_matrix) {
  n <- length(y[,1])
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2])
  XZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = XZ_lasso, "lambda.min")[,1][-1])
  YZ_residual <- y[,d+1] - YZ_fitted
  XZ_residual <- y[,d+2] - XZ_fitted
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
  test_stat <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  ## resampling procedure
  dcrt <- rep(0, 1000)
  for (m in 1:1000) {
    set.seed(m)
    result_YZ <- y[, d+1] - YZ_fitted
    result_XZ <- MASS::mvrnorm(1, mu = XZ_fitted, Sigma = diag(n)) - XZ_fitted
    variance_YZ <- result_YZ**2
    variance_XZ <- result_XZ**2
    sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
    re_test_stat <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
    if (re_test_stat >= test_stat){
      dcrt[m] <- 1
    }
  }
  p_value <- (sum(dcrt)+1)/1001
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"), value = c(p_value, re_test_stat))
}


GCM_Gaussf <- function(y, design_matrix) {
  n <- length(y)
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_residual <- y[,d+1] - design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2])
  XZ_residual <- y[,d+2] - design_matrix %*% (glmnet::coef.glmnet(object = XZ_lasso, "lambda.min")[,1][-1])
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  ## only difference with MX(2)
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n-(sum(result_XZ*result_YZ)/n)**2)
  prod_res <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"),
             value = c(stats::pnorm(prod_res), prod_res))
}

MX2_logisf <- function(y, design_matrix) {
  n <- length(y)
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_residual <- y[,d+1] - design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_logis <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2], family = "binomial")
  XZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = XZ_logis, "lambda.min")[,1][-1])
  XZ_residual <- y[,d+2] - exp(XZ_fitted)/(1+exp(XZ_fitted))
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
  prod_res <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"),
             value = c(stats::pnorm(prod_res), prod_res))
}

# OLS regression & logistic-regression
## known variance for X|Z
# Note: y is the response vector of length 2n, the first column is the response of ols and the second is glm
# the output is a data frame with columns "parameter," "target," "value"
dCRT_logisf <- function(y, design_matrix) {
  n <- length(y[,1])
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_logis <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2], family = "binomial")
  XZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = XZ_logis, "lambda.min")[,1][-1])
  YZ_residual <- y[,d+1] - YZ_fitted
  XZ_residual <- y[,d+2] - exp(XZ_fitted)/(1+exp(XZ_fitted))
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
  test_stat <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  ## resampling procedure
  dcrt <- rep(0, 10000)
  for (m in 1:10000) {
    set.seed(m)
    result_YZ <- y[, d+1] - YZ_fitted
    prob_fitted <- exp(XZ_fitted)/(1+exp(XZ_fitted))
    result_XZ <- stats::rbinom(n, 1, prob = prob_fitted) - exp(XZ_fitted)/(1+exp(XZ_fitted))
    variance_YZ <- result_YZ**2
    variance_XZ <- result_XZ**2
    sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
    re_test_stat <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
    if (re_test_stat >= test_stat){
      dcrt[m] <- 1
    }
  }
  p_value <- (sum(dcrt)+1)/10001
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"), value = c(p_value, re_test_stat))
}


GCM_logisf <- function(y, design_matrix) {
  n <- length(y)
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_residual <- y[,d+1] - design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_logis <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2], family = "binomial")
  XZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = XZ_logis, "lambda.min")[,1][-1])
  XZ_residual <- y[,d+2] - exp(XZ_fitted)/(1+exp(XZ_fitted))
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  ## only difference with MX(2)
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n-(sum(result_XZ*result_YZ)/n)**2)
  prod_res <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"),
             value = c(stats::pnorm(prod_res), prod_res))
}


dCRT_lapf <- function(y, design_matrix) {
  n <- length(y[,1])
  d <- dim(design_matrix)[2]
  design_matrix <- y[, 1:d]
  YZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+1])
  YZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = YZ_lasso, "lambda.min")[,1][-1])
  XZ_lasso <- glmnet::cv.glmnet(x = design_matrix, y = y[, d+2])
  XZ_fitted <- design_matrix %*% (glmnet::coef.glmnet(object = XZ_lasso, "lambda.min")[,1][-1])
  YZ_residual <- y[,d+1] - YZ_fitted
  XZ_residual <- y[,d+2] - XZ_fitted
  result_YZ <- YZ_residual
  result_XZ <- XZ_residual
  variance_YZ <- result_YZ**2
  variance_XZ <- result_XZ**2
  sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
  test_stat <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
  ## resampling procedure
  dcrt <- rep(0, 1000)
  for (m in 1:1000) {
    set.seed(m)
    result_YZ <- y[, d+1] - YZ_fitted
    result_XZ <- ExtDist::rLaplace(n, mu = XZ_fitted, b = 1) - XZ_fitted
    variance_YZ <- result_YZ**2
    variance_XZ <- result_XZ**2
    sd_residual <- sqrt(sum(variance_XZ*variance_YZ)/n)
    re_test_stat <- sum(result_XZ*result_YZ)/(sqrt(n)*sd_residual)
    if (re_test_stat >= test_stat){
      dcrt[m] <- 1
    }
  }
  p_value <- (sum(dcrt)+1)/1001
  data.frame(parameter = c("p_value", "statistic"), target = c("estimate"), value = c(p_value, re_test_stat))
}
