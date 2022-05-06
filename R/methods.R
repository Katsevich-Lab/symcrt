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
    E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)$conditional_mean
  }

  # fit conditional variance of X given Z
  if (X_on_Z_reg$var_method_type == "oracle") {
    Var_X_given_Z <- data$Var_X_given_Z_oracle
    if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  } else {
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_reg$var_method_type)
  }

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)$conditional_mean

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
  E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)$conditional_mean

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)$conditional_mean

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
  E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)$conditional_mean # hat{E(X|Z)}

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)$conditional_mean # hat{E(Y|Z)}

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
    E_X_given_Z <- fit_conditional_mean(X, Z, X_on_Z_reg)$conditional_mean
  }

  # fit conditional variance of X given Z
  if (X_on_Z_reg$var_method_type == "oracle") {
    Var_X_given_Z <- data$Var_X_given_Z_oracle
    if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  } else {
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_reg$var_method_type)
  }

  # fit conditional mean of Y given Z
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)$conditional_mean

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




#' Maxway CRT proposed in Li & Liu 2022.
#'
#' @param data A named list with fields X, Y, Z.
#' @param X_on_Z_reg The regression method to apply for X|Z.
#' @param Y_on_Z_reg The regression method to apply for Y|Z.
#' @param test_hyperparams Additional test hyperparameters

#' @return A data frame with columns "parameter," "target," "value" with p-value.
#' @export
MaxwayCRT <- function(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams) {
  # TODO: Write a test for this function as follows: Give it a Gaussian resampling
  # distribution, and then compare the output to Molei's implementation
  test_hyperparams <- set_default_test_hyperparams("dCRT", test_hyperparams)
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)
  
  # extract the proportion of the number of unlabelled data and sample 1:n
  no_unlabel <- round(n*test_hyperparams$unlabel_prop)
  index_unlabel <- sample(1:n, no_unlabel)
  index_rest <- setdiff(1:n, index_unlabel)
  
  # extract the proportion of the number of holdout data and sample 1:n
  no_holdout <- round(n*test_hyperparams$holdout_prop)
  index_holdout <- sample(index_rest, no_holdout)
  
  # extract the index that used for testing
  index_test <- setdiff(index_rest, index_holdout)
  
  # fit conditional mean of X given Z with unlablled data
  if (X_on_Z_reg$mean_method_type == "oracle") {
    E_X_given_Z_full <- data$E_X_given_Z_oracle
    E_X_given_Z <- E_X_given_Z_full[index_unlabel]
    if (is.null(E_X_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    fit_X_given_Z <- fit_conditional_mean(X[index_unlabel], 
                                          Z[index_unlabel,], 
                                          X_on_Z_reg)
    E_X_given_Z <- fit_X_given_Z$conditional_mean
    coef_X_given_Z <- fit_X_given_Z$coef_vec
  }
  
  # # fit conditional variance of X given Z
  # if (X_on_Z_reg$var_method_type == "oracle") {
  #   Var_X_given_Z <- data$Var_X_given_Z_oracle
  #   if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  # } else {
  #   Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, X_on_Z_reg$var_method_type)
  # }
  
  # fit conditional mean of Y given Z; construct g(Z)
  if(no_holdout == 0){
    fit_Y_given_Z <- fit_conditional_mean(Y[index_test], 
                                          Z[index_test, ],
                                          Y_on_Z_reg)
  }else{
    fit_Y_given_Z <- fit_conditional_mean(Y[index_holdout], 
                                          Z[index_holdout, ],
                                          Y_on_Z_reg)
  }
  E_Y_given_Z <- fit_Y_given_Z$conditional_mean
  coef_Y_given_Z <- fit_Y_given_Z$coef_vec
  
  
  # extract the active set part in g(Z)
  if(test_hyperparams$g_Z_support == "support of beta_hat"){
    act_set_Y_given_Z <- which(coef_Y_given_Z != 0)
  }else{
    act_set_Y_given_Z <- test_hyperparams$g_Z_support
  }
  
  # concatenate predictor and active Z together
  predictor_Y_Z_unlabel <- Z[index_unlabel, ]%*%coef_Y_given_Z
  predictor_Y_Z_test <- Z[index_test, ]%*%coef_Y_given_Z
  g_Z_unlabel <- cbind(predictor_Y_Z_unlabel, 
                       Z[index_unlabel, act_set_Y_given_Z])
  g_Z_test <- cbind(predictor_Y_Z_test,
                    Z[index_test, act_set_Y_given_Z])
  
  
  # fit the conditional expectation of Y|Z on label data
  Y_on_Z_reg_hyperparams <- Y_on_Z_reg$mean_method_hyperparams
  switch(Y_on_Z_reg_hyperparams$family,
         binomial = {
           E_Y_given_Z_test <- 1/exp(-Z[index_test, ]%*%coef_Y_given_Z)
         },
         gaussian = {
           E_Y_given_Z_test <- Z[index_test, ]%*%coef_Y_given_Z
         },
         {
           stop("The rsampling distribution of X given Z is invalid!")
         }
  )
  
  # learning the resampling distribution with unlabelled data
  switch(test_hyperparams$resample_family,
         binomial = {
           
           # construct h(Z)=Z%*%gamma_hat and fit conditional mean on label data
           h_Z_unlabel <- Z[index_unlabel, ]%*%coef_X_given_Z
           h_Z_test <- Z[index_test, ]%*%coef_X_given_Z
           g_h_Z_unlabel <- cbind(g_Z_unlabel, h_Z_unlabel)
           E_X_given_Z_test <- 1/exp(-Z[index_test, ]%*%coef_X_given_Z)
           
           # fit X on g and h
           fit_X_on_g_h <- fit_conditional_mean(X[index_unlabel], g_h_Z_unlabel, X_on_Z_reg)
           Var_X_on_g_h <- fit_X_on_g_h$conditional_mean*(1 - fit_X_on_g_h$conditional_mean)
           
           # resample matrix from the specified distribution
           g_h_Z_test <- cbind(g_Z_test, h_Z_test)
           resample_mean <- g_h_Z_test%*%(fit_X_on_g_h$coef_vec)
           resample_var <- Var_X_on_g_h
           resample_matrix <- resample_dCRT(conditional_mean = resample_mean,
                                            conditional_variance = resample_var,
                                            no_resample = test_hyperparams$no_resample,
                                            resample_dist = test_hyperparams$resample_family)
           resample_X_residuals <- resample_matrix - 1/(1+exp(-Z[index_test,]%*%coef_X_given_Z))
         },
         gaussian = {
           
           # compute the residual and fit conditional mean on label data
           r_unlabel <- X[index_unlabel] - Z[index_unlabel, ]%*%coef_X_given_Z
           r_label <- X[index_test] - Z[index_test, ]%*%coef_X_given_Z
           E_X_given_Z_test <- Z[index_test, ]%*%coef_X_given_Z
           
           # fit the residual and g(Z) with lasso
           fit_residual_g_Z <- fit_conditional_mean(r_unlabel, g_Z_unlabel, X_on_Z_reg)
           Var_residual_g_Z <- fit_conditional_variance(r_unlabel, 
                                                        g_Z_unlabel, 
                                                        fit_residual_g_Z$conditional_mean, 
                                                        X_on_Z_reg$var_method_type)
           
           # resample matrix from the specified distribution
           resample_mean <- g_Z_test%*%(fit_residual_g_Z$coef_vec)
           resample_var <- Var_residual_g_Z
           resample_matrix <- resample_dCRT(conditional_mean = resample_mean,
                                            conditional_variance = resample_var,
                                            no_resample = test_hyperparams$no_resample,
                                            resample_dist = test_hyperparams$resample_family)
           resample_X_residuals <- resample_matrix
         },
         {
           stop("The rsampling distribution of X given Z is invalid!")
         }
  )
  
  
  
  # define the test statistic
  X_residuals <- c(X[index_test] - E_X_given_Z_test)
  Y_residuals <- c(Y[index_test] - E_Y_given_Z_test)
  test_statistic <- 1 / (sqrt(n)) * sum(X_residuals * Y_residuals)
  
  
  # compute the residuals and variance vector for each resample
  residual_star <- apply(resample_X_residuals * Y_residuals, 2, sum)
  
  # compute the resample test statistic and quantile
  resample_test_statistic <- 1 / (sqrt(n)) * residual_star
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

