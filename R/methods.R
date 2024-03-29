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
  if (test_hyperparams$var_method_type == "oracle") {
    Var_X_given_Z <- data$Var_X_given_Z_oracle
    if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  } else {
    Var_X_given_Z <- fit_conditional_variance(X, Z, E_X_given_Z, test_hyperparams$var_method_type)
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
  # set the default hyperparameters
  test_hyperparams <- set_default_test_hyperparams("GCM", test_hyperparams)
  # extract the unlabel data
  if(test_hyperparams$way_to_learn == "semi_supervised"){
    if(is.null(data$X_unlabel)){
      stop("Unlabel data should be provided!")
    }else{
      X_unlabel <- as.matrix(data$X_unlabel)
      Z_unlabel <- as.matrix(data$Z_unlabel)
    }
  }else{
    X_unlabel <- NULL
    Z_unlabel <- NULL
  }
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)
  d <- ncol(Z)

  # union the unlabel data and label data
  combine_X <- rbind(as.matrix(X), X_unlabel)
  combine_Z <- rbind(as.matrix(Z), Z_unlabel)

  # fit conditional mean of X given Z
  if (X_on_Z_reg$mean_method_type == "oracle") {
    E_X_given_Z <- data$E_X_given_Z_oracle
    E_X_given_Z_label <- E_X_given_Z[1:length(X)]
    if (is.null(E_X_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    fit_X_given_Z <- fit_conditional_mean(combine_X, combine_Z, X_on_Z_reg)
    E_X_given_Z <- fit_X_given_Z$conditional_mean
    E_X_given_Z_label <- E_X_given_Z[1:length(X)]
  }

  # fit conditional mean of Y given Z
  if (Y_on_Z_reg$mean_method_type == "oracle") {
    E_Y_given_Z <- data$E_Y_given_Z_oracle
    E_Y_given_Z_label <- E_Y_given_Z[1:length(Y)]
    if (is.null(E_Y_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    fit_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)
    E_Y_given_Z <- fit_Y_given_Z$conditional_mean
  }

  # define the test statistic
  X_residuals <- X - E_X_given_Z_label
  Y_residuals <- Y - E_Y_given_Z
  R <- X_residuals * Y_residuals
  test_statistic <- sqrt(n) * mean(R) / sqrt(mean(R^2) - (mean(R))^2)

  # define the p-value (for now, define as two-sided)
  p_value <- 2 * stats::pnorm(abs(test_statistic), lower.tail = FALSE)

  # compute MSE only for Gaussian
  if(test_hyperparams$MSE){
    if(is.null(data$beta)){
      stop("Must specify the orcale beta coefficient when computing MSE")
    }
    if(is.null(data$gamma)){
      stop("Must specify the orcale gamma coefficient when computing MSE")
    }
    # extract the fitted coefficient
    gamma_hat <- as.vector(fit_X_given_Z$coef_vec)[-1]
    beta_hat <- as.vector(fit_Y_given_Z$coef_vec)[-1]
    # extract oracle conditional expectation
    oracle_X_given_Z <- data$E_X_given_Z_oracle
    oracle_Y_given_Z <- data$E_Y_given_Z_oracle
    # compute MSE on shared/total variables
    MSE_shared_X_on_Z <- MSE_shared(data = data, coef_hat = gamma_hat, type = "null", cond_distr = "X_on_Z")
    MSE_total_X_on_Z <- MSE_total(oracle_pred = oracle_X_given_Z, covariate = data$Z, coef_hat = gamma_hat)
    MSE_shared_Y_on_Z <- MSE_shared(data = data, coef_hat = beta_hat, type = data$setting, cond_distr = "Y_on_Z")
    MSE_total_Y_on_Z <- MSE_total(oracle_pred = oracle_Y_given_Z, covariate = data$Z, coef_hat = beta_hat)
  }

  # output the results
  if(!test_hyperparams$MSE){
    return(
      result <- data.frame(
        parameter = c("test_statistic", "p_value"),
        target = "conditional_independence",
        value = c(test_statistic, p_value))
    )
  }else{
    result <- data.frame(
      parameter = c("test_statistic", "p_value",
                    "MSE_shared_X_Z", "MSE_total_X_Z",
                    "MSE_shared_Y_Z", "MSE_total_Y_Z"),
      target = "conditional_independence",
      value = c(test_statistic, p_value,
                MSE_shared_X_on_Z, MSE_total_X_on_Z,
                MSE_shared_Y_on_Z, MSE_total_Y_on_Z))
  }
  result
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
  X_given_Z <- symcrt::fit_conditional_mean(X, Z, X_on_Z_reg)
  E_X_given_Z <- X_given_Z$conditional_mean # hat{E(X|Z)}
  
  # fit conditional mean of Y given Z
  Y_given_Z <- symcrt::fit_conditional_mean(Y, Z, Y_on_Z_reg)
  E_Y_given_Z <- Y_given_Z$conditional_mean # hat{E(Y|Z)}

  # compute the oracle residuals
  X_oracle_residuals <- X - cond_mean_X_Z # X - E(X|Z)
  Y_oracle_residuals <- Y - cond_mean_Y_Z # Y - E(Y|Z)
  cond_mean_diff_X_Z <- cond_mean_X_Z - E_X_given_Z # E(X|Z) - hat{E(X|Z)}
  cond_mean_diff_Y_Z <- cond_mean_Y_Z - E_Y_given_Z # E(Y|Z) - hat{E(Y|Z)}

  # compute three bias terms
  R <- (X-E_X_given_Z)*(Y-E_Y_given_Z)
  cross_bias <- sqrt(n)*mean(cond_mean_diff_X_Z*cond_mean_diff_Y_Z)/sqrt(mean(R^2) - (mean(R))^2)
  X_Z_bias <- sqrt(n)*mean(Y_oracle_residuals*cond_mean_diff_X_Z)/sqrt(mean(R^2) - (mean(R))^2)
  Y_Z_bias <- sqrt(n)*mean(X_oracle_residuals*cond_mean_diff_Y_Z)/sqrt(mean(R^2) - (mean(R))^2)

  # output the results
  out <- data.frame(parameter = c("cross_bias", "X_Z_bias", "Y_Z_bias"),
                    target = "bias",
                    value = c(cross_bias, X_Z_bias, Y_Z_bias))
  out$coef_X_on_Z <- list(as.vector(X_given_Z$coef_vec)[-1])
  out$coef_Y_on_Z <- list(as.vector(Y_given_Z$coef_vec)[-1])
  return(out)
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
  # set default test hyperparameters
  test_hyperparams <- set_default_test_hyperparams("dCRT", test_hyperparams)

  if(test_hyperparams$way_to_learn == "semi_supervised"){
    if(is.null(data$X_unlabel)){
      stop("Unlabel data should be provided!")
    }else{
      X_unlabel <- as.matrix(data$X_unlabel)
      Z_unlabel <- as.matrix(data$Z_unlabel)
    }
  }else{
    X_unlabel <- NULL
    Z_unlabel <- NULL
  }

  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)

  # union the unlabel data and label data
  combine_X <- rbind(as.matrix(X), X_unlabel)
  combine_Z <- rbind(as.matrix(Z), Z_unlabel)

  # fit conditional mean of X given Z
  if (X_on_Z_reg$mean_method_type == "oracle") {
    E_X_given_Z <- data$E_X_given_Z_oracle
    E_X_given_Z_label <- E_X_given_Z[1:length(X)]
    if (is.null(E_X_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    E_X_given_Z <- fit_conditional_mean(combine_X, combine_Z, X_on_Z_reg)$conditional_mean
    E_X_given_Z_label <- E_X_given_Z[1:length(X)]
  }

  # fit conditional variance of X given Z
  if (test_hyperparams$var_method_type == "oracle") {
    Var_X_given_Z <- data$Var_X_given_Z_oracle
    Var_X_given_Z_label <- Var_X_given_Z[1:length(X)]
    if (is.null(Var_X_given_Z)) stop("Must specify cond_var if var_method_type is oracle")
  } else {
    Var_X_given_Z <- fit_conditional_variance(combine_X, combine_Z, E_X_given_Z, test_hyperparams$var_method_type)
    Var_X_given_Z_label <- Var_X_given_Z[1:length(X)]
  }

  # fit conditional mean of Y given Z with labeled data
  E_Y_given_Z <- fit_conditional_mean(Y, Z, Y_on_Z_reg)$conditional_mean

  # compute residuals
  X_residuals <- c(X - E_X_given_Z_label)
  Y_residuals <- c(Y - E_Y_given_Z)

  # resample matrix from the specified distribution
  resample_matrix <- resample_dCRT(conditional_mean = E_X_given_Z_label,
                                   conditional_variance = Var_X_given_Z_label,
                                   no_resample = test_hyperparams$no_resample,
                                   resample_dist = test_hyperparams$resample_family)

  # compute the residuals and variance vector for each resample
  resample_X_residuals <- resample_matrix - E_X_given_Z_label
  residual_star <- apply(resample_X_residuals * Y_residuals, 2, sum)

  # define the test statistic with labeled data (with/without normalized GCM sd)
  if(!test_hyperparams$normalize){
    S_hat <- 1
    test_statistic <- 1 / (sqrt(n) * S_hat) * sum(X_residuals * Y_residuals)
  }else{
    R <- X_residuals * Y_residuals
    GCM_sd <- sqrt(mean(R^2) - (mean(R))^2)
    test_statistic <- 1 / (sqrt(n) * GCM_sd) * sum(X_residuals * Y_residuals)
    R_hat <- resample_X_residuals * Y_residuals
    S_hat <- sqrt(apply(R_hat^2, 2, mean) - apply(R_hat, 2, mean)^2)
  }

  # compute the resample test statistic and quantile
  resample_test_statistic <- 1 / (sqrt(n) * S_hat) * residual_star
  no_exceed <- length(which(abs(resample_test_statistic) >= abs(test_statistic)))

  # compute the p-value (two sided)
  p_value <- (no_exceed + 1) / (test_hyperparams$no_resample + 1)

  # output the results
  data.frame(
    parameter = c("test_statistic", "p_value"),
    target = "conditional_independence",
    value = c(test_statistic, p_value)
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
  test_hyperparams <- set_default_test_hyperparams("MaxwayCRT", test_hyperparams)
  # extract X, Y, Z from first input argument
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- nrow(Z)
  d <- ncol(Z)
  if(test_hyperparams$way_to_learn == "semi_supervised"){
    X_unlabel <- data$X_unlabel
    Z_unlabel <- data$Z_unlabel
    no_unlabel <- nrow(Z_unlabel)
    index_unlabel <- sample(1:n, no_unlabel)
    
    if(is.null(X_unlabel)){
      stop("Unlabel data should be provided!")
    }
    
    # extract the proportion of the number of holdout data and sample 1:n
    no_holdout <- round(n*test_hyperparams$holdout_prop)
    index_holdout <- sample(1:n, no_holdout)
    
    # extract the index that used for testing
    index_test <- setdiff(1:n, index_holdout)
  }else{
    if(test_hyperparams$unlabel_prop == 0){
      stop("Must specify the proportion of unlabel data in supervised setting!")
    }
    # extract the proportion of the number of unlabelled data and sample 1:n
    no_unlabel <- round(n*test_hyperparams$unlabel_prop)
    index_unlabel <- sample(1:n, no_unlabel)

    # assign the unlabel data
    X_unlabel <- X[index_unlabel]
    Z_unlabel <- Z[index_unlabel, ]

    # find the rest data index
    index_rest <- setdiff(1:n, index_unlabel)

    # extract the proportion of the number of holdout data and sample 1:n
    no_holdout <- round(n*test_hyperparams$holdout_prop)
    index_holdout <- sample(index_rest, no_holdout)

    # extract the index that used for testing
    index_test <- setdiff(index_rest, index_holdout)
  }








  # fit conditional mean of X given Z with unlablled data
  if (X_on_Z_reg$mean_method_type == "oracle") {
    E_X_given_Z_full <- data$E_X_given_Z_oracle
    E_X_given_Z <- E_X_given_Z_full[index_unlabel]
    if (is.null(E_X_given_Z)) stop("Must specify cond_mean if mean_method_type is oracle")
  } else {
    fit_X_given_Z <- fit_conditional_mean(X_unlabel,
                                          Z_unlabel,
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
  coef_Y_given_Z <- as.vector(fit_Y_given_Z$coef_vec)


  # extract the active set part in g(Z)
  if(test_hyperparams$g_Z_support == "support of beta_hat"){
    act_set_Y_given_Z <- which(coef_Y_given_Z[-1] != 0)
  }else if(test_hyperparams$g_Z_support == "default"){
    no_act_set_Y_given_Z <- as.integer(min(2 * log(d), d - 1, (no_unlabel)^(1/3)))
    act_set_Y_given_Z <- order(abs(coef_Y_given_Z[-1]), decreasing=TRUE)[1:no_act_set_Y_given_Z]
  }


  # concatenate predictor and active Z together
  predictor_Y_Z_unlabel <- cbind(1, Z_unlabel)%*%coef_Y_given_Z
  predictor_Y_Z_test <- cbind(1, Z[index_test, ])%*%coef_Y_given_Z
  if(length(act_set_Y_given_Z) == 0){
    g_Z_unlabel <- "zero"
    # g_Z_unlabel <- stats::rnorm(no_unlabel)
    # g_Z_test <- stats::rnorm(length(index_test))
  }else{
    g_Z_unlabel_act <- Z_unlabel[, act_set_Y_given_Z]
    g_Z_test_act <- Z[index_test, act_set_Y_given_Z]
    g_Z_unlabel <- orthogonalize(predictor_Y_Z_unlabel, Z_unlabel[, act_set_Y_given_Z])
    g_Z_test <- orthogonalize(predictor_Y_Z_test, Z[index_test, act_set_Y_given_Z])
  }
  # print(g_Z_unlabel[,ncol(g_Z_unlabel)])

  # fit the conditional expectation of Y|Z on label data
  Y_on_Z_reg_hyperparams <- Y_on_Z_reg$mean_method_hyperparams
  switch(Y_on_Z_reg_hyperparams$family,
         binomial = {
           E_Y_given_Z_test <- expit(cbind(1, Z[index_test, ])%*%coef_Y_given_Z)
         },
         gaussian = {
           E_Y_given_Z_test <- cbind(1, Z[index_test, ])%*%coef_Y_given_Z
         },
         {
           stop("The regression type is invalid!")
         }
  )
  

  # learning the resampling distribution with unlabelled data
  switch(test_hyperparams$resample_family,
         binomial = {
           # compute the labeled E(X|Z)
           E_X_given_Z_test <- as.vector(expit(cbind(1, Z[index_test, ])%*%coef_X_given_Z))

           # calibration
           if(any(g_Z_unlabel == "zero")){
             resample_mean <- E_X_given_Z_test
           }else{
             fit_X_on_g_h <- stats::glm(X_unlabel ~ g_Z_unlabel, offset = logit(E_X_given_Z), family = "binomial")
             if(any(is.na(fit_X_on_g_h$coefficients))){
               fit_X_on_g_h <- stats::glm(X_unlabel ~ g_Z_unlabel_act, offset = logit(E_X_given_Z), family = "binomial")
               # resample matrix from the specified distribution
               resample_mean <- expit(logit(E_X_given_Z_test) + as.vector(cbind(1, g_Z_test_act)%*%(fit_X_on_g_h$coefficients)))
             }else{
               # resample matrix from the specified distribution
               resample_mean <- expit(logit(E_X_given_Z_test) + as.vector(cbind(1, g_Z_test)%*%(fit_X_on_g_h$coefficients)))
             }
           }
           resample_var <- resample_mean*(1-resample_mean)
           resample_matrix <- resample_dCRT(conditional_mean = resample_mean,
                                            conditional_variance = resample_var,
                                            no_resample = test_hyperparams$no_resample,
                                            resample_dist = test_hyperparams$resample_family)
           resample_X_residuals <- resample_matrix - E_X_given_Z_test

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
         },
         gaussian = {

           # compute the residual and fit conditional mean on label data
           r_unlabel <- X_unlabel - cbind(1, Z_unlabel)%*%coef_X_given_Z
           r_label <- X[index_test] - cbind(1, Z[index_test, ])%*%coef_X_given_Z
           E_X_given_Z_test <- cbind(1, Z[index_test, ])%*%coef_X_given_Z
        

           if(any(g_Z_unlabel == "zero")){
             X_residuals <- r_label - rep(mean(r_label), length(r_label))
           }else{
             # fit the residual and g(Z) with lasso
             residual_on_gZ <- list(mean_method_type = "MLE",
                                    mean_method_hyperparams = list(family = "gaussian"),
                                    var_method_type = "homoskedastic")

             fit_residual_g_Z <- fit_conditional_mean(r_unlabel, g_Z_unlabel, residual_on_gZ)
             Var_residual_g_Z <- fit_conditional_variance(r_unlabel,
                                                          g_Z_unlabel,
                                                          fit_residual_g_Z$conditional_mean,
                                                          test_hyperparams$var_method_type)
             coef_residual_on_g <- fit_residual_g_Z$coef_vec
             if(any(is.na(coef_residual_on_g))){
               fit_residual_g_Z <- fit_conditional_mean(r_unlabel, g_Z_unlabel_act, residual_on_gZ)
               Var_residual_g_Z <- fit_conditional_variance(r_unlabel,
                                                            g_Z_unlabel_act,
                                                            fit_residual_g_Z$conditional_mean,
                                                            test_hyperparams$var_method_type)
               coef_residual_on_g <- fit_residual_g_Z$coef_vec
               # compute residuals
               X_residuals <- r_label - cbind(1, g_Z_test_act)%*%coef_residual_on_g
             }else{
               # compute residuals
               X_residuals <- r_label - cbind(1, g_Z_test)%*%coef_residual_on_g
             }
           }
           

           # compute p_value
           X_residuals <- X_residuals/(stats::sd(X_residuals))
           Y_residuals <- c(Y[index_test] - E_Y_given_Z_test)
           print(any(is.na(X_residuals)))
           n_test <- length(X_residuals)
           imp_obe <- mean(X_residuals * Y_residuals)
           emp_var <- mean(Y_residuals^2)
           test_statistic <- sqrt(n_test) * imp_obe / sqrt(emp_var)
           p_value <- 2 * stats::pnorm(- abs(test_statistic))

         },
         {
           stop("The rsampling distribution of X given Z is invalid!")
         }
  )

  # output the results
  data.frame(
    parameter = c("test_statistic", "p_value"),
    target = "conditional_independence",
    value = c(test_statistic, p_value)
  )
}
