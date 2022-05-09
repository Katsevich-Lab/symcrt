#' Fit conditional mean.
#'
#' Function that fits some regression method to get a conditional mean of a
#' response (could be Y or X) given features (usually Z).
#'
#' @param response A response vector of length n
#' @param features A design matrix of dimension nxp
#' @param method   A string that specifies which regression method to use, e.g. \code{OLS} or \code{LASSO}
#' @return A vector of fitted conditional means
#' @export
fit_conditional_mean <- function(response, features, method) {
  method_type <- method$mean_method_type
  hyperparams <- method$mean_method_hyperparams
  hyperparams <- set_default_method_hyperparams(method_type, hyperparams)
  switch(method_type,
    MLE = {
      GLM_fit <- stats::glm(response ~ features, family = hyperparams$family)
      list(conditional_mean = GLM_fit$fitted.values,
           coef_vec = as.vector(GLM_fit$coefficients))
    },
    LASSO = {
      lasso_fit <- fit_lasso(response, features, hyperparams)
      # generate predictions
      list(conditional_mean = as.vector(stats::predict(lasso_fit,
                                             newx = features,
                                             type = "response",
                                             s = hyperparams$s)),
           coef_vec = as.vector(stats::coef(lasso_fit, s = hyperparams$s))
           )
      
    },
    PLASSO = {
      lasso_fit <- fit_lasso(response, features, hyperparams)
      coefs <- stats::coef(lasso_fit, s = hyperparams$s)
      act_set <- which(coefs[-1, 1] != 0) |> unname()
      if (length(act_set) == 0) {
        glm_fit <- stats::glm(response ~ 1, family = hyperparams$family)
      } else {
        glm_fit <- stats::glm(response ~ features[, act_set], family = hyperparams$family)
      }
      list(conditional_mean = glm_fit$fitted.values,
           coef_vec = as.vector(glm_fit$coefficients))
    },
    zero = {
      # wrong estimate!
      n <- nrow(features)
      numeric(n)
    },
    {
      stop("Invalid specification of regression method.")
    }
  )
}

# TODO: document this function
fit_lasso <- function(response, features, hyperparams){
  # run lasso
  if(hyperparams$nfolds*3 > length(response)){
    glmnet::cv.glmnet(x = features,
                      y = response,
                      alpha = hyperparams$alpha,
                      family = hyperparams$family)
  }else{
    glmnet::cv.glmnet(x = features,
                      y = response,
                      nfolds = hyperparams$nfolds,
                      alpha = hyperparams$alpha,
                      family = hyperparams$family) 
  }
}

# helper function to fit a conditional variance of a response (could be Y or X) on
# a set of predictors (usually Z)
#' Title
#'
#' @param response A response vector of length n
#' @param features A design matrix of dimension nxp
#' @param conditional_mean The conditional mean vector of length n
#' @param method_type A string that specifies the variance estimation method to use, e.g. \code{squared_residual}
#'
#' @return A vector of estimated conditional variances of length n
#' @export

fit_conditional_variance <- function(response, features, conditional_mean, method_type) {
  switch(method_type,
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

#' resampling function to get resample with fitted mean and variance
#'
#' @param conditional_mean A vector of result from \code{fit_conditional_mean}
#' @param conditional_variance A vector of result from \code{fit_conditional_variance}
#' @param no_resample A number that specifies the size of resamples
#' @param resample_dist A distribution that specifies the resampling distribution \code{Binom} or \code{Gaussian}
#'
#' @return A matrix of nxno_resample including the resamples
#' @export
resample_dCRT <- function(conditional_mean, conditional_variance = NULL, no_resample = 1000, resample_dist) {
  # TODO: Write unit tests for this function
  switch(resample_dist,
    binomial = {
      n <- length(conditional_mean)
      matrix(stats::rbinom(n * no_resample, 1, prob = rep(conditional_mean, no_resample)),
        nrow = n,
        ncol = no_resample
      )
    },
    gaussian = {
      n <- length(conditional_mean)
      matrix(stats::rnorm(n * no_resample,
        mean = rep(conditional_mean, no_resample),
        sd = sqrt(rep(conditional_variance, no_resample))
      ),
      nrow = n,
      ncol = no_resample
      )
    },
    {
      stop("Invalid specification of resampling method.")
    }
  )
}

#' Populate default hyperparameters for each method
#'
#' @param method_type The method type (MLE, LASSO, etc.)
#' @param hyperparams The named list of hyperparams
#'
#' @return The same named list that was input, with default values added if necessary
#' @export
set_default_method_hyperparams <- function(method_type, hyperparams){
  switch(method_type,
         MLE={
           # currently no hyperparams
         },
         LASSO={
           if (is.null(hyperparams$s)) {
             hyperparams$s <- "lambda.1se"
           }
           if (is.null(hyperparams$nfolds)) {
             hyperparams$nfolds <- 10
           }
           if (is.null(hyperparams$alpha)) {
             hyperparams$alpha <- 1
           }
           if (is.null(hyperparams$family)) {
             hyperparams$family <- "gaussian"
           }
         },
         PLASSO={
           if (is.null(hyperparams$s)) {
             hyperparams$s <- "lambda.1se"
           }
           if (is.null(hyperparams$nfolds)) {
             hyperparams$nfolds <- 10
           }
           if (is.null(hyperparams$alpha)) {
             hyperparams$alpha <- 1
           }
           if (is.null(hyperparams$family)) {
             hyperparams$family <- "gaussian"
           }
         },
         zero={
           # currently no hyperparams
         }
         )
  hyperparams
}

#' Populate default hyperparameters for each test
#'
#' @param method_type The method type (GCM, dCRT, etc.)
#' @param hyperparams The named list of hyperparams
#'
#' @return The same named list that was input, with default values added if necessary
#' @export
set_default_test_hyperparams <- function(method_type, hyperparams){
  switch(method_type,
         GCM={
           # currently no hyperparams
         },
         GCM_debug={
           # currently no hyperparams
         },
         dCRT={
           if(is.null(hyperparams$resample_family)){
             hyperparams$resample_family <- "gaussian"
           }
           if(is.null(hyperparams$no_resample)){
             hyperparams$no_resample <- 2000
           }
           if(is.null(hyperparams$unlabel_prop)){
             hyperparams$unlabel_prop <- 0
           }
           if(is.null(hyperparams$holdout_prop)){
             hyperparams$holdout_prop <- 0
           }
           if(is.null(hyperparams$g_Z_support)){
             hyperparams$g_Z_support <- "default"
           }
           if(hyperparams$holdout_prop + hyperparams$unlabel_prop >= 1){
             stop("The number of unlablled and holdout data exceeds the total number of data!")
           }
         },
         MX2_F_test = {
           # currently no hyperparams
         }
  )
  hyperparams
}


#' Title
#'
#' @param gZ The predictor vector
#' @param Z_sub The submatrix selected by user for g(Z)
#'
#' @return A converted, stacked g(Z)
#' @export

orthogonalize <- function(gZ, Z_sub){
  model.fit <- stats::lm(gZ ~ Z_sub)
  Z_convert <- cbind(Z_sub, model.fit$residual)
  return(Z_convert)
}

#' Logit function
#'
#' @param x Input
#'
#' @return Logit value
#' @export

logit <- function(x){
  return(log(x / (1 - x)))
}

#' Inverse logit function (expit)
#'
#' @param x Input
#'
#' @return Expit output
#' @export

glogit <- function(x){
  return(1 / (1 + exp(-x)))
}

