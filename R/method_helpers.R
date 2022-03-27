#' Fit conditional mean.
#'
#' Function that fits some regression method to get a conditional mean of a
#' response (could be Y or X) given features (usually Z).
#'
#' @param response A response vector of length n
#' @param features A design matrix of dimension nxp
#' @param method   A string that specifies which regression method to use, e.g. \code{OLS} or \code{LASSO} or \code{PLASSO}
#'
#' @return A vector of fitted conditional means
#' @export
fit_conditional_mean <- function(response, features, method) {
  switch(method,
    OLS = {
      lm_fit <- stats::lm(response ~ features)
      lm_fit$fitted.values |> unname()
    },
    LASSO = {
      lasso_fit <- glmnet::cv.glmnet(x = features, y = response)
      stats::predict(lasso_fit, newx = features, s = "lambda.1se") |> as.vector()
    },
    PLASSO = {
      lasso_fit <- glmnet::cv.glmnet(x = features, y = response)
      act_set <- which(coef(lasso_fit, s = 'lambda.1se')[-1,]!= 0) |> unname()
      lm_fit <- stats::lm(response ~ features[, act_set])
      lm_fit$fitted.values |> unname()
    },
    zero = {
      # wrong estimate!
      n <- nrow(features)
      rep(0, n)
    },
    {
      stop("Invalid specification of regression method.")
    }
  )
}

# helper function to fit a conditional variance of a response (could be Y or X) on
# a set of predictors (usually Z)
#' Title
#'
#' @param response A response vector of length n
#' @param features A design matrix of dimension nxp
#' @param conditional_mean The conditional mean vector of length n 
#' @param method A string that specifies the variance estimation method to use, e.g. \code{squared_residual}
#'
#' @return A vector of estimated conditional variance of length n
#' @export

fit_conditional_variance <- function(response, features, conditional_mean, method) {
  switch(method,
    squared_residual = {
      (response - conditional_mean)^2
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
resample_dCRT <- function(conditional_mean, conditional_variance = NULL, no_resample = 1000, resample_dist){
  switch(resample_dist,
         Binom = {
           n <- length(conditional_mean)
           matrix(rbinom(n*no_resample, 1, prob = rep(conditional_mean, no_resample)),
                               nrow = n, 
                               ncol = no_resample)
         },
         Heter_Gaussian = {
           t(fast_data_generate(mean = conditional_mean, covariance = diag(conditional_variance), 
                                B = no_resample, 
                                n = 1))
         },
         Homo_Gaussian = {
           variance_estimate <- mean(conditional_variance)
           variance <- rep(variance_estimate, length(conditional_variance))
           t(fast_data_generate(mean = conditional_mean, covariance = diag(variance), 
                                B = no_resample, 
                                n = 1))
         },
         {
           stop("Invalid specification of resampling method.")
         }
  )
}
