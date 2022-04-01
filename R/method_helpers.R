#' Fit conditional mean.
#'
#' Function that fits some regression method to get a conditional mean of a
#' response (could be Y or X) given features (usually Z).
#'
#' @param response A response vector of length n
#' @param features A design matrix of dimension nxp
#' @param method   A string that specifies which regression method to use, e.g. \code{OLS} or \code{LASSO}
#' @param family   A string that specifies the link function e.g.\code{Gaussian} or \code{Binomial}
#'
#' @return A vector of fitted conditional means
#' @export
fit_conditional_mean <- function(response, features, method) {
  method_type <- method$method_type
  hyperparams <- method$hyperparams
  family <- method$family
  switch(method_type,
    MLE = {
      if (is.null(family)) {
        stop("A family in GLM should be specified!")
      }
      GLM_fit <- stats::glm(response ~ features, family = family)
      GLM_fit$fitted.values |> unname()
    },
    LASSO = {
      lasso_fit <- fit_lasso(response, features, hyperparams)
      # generate predictions
      stats::predict(lasso_fit,
        newx = features,
        type = "response",
        s = s,
      ) |> as.vector()
    },
    PLASSO = {
      # parse s hyperparameter
      if (is.null(hyperparams$s)) {
        s <- "lambda.1se"
      } else {
        s <- hyperparams$s
      }
      lasso_fit <- fit_lasso(response, features, hyperparams)
      coefs <- stats::coef(lasso_fit, s = s)
      act_set <- which(coefs[-1, 1] != 0) |> unname()
      if (length(act_set) == 0) {
        glm_fit <- stats::glm(response ~ 1, family = family)
      } else {
        glm_fit <- stats::glm(response ~ features[, act_set], family = family)
      }
      glm_fit$fitted.values |> unname()
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
  # set default values of hyperparameters
  if (is.null(hyperparams$s)) {
    s <- "lambda.1se"
  } else {
    s <- hyperparams$s
  }
  if (is.null(hyperparams$nfolds)) {
    nfolds <- 10
  } else {
    nfolds <- hyperparams$nfolds
  }
  if (is.null(hyperparams$alpha)) {
    alpha <- 1
  } else {
    alpha <- hyperparams$alpha
  }
  if (is.null(family)) {
    family <- "gaussian"
  } else{
    family <- hyperparams$family
  }
  # run lasso
  glmnet::cv.glmnet(x = features, y = response, nfolds = nfolds, alpha = alpha, family = family)
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
    Binom = {
      n <- length(conditional_mean)
      matrix(stats::rbinom(n * no_resample, 1, prob = rep(conditional_mean, no_resample)),
        nrow = n,
        ncol = no_resample
      )
    },
    Gaussian = {
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
