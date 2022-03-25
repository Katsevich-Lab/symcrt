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
