# helper function to fit a conditional mean of a response (could be Y or X) on
# a set of predictors (usually Z)
fit_conditional_mean <- function(response, features, method) {
  switch(method,
    OLS = {
      lm_fit <- stats::lm(response ~ features)
      lm_fit$fitted.values |> unname()
    },
    LASSO = {
      lasso_fit <- glmnet::cv.glmnet(x = features, y = response)
      glmnet::predict.glmnet(lasso_fit, newx = features, s = "lambda.1se") |> as.vector()
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
