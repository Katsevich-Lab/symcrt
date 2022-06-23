
test_that("fit_conditional_variance works", {
  # check whether the fit_conditional_variance generates right sequence with
  # homoskedastic argument in OLS setting
  set.seed(1)
  
  # generate data
  n <- 1000
  features <- rnorm(n)
  response <- features * 1 + rnorm(n)
  lm_fit <- stats::lm(response ~ features)
  conditional_mean <- lm_fit$fitted.values |> unname()
  method <- "homoskedastic"
  true_cond_var <- rep(1, n)
  fit_cond_var <- fit_conditional_variance(
    response,
    features,
    conditional_mean,
    method
  )
  expect_lt(
    max(abs(true_cond_var - fit_cond_var)),
    0.1
  )
})



test_that("resample_dCRT works", {
  # check whether the resample_dCRT generates right samples with
  # Gaussian argument within OLS setting
  features <- rnorm(100)
  response <- features * 1 + rnorm(100)
  lm_fit <- stats::lm(response ~ features)
  conditional_mean <- lm_fit$fitted.values |> unname()
  method <- "homoskedastic"
  conditional_variance <- fit_conditional_variance(
    response,
    features,
    conditional_mean,
    method
  )
  no_resample <- 5000
  resample_dist <- "gaussian"
  resample_matrix <- resample_dCRT(
    conditional_mean = conditional_mean,
    conditional_variance = conditional_variance,
    no_resample = no_resample,
    resample_dist = resample_dist
  )
  resample_mean <- apply(resample_matrix, 1, mean)
  expect_lt(
    max(abs(resample_mean - conditional_mean)),
    0.1
  )
  resample_var <- apply(resample_matrix, 1, var)
  expect_lt(
    max(abs(resample_var - conditional_variance)),
    0.1
  )
})



test_that("fit_conditional_mean works for PLASSO", {
  # check whether the fit_conditional_mean works with PLASSO
  set.seed(100)
  # test high-dimensional logistic regression with post fitting
  # specify method list with mean_method_type, mean_method_hyperparams
  # set the all the hyperparams but family to be default
  family <- "binomial"
  s <- "lambda.min"
  nfolds <- 10
  alpha <- 1
  method <- list(mean_method_type = "PLASSO",
                 mean_method_hyperparams = list(
                   family = family,
                   s = s,
                   nfolds = nfolds,
                   alpha = alpha
                 ))
  
  # generate the data
  n <- 400
  d <- 200
  no_act <- 5
  if (no_act == 0) {
    beta <- rep(0, d)
  } else {
    beta <- rep(0, d)
    beta[1:no_act] <- 2 * rbinom(no_act, 1, 0.5) - 1
  }
  features <- matrix(rnorm(n * d), nrow = n, ncol = d)
  cond_mean <- exp(features %*% beta) / (1 + exp(features %*% beta))
  response <- rbinom(n, 1, prob = cond_mean)
  
  # run the oracle post PMLE
  glm_fit <- glmnet::cv.glmnet(x = features, y = response, family = family)
  coef_glm <- as.vector(stats::coef(glm_fit, s = "lambda.min"))
  act_set <- which(coef_glm[-1] != 0) |> unname()
  if (length(act_set) == 0) {
    feature_intercept <- rep(coef_glm[1, 1], n)
    p_glm <- stats::glm(response ~ feature_intercept, family = family)
  } else {
    p_glm <- stats::glm(response ~ features[, act_set], family = family)
  }
  true_fitted <- p_glm$fitted.values |> unname()
  coef_glm[act_set] <- as.vector(p_glm$coefficients)[-1]
  coef_glm[1] <- as.vector(p_glm$coefficients)[1]
  
  # try symcrt function
  result <- fit_conditional_mean(
    response = response,
    features = features,
    method = method
  )
  conditional_mean <- result$conditional_mean
  coeffi <- as.vector(result$coef_vec)
  expect_lt(
    mean(abs(true_fitted - conditional_mean)^2),
    0.05
  )
  expect_lt(
    sum(abs(coeffi - coef_glm)^2),
    0.05
  )
})


test_that("fit_conditional_mean works for MLE", {
  # check whether the fit_conditional_mean works with GLM method
  set.seed(100)
  # test logistic regression
  # specify method list with mean_method_type, mean_method_hyperparams
  family <- "binomial"
  method <- list(mean_method_type = "MLE",
                 mean_method_hyperparams = list(
                   family = family
                 ))
  # generate the data
  n <- 400
  d <- 20
  beta <- 2 * rbinom(d, 1, 0.5) - 1
  features <- matrix(rnorm(n * d), nrow = n, ncol = d)
  cond_mean <- exp(features %*% beta) / (1 + exp(features %*% beta))
  response <- rbinom(n, 1, prob = cond_mean)
  
  # run the oracle MLE and get fitted values
  glm_fit <- stats::glm(response ~ features, family = family)
  true_fitted <- glm_fit$fitted.values |> unname()
  
  # run the symcrt function
  conditional_mean <- fit_conditional_mean(
    response = response,
    features = features,
    method = method
  )$conditional_mean
  expect_lt(
    mean(abs(true_fitted - conditional_mean)^2),
    0.05
  )
})


test_that("fit_conditional_mean works for PMLE", {
  # check whether the fit_conditional_mean works with penalized GLM method
  set.seed(100)
  # test high-dimensional logistic regression
  # specify method list with mean_method_type, mean_method_hyperparams
  # set the all the hyperparams but family to be default
  family <- "binomial"
  method <- list(mean_method_type = "LASSO",
                 mean_method_hyperparams = list(
                   family = family
                 ))
  
  # generate the data
  n <- 400
  d <- 200
  no_act <- 5
  if (no_act == 0) {
    beta <- rep(0, d)
  } else {
    beta <- rep(0, d)
    beta[1:no_act] <- 2 * rbinom(no_act, 1, 0.5) - 1
  }
  features <- matrix(rnorm(n * d), nrow = n, ncol = d)
  cond_mean <- exp(features %*% beta) / (1 + exp(features %*% beta))
  response <- rbinom(n, 1, prob = cond_mean)
  
  # obtain the oracle fitted values
  glm_fit <- glmnet::cv.glmnet(x = features, y = response, family = family)
  true_fitted <- stats::predict(glm_fit,
    newx = features,
    type = "response",
    s = "lambda.1se"
  ) |> unname()
  
  # try symcrt function
  conditional_mean <- fit_conditional_mean(
    response = response,
    features = features,
    method = method
  )$conditional_mean
  expect_lt(
    mean(abs(true_fitted - conditional_mean)^2),
    0.05
  )
})
