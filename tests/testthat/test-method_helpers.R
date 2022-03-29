
test_that("fit_conditional_variance works",{
  # check whether the fit_conditional_variance generates right sequence with 
  # homoskedastic argument in OLS setting
  set.seed(1)
  n <- 1000
  features <- rnorm(n)
  response <- features*1 + rnorm(n)
  lm_fit <- stats::lm(response ~ features)
  conditional_mean <- lm_fit$fitted.values |> unname()
  method <- "homoskedastic"
  true_cond_var <- rep(1, n)
  fit_cond_var <- fit_conditional_variance(response, 
                                           features, 
                                           conditional_mean, 
                                           method)
  expect_lt(max(abs(true_cond_var - fit_cond_var)),
            0.1)
})



test_that("resample_dCRT works",{
  # check whether the resample_dCRT generates right samples with 
  # Gaussian argument within OLS setting
  features <- rnorm(100)
  response <- features*1 + rnorm(100)
  lm_fit <- stats::lm(response ~ features)
  conditional_mean <- lm_fit$fitted.values |> unname()
  method <- "homoskedastic"
  conditional_variance <- fit_conditional_variance(response, 
                                                   features, 
                                                   conditional_mean, 
                                                   method)
  no_resample <- 5000
  resample_dist <- "Gaussian"
  resample_matrix <- resample_dCRT(conditional_mean = conditional_mean,
                                   conditional_variance = conditional_variance,
                                   no_resample = no_resample,
                                   resample_dist = resample_dist)
  resample_mean <- apply(resample_matrix, 1, mean)
  expect_lt(max(abs(resample_mean - conditional_mean)),
            0.1)
  resample_var <- apply(resample_matrix, 1, var)
  expect_lt(max(abs(resample_var - conditional_variance)),
            0.1)
})



