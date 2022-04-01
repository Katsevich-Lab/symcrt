test_that("dCRT_f works", {

  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  n <- 100
  d <- 5
  gamma <- 2 * rbinom(d, 1, 0.5) - 1
  Z <- matrix(rnorm(n = n * d), nrow = n, ncol = d)
  noise <- rnorm(2 * n)
  cond_mean <- Z %*% gamma |> as.vector()
  cond_var <- rep(1, n)
  X <- cond_mean + noise[1:n]
  Y <- Z %*% gamma + noise[(n + 1):(2 * n)]
  data <- list(
    Z = Z,
    X = X,
    Y = Y,
    cond_mean = cond_mean,
    cond_var = cond_var
  )
  X_on_Z_reg <- "oracle"
  Y_on_Z_reg <- "OLS"
  X_on_Z_var <- "oracle"
  resample_dist <- "Gaussian"
  no_resample <- 10000
  # test whether we get similar p-value from MX2 and dCRT
  # run MX2
  MX2_result <- MX2_F_test_f(
    data = data,
    X_on_Z_reg = X_on_Z_reg,
    Y_on_Z_reg = Y_on_Z_reg,
    X_on_Z_var = X_on_Z_var
  )
  p_value_MX2 <- MX2_result |>
    dplyr::filter(parameter == "p_value") |>
    dplyr::select(value) |>
    unname() |>
    as.numeric()
  # run dCRT
  dCRT_result <- dCRT_f(
    data = data,
    X_on_Z_reg = X_on_Z_reg,
    Y_on_Z_reg = Y_on_Z_reg,
    X_on_Z_var = X_on_Z_var,
    resample_dist = resample_dist,
    no_resample = no_resample
  )
  p_value_dCRT <- dCRT_result |>
    dplyr::select(value) |>
    unname() |>
    as.numeric()
  expect_lt(
    abs(p_value_MX2 - p_value_dCRT),
    0.01
  )
})

test_that("GCM_debug works", {
  
  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  set.seed(1)
  n <- 2000
  d <- 200
  no_act <- 5
  gamma <- rep(0, d)
  beta <- rep(0, d)
  gamma[1:no_act] <- 2*rbinom(no_act, 1, 0.5) - 1
  beta[1:no_act] <- 2*rbinom(no_act, 1, 0.5) - 1
  Z <- matrix(rnorm(n = n*d), nrow = n, ncol = d)
  noise <- rnorm(2*n)
  cond_mean_X_Z <- Z %*% gamma |> as.vector()
  cond_mean_Y_Z <- Z %*% beta |> as.vector()
  X <- cond_mean_X_Z + noise[1:n]
  Y <- cond_mean_Y_Z + noise[(n+1):(2*n)]
  data <- list(Z = Z, 
               X = X, 
               Y = Y,
               cond_mean_X_Z = cond_mean_X_Z,
               cond_mean_Y_Z = cond_mean_Y_Z)
  X_on_Z_reg <- "LASSO"
  Y_on_Z_reg <- "LASSO"
  lasso_fit_X_Z <- glmnet::cv.glmnet(x = Z, y = X)
  lasso_fit_Y_Z <- glmnet::cv.glmnet(x = Z, y = Y)
  fitted_X_Z <- stats::predict(lasso_fit_X_Z, newx = Z, s = "lambda.1se")
  fitted_Y_Z <- stats::predict(lasso_fit_Y_Z, newx = Z, s = "lambda.1se")
  
  cross_bias <- sqrt(n)*mean((cond_mean_X_Z - fitted_X_Z)*(cond_mean_Y_Z - fitted_Y_Z))
  X_Z_bias <- sqrt(n)*mean((Y - cond_mean_Y_Z)*(cond_mean_X_Z - fitted_X_Z))
  Y_Z_bias <- sqrt(n)*mean((X - cond_mean_X_Z)*(cond_mean_Y_Z - fitted_Y_Z)) 
  oracle_bias <- matrix(c(cross_bias, X_Z_bias, Y_Z_bias), nrow = 3, ncol = 1)
  # test whether we get similar bias
  # run GCM_debug
  GCM_debug_result <- GCM_debug(data = data,
                        X_on_Z_reg = X_on_Z_reg,
                        Y_on_Z_reg = Y_on_Z_reg)
  GCM_bias <- GCM_debug_result |> 
    dplyr::select(value) |> 
    unname() |> as.matrix()
  expect_lt(max(abs(GCM_bias - oracle_bias)),
            0.1)
})

