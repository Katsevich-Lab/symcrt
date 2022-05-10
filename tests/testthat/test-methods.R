test_that("dCRT works", {

  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  n <- 100
  d <- 5
  gamma <- 2 * rbinom(d, 1, 0.5) - 1
  Z <- matrix(rnorm(n = n * d), nrow = n, ncol = d)
  noise <- rnorm(2 * n)
  E_X_given_Z_oracle <- Z %*% gamma |> as.vector()
  Var_X_given_Z_oracle <- rep(1, n)
  X <- E_X_given_Z_oracle + noise[1:n]
  Y <- Z %*% gamma + noise[(n + 1):(2 * n)]
  data <- list(
    Z = Z,
    X = X,
    Y = Y,
    E_X_given_Z_oracle = E_X_given_Z_oracle,
    Var_X_given_Z_oracle = Var_X_given_Z_oracle
  )
  X_on_Z_reg <- list(
    mean_method_type = "oracle",
    mean_method_hyperparams = list(
      family = "gaussian"
    ),
    var_method_type = "oracle"
  )
  
  Y_on_Z_reg <- list(
    mean_method_type = "MLE",
    mean_method_hyperparams = list(
      family = "gaussian"
    ),
    var_method_type = "squared_residuals"
  )
  
  # test whether we get similar p-value from MX2 and dCRT
  # run MX2
  MX2_result <- MX2_F_test(
    data = data,
    X_on_Z_reg = X_on_Z_reg,
    Y_on_Z_reg = Y_on_Z_reg
  )
  p_value_MX2 <- MX2_result |>
    dplyr::filter(parameter == "p_value") |>
    dplyr::select(value) |>
    unname() |>
    as.numeric()
  # run dCRT
  dCRT_result <- dCRT(
    data = data,
    X_on_Z_reg = X_on_Z_reg,
    Y_on_Z_reg = Y_on_Z_reg,
    test_hyperparams = list()
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
  set.seed(1)
  # test whether a specific example gives exactly the right answer
  # (generated using dput)
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
  
  X_on_Z_reg <- list(
    mean_method_type = "LASSO",
    mean_method_hyperparams = list(
      family = "gaussian",
      s = "lambda.min",
      nfols = 10,
      alpha = 1
    )
  )
  
  Y_on_Z_reg <- list(
    mean_method_type = "LASSO",
    mean_method_hyperparams = list(
      family = "gaussian",
      s = "lambda.min",
      nfols = 10,
      alpha = 1
    )
  )
  
  # obtain the oracle fitted mean
  lasso_fit_X_Z <- glmnet::cv.glmnet(x = Z, y = X)
  lasso_fit_Y_Z <- glmnet::cv.glmnet(x = Z, y = Y)
  fitted_X_Z <- stats::predict(lasso_fit_X_Z, newx = Z, s = "lambda.min")
  fitted_Y_Z <- stats::predict(lasso_fit_Y_Z, newx = Z, s = "lambda.min")
  
  # compute the three bias term
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
            0.2)
})




test_that("MaxwayCRT works for Gaussian", {
  
  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  source("function.R")
  set.seed(1)
  n <- 100
  d <- 1000
  s <- 5
  gamma <- numeric(d)
  gamma[1:s] <- 2 * rbinom(s, 1, 0.5) - 1
  Z <- matrix(rnorm(n = n * d), nrow = n, ncol = d)
  noise <- rnorm(2 * n)
  E_X_given_Z_oracle <- Z %*% gamma |> as.vector()
  Var_X_given_Z_oracle <- rep(1, n)
  X <- E_X_given_Z_oracle + noise[1:n]
  Y <- Z %*% gamma + noise[(n + 1):(2 * n)]
  data <- list(
    Z = Z,
    X = X,
    Y = Y,
    E_X_given_Z_oracle = E_X_given_Z_oracle,
    Var_X_given_Z_oracle = Var_X_given_Z_oracle
  )
  X_on_Z_reg <- list(
    mean_method_type = "LASSO",
    mean_method_hyperparams = list(
      s = "lambda.min",
      family = "gaussian"
    ),
    var_method_type = "homoskedastic"
  )
  
  Y_on_Z_reg <- list(
    mean_method_type = "LASSO",
    mean_method_hyperparams = list(
      s = "lambda.min",
      family = "gaussian"
    ),
    var_method_type = "homoskedastic"
  )
  
  test_hyperparams <- list(resample_family = "gaussian", 
                           no_resample = 5000,
                           unlabel_prop = 0.5)
  
  
  # test whether we get similar p-value from MaxwayCRT and code by Molei
  # output the results
  MaxwayCRT_result <- MaxwayCRT(
    data = data,
    X_on_Z_reg = X_on_Z_reg,
    Y_on_Z_reg = Y_on_Z_reg,
    test_hyperparams = test_hyperparams
  )
  p_value_MaxwayCRT <- MaxwayCRT_result |>
    dplyr::select(value) |>
    unname() |>
    as.numeric()
  
  # run Molei's code
  no_unlabel <- round(n*0.5)
  set.seed(1)
  index_unlabel <- sample(1:n, no_unlabel)
  index_test <- setdiff(1:n, index_unlabel)
  A_add <- X[index_unlabel]
  Z_add <- Z[index_unlabel, ]
  A_label <- X[index_test]
  Z_label <- Z[index_test, ]
  y_label <- Y[index_test]
  MaxwayCRT <- Maxway_CRT(A_label, Z_label, y_label, A_add, Z_add, 
                          Z_label_add = NULL, y_label_add = NULL,
                          model_x = 'Gaussian_lasso', model_y = 'Gaussian_lasso',
                          lambda.seq = NULL, k = NULL, seed = 1, M = 5000)
  
  Molei_MaxwayCRT <- MaxwayCRT$Maxway_d0CRT_pvl
  expect_lt(
    abs(Molei_MaxwayCRT - p_value_MaxwayCRT),
    0.05
  )
})



test_that("MaxwayCRT works for Gaussian", {
  
  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  source("function.R")
  set.seed(1)
  n <- 100
  d <- 1000
  s <- 5
  rho <- 0.5
  N <- n
  gamma <- numeric(d)
  gamma[1:s] <- 2 * rbinom(s, 1, 0.5) - 1
  sig_Z <- katlabutils::generate_cov_ar1(rho, d)
  Z <- katlabutils::fast_generate_mvn(numeric(d), sig_Z, N)
  noise <- rnorm(2 * n)
  E_X_given_Z_oracle <- Z %*% gamma |> as.vector()
  Var_X_given_Z_oracle <- rep(1, n)
  X <- E_X_given_Z_oracle + noise[1:n]
  Y <- Z %*% gamma + noise[(n + 1):(2 * n)]
  Z_unlabel <- matrix(rnorm(n = N * d), nrow = N, ncol = d)
  E_X_given_Z_unlabel <- Z_unlabel %*% gamma |> as.vector()
  X_unlabel <- E_X_given_Z_unlabel + rnorm(N)
  data <- list(
    Z = Z,
    X = X,
    Y = Y,
    E_X_given_Z_oracle = E_X_given_Z_oracle,
    Var_X_given_Z_oracle = Var_X_given_Z_oracle,
    Z_unlabel = Z_unlabel,
    X_unlabel = X_unlabel
  )
  X_on_Z_reg <- list(
    mean_method_type = "LASSO",
    mean_method_hyperparams = list(
      s = "lambda.min",
      family = "gaussian"
    ),
    var_method_type = "homoskedastic"
  )
  
  Y_on_Z_reg <- list(
    mean_method_type = "LASSO",
    mean_method_hyperparams = list(
      s = "lambda.min",
      family = "gaussian"
    ),
    var_method_type = "homoskedastic"
  )
  
  test_hyperparams <- list(resample_family = "gaussian", 
                           no_resample = 5000,
                           unlabel_prop = 0.5)
  
  
  # test whether we get similar p-value from MaxwayCRT and code by Molei
  # output the results
  MaxwayCRT_result <- MaxwayCRT(
    data = data,
    X_on_Z_reg = X_on_Z_reg,
    Y_on_Z_reg = Y_on_Z_reg,
    test_hyperparams = test_hyperparams
  )
  p_value_MaxwayCRT <- MaxwayCRT_result |>
    dplyr::select(value) |>
    unname() |>
    as.numeric()
  
  # run Molei's code
  no_unlabel <- round(n*0.5)
  set.seed(1)
  index_unlabel <- sample(1:n, no_unlabel)
  index_test <- setdiff(1:n, index_unlabel)
  A_add <- X[index_unlabel]
  Z_add <- Z[index_unlabel, ]
  A_label <- X[index_test]
  Z_label <- Z[index_test, ]
  y_label <- Y[index_test]
  MaxwayCRT <- Maxway_CRT(A_label, Z_label, y_label, A_add, Z_add, 
                          Z_label_add = NULL, y_label_add = NULL,
                          model_x = 'Gaussian_lasso', model_y = 'Gaussian_lasso',
                          lambda.seq = NULL, k = NULL, seed = 1, M = 5000)
  
  Molei_MaxwayCRT <- MaxwayCRT$Maxway_d0CRT_pvl
  expect_lt(
    abs(Molei_MaxwayCRT - p_value_MaxwayCRT),
    0.05
  )
})


