# test_that("generate_method_list works", {
#   # check whether function works with different methods input
#   n <- 200
#   d <- 150
#   no_act <- 20
#   B <- 1
#   sig_Z <- diag(d)
# 
#   # fake p_grid
#   p_grid <- expand.grid(b = c(1))
# 
#   # fixed parameter
#   fixed_params <- list(
#     n = n,
#     d = d,
#     sig_Z = sig_Z,
#     B = B,
#     seed = 4,
#     n_processors = 5
#   )
# 
#   # generate data function
#   generate_data_f <- function(n, d, B, sig_Z) {
#     # in order to compare different methods, we expect to use the same data matrix
#     set.seed(1)
#     no_act <- 5
#     nu <- numeric(d)
#     nu[1:no_act] <- 2*rbinom(no_act, 1, 0.5) - 1
#     beta <- 0.3 * nu # note that beta = gamma = 0.3*nu
#     gamma <- 0.3 * nu
# 
#     # sample Z from multivariate normal
#     Z <- katlabutils::fast_generate_mvn(mean = numeric(d), covariance = sig_Z, num_samples = B * n)
# 
#     # sample X|Z, Y|Z from linear regression model
#     X <- katlabutils::generate_glm_response_data(Z, gamma, "gaussian")
#     Y <- katlabutils::generate_glm_response_data(Z, beta, "gaussian")
# 
#     # generate ground truth mean vectors for oracle methods and debugging
#     E_X_given_Z_oracle <- Z %*% gamma
#     E_Y_given_Z_oracle <- Z %*% beta
#     Var_X_given_Z_oracle <- rep(1, n)
# 
#     # split the data into B lists
#     data_list <- sapply(seq(1, B), function(i) {
#       start <- (i - 1) * n + 1
#       end <- i * n
#       df <- list(
#         X = X[start:end],
#         Y = Y[start:end],
#         Z = Z[start:end, ],
#         Var_X_given_Z_oracle = Var_X_given_Z_oracle,
#         E_X_given_Z_oracle = E_X_given_Z_oracle[start:end],
#         E_Y_given_Z_oracle = E_Y_given_Z_oracle[start:end]
#       )
#       return(df)
#     }, simplify = FALSE)
# 
#     # return the list of datasets
#     return(data_list)
#   }
# 
#   # translate data generation function into a simulatr function
#   generate_data_spec_f <- simulatr::simulatr_function(
#     f = generate_data_f,
#     arg_names = formalArgs(generate_data_f),
#     loop = FALSE
#   )
# 
# 
#   methods_df <- tibble::tribble(
#     ~test_type,
#     ~X_on_Z_reg,
#     ~Y_on_Z_reg,
#     ~test_hyperparams,
#     ##############################
#     "GCM",
#     list(mean_method_type = "LASSO",
#          mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
#          var_method_type = "squared_residual"),
#     list(mean_method_type = "LASSO",
#          mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
#          var_method_type = "squared_residual"),
#     NULL,
#     ##############################
#     "MX2_F_test",
#     list(mean_method_type = "MLE",
#          mean_method_hyperparams = list(family = "gaussian"),
#          var_method_type = "homoskedastic"),
#     list(mean_method_type = "MLE",
#          mean_method_hyperparams = list(family = "gaussian"),
#          var_method_type = "homoskedastic"),
#     NULL
#   )
# 
# 
#   # translate method list into a list of simulatr functions
#   method_list <- generate_method_list(methods_df)
# 
#   # create simulatr specifier object
#   sim_spec <- simulatr::simulatr_specifier(
#     parameter_grid = p_grid,
#     fixed_parameters = fixed_params,
#     generate_data_function = generate_data_spec_f,
#     run_method_functions = method_list
#   )
# 
#   output <- simulatr::check_simulatr_specifier_object(
#     simulatr_spec = sim_spec,
#     B = B,
#     parallel = FALSE
#   )
# 
#   data <- generate_data_f(n = n,
#                           d = d,
#                           B = B,
#                           sig_Z = sig_Z)
# 
#   # specify the argument of X_on_Z_reg/Y_on_Z_reg for GCM
#   X_on_Z_reg_GCM <- list(mean_method_type = "LASSO",
#                          mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
#                          var_method_type = "squared_residual")
#   Y_on_Z_reg_GCM <- list(mean_method_type = "LASSO",
#                          mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
#                          var_method_type = "squared_residual")
# 
#   # specify the argument of X_on_Z_reg/Y_on_Z_reg for MX2
#   X_on_Z_reg_MX2 <- list(mean_method_type = "MLE",
#                          mean_method_hyperparams = list(family = "gaussian"),
#                          var_method_type = "homoskedastic")
#   Y_on_Z_reg_MX2 <- list(mean_method_type = "MLE",
#                          mean_method_hyperparams = list(family = "gaussian"),
#                          var_method_type = "homoskedastic")
# 
#   true_output_GCM <- GCM(data = data[[1]],
#                                  X_on_Z_reg = X_on_Z_reg_GCM,
#                                  Y_on_Z_reg = Y_on_Z_reg_GCM,
#                                  test_hyperparams = NULL)
#   true_output_MX2 <- GCM(data = data[[1]],
#                                  X_on_Z_reg = X_on_Z_reg_MX2,
#                                  Y_on_Z_reg = Y_on_Z_reg_MX2,
#                                  test_hyperparams = NULL)
#   true_output <- rbind(true_output_GCM, true_output_MX2)
#   expect_lt(
#     max(abs(true_output$value - output$results$value)),
#     0.1
#   )
# })


test_that("magnitude_detect works", {
  # check whether the correct magnitude can be detected
  set.seed(1)
  n <- 250
  B <- 200*n
  p <- 500
  beta <- numeric(p)
  # magnitude <- 5
  magnitude <- 0.3
  beta[1:5] <- magnitude*(2*rbinom(5, 1, 0.5) - 1)
  gamma <- beta
  rho <- 0.5
  sig_Z <- katlabutils::generate_cov_ar1(rho = rho, d = p)
  Z <- katlabutils::fast_generate_mvn(mean = numeric(p),
                                 covariance = sig_Z,
                                 num_samples = B)
  res_X_Z <- rnorm(B)
  res_Y_Z <- rnorm(B)
  X <- Z %*% gamma + res_X_Z
  Y <- Z %*% beta + res_Y_Z
  data <- list(res_X_Z = res_X_Z, res_Y_Z = res_Y_Z, Z = Z)
  c <- simulate_confounding(n, X, Y)
  k_true = 5
  magnitude_linsearch <- magnitude_detect(n = n, 
                                          data = data,
                                          c = c,
                                          alpha = 0.05,
                                          beta = beta/k_true,
                                          gamma = gamma/k_true,
                                          eps = 0.2)
  expect_lt(
    max(abs(magnitude_linsearch - k_true)),
    0.2
  )
})




test_that("magnitude_detect works for Binary data", {
  # check whether the correct magnitude can be detected
  set.seed(1)
  n <- 250
  B <- 200*n
  p <- 500
  beta <- numeric(p)
  magnitude <- 5
  beta[1:5] <- magnitude*(2*rbinom(5, 1, 0.5) - 1)
  gamma <- beta
  rho <- 0.5
  sig_Z <- katlabutils::generate_cov_ar1(rho = rho, d = p)
  Z <- katlabutils::fast_generate_mvn(mean = numeric(p), 
                                 covariance = sig_Z, 
                                 num_samples = B)
  X <- rbinom(B, 1, exp(Z %*% gamma)/(1+exp(Z %*% gamma))) 
  Y <- rbinom(B, 1, exp(Z %*% beta)/(1+exp(Z %*% beta)))
  data <- list(Z = Z)
  c <- simulate_confounding(n, X, Y)
  magnitude_linsearch <- magnitude_detect(data = data,
                                                  n = n,
                                                  c = c,
                                                  alpha = 0.01,
                                                  beta = beta/5,
                                                  gamma = gamma/5,
                                                  eps = 0.1,
                                                  type = "binomial")
  expect_lt(
    max(abs(magnitude_linsearch - magnitude)),
    0.4
  )
})




test_that("compute_nu works for Gaussian data", {
  # check whether the correct magnitude can be detected
  set.seed(1)
  n <- 50
  B <- 10000
  d <- 500
  rho <- 0.5
  s <- 5
  c <- seq(0,4,0.2)
  grid <- data.frame(n = n, d = d, s = s, rho = rho)
  grid_new <- symcrt::compute_nu(grid = grid, 
                                 c = c, 
                                 B = B, 
                                 no_nu_grid = 2000, 
                                 response_type = "gaussian")
  beta_base <- numeric(d)
  beta_base[grid_new$coef_pos[[1]]] <- 1
  beta_base[grid_new$coef_neg[[1]]] <- -1
  c_recover <- numeric(length(c))
  # generate covariance matrix and design matrix
  sig <- katlabutils::generate_cov_ar1(rho, d)
  Z <- katlabutils::fast_generate_mvn(mean = numeric(d), 
                                      covariance = sig, 
                                      num_samples = B)
  for (k in 1:length(c)) {
    beta <- beta_base*grid_new$nu[k]
    gamma <- beta
    X <- Z%*%gamma + rnorm(B)
    Y <- Z%*%beta + rnorm(B)
    c_recover[k] <- symcrt::simulate_confounding(n, X, Y)
  }
  expect_lt(
    max(abs(c - c_recover)),
    0.3
  )
})



