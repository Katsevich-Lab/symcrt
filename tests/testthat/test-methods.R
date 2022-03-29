test_that("dCRT_f works", {
  
  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  n <- 100
  d <- 5
  gamma <- 2*rbinom(d, 1, 0.5) - 1
  Z <- matrix(rnorm(n = n*d), nrow = n, ncol = d)
  noise <- rnorm(2*n)
  cond_mean <- Z %*% gamma |> as.vector()
  cond_var <- rep(1, n)
  X <- cond_mean + noise[1:n]
  Y <- Z %*% gamma + noise[(n+1):(2*n)]
  data <- list(Z = Z, 
               X = X, 
               Y = Y,
               cond_mean = cond_mean,
               cond_var = cond_var)
  X_on_Z_reg <- "oracle"
  Y_on_Z_reg <- "OLS"
  X_on_Z_var <- "oracle"
  resample_dist <- "Gaussian"
  no_resample <- 10000
  # test whether we get similar p-value from MX2 and dCRT
  # run MX2
  MX2_result <- MX2_F_test_f(data = data,
                              X_on_Z_reg = X_on_Z_reg,
                              Y_on_Z_reg = Y_on_Z_reg,
                              X_on_Z_var = X_on_Z_var)
  p_value_MX2 <- MX2_result |> 
                 dplyr::filter(parameter == "p_value") |> 
                 dplyr::select(value) |> unname() |> as.numeric()
  # run dCRT
  dCRT_result <- dCRT_f(data = data,
                        X_on_Z_reg = X_on_Z_reg,
                        Y_on_Z_reg = Y_on_Z_reg,
                        X_on_Z_var = X_on_Z_var,
                        resample_dist = resample_dist,
                        no_resample = no_resample)
  p_value_dCRT <- dCRT_result |> 
                  dplyr::select(value) |> 
                  unname() |> as.numeric()
  expect_lt(abs(p_value_MX2 - p_value_dCRT),
            0.01)
})
