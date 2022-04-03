test_that("generate_cov_ar1 works", {

  # test whether a specific example gives exactly the right answer
  # (generated using dput)
  rho <- 0.5
  d <- 5
  expect_identical(
    generate_cov_ar1(rho, d),
    structure(c(
      1, 0.5, 0.25, 0.125, 0.0625, 0.5, 1, 0.5, 0.25, 0.125,
      0.25, 0.5, 1, 0.5, 0.25, 0.125, 0.25, 0.5, 1, 0.5, 0.0625, 0.125,
      0.25, 0.5, 1
    ), .Dim = c(5L, 5L))
  )

  # test whether we get identity matrix for rho = 0
  rho <- 0
  d <- 100
  expect_identical(
    generate_cov_ar1(rho, d),
    diag(d)
  )

  # test whether dimensions of output are correct
  rho <- 0.2
  d <- 53
  expect_equal(
    dim(generate_cov_ar1(rho, d)),
    c(d, d)
  )
})

test_that("fast_generate_mvn works", {
  # check whether the empirical covariance roughly matches the true covariance
  # for a large number of samples from a MVN
  rho <- 0.3
  d <- 5
  num_samples <- 1000
  mean <- numeric(d)
  true_cov <- generate_cov_ar1(rho, d)
  set.seed(1)
  mvn_data <- fast_generate_mvn(mean, true_cov, num_samples)
  sample_cov <- var(mvn_data)
  expect_lt(
    max(abs(sample_cov - true_cov)),
    0.1
  )
  expect_equal(dim(mvn_data), c(num_samples, d))
})


test_that("generate_method_list works", {
  # check whether function works with different methods input
  n <- 200
  d <- 150
  no_act <- 20
  B <- 1
  sig_Z <- diag(d)
  
  # fake p_grid
  p_grid <- expand.grid(b = c(1)) 
  
  # fixed parameter
  fixed_params <- list(
    n = n,
    d = d,
    sig_Z = sig_Z,
    B = B,
    seed = 4,
    n_processors = 5
  )
  
  # generate data function
  generate_data_f <- function(n, d, B, sig_Z) {
    # in order to compare different methods, we expect to use the same data matrix
    set.seed(1)
    no_act <- 5
    nu <- numeric(d)
    nu[1:no_act] <- 2*rbinom(no_act, 1, 0.5) - 1
    beta <- 0.3 * nu # note that beta = gamma = 0.3*nu
    gamma <- 0.3 * nu
    
    # sample Z from multivariate normal
    Z <- symcrt::fast_generate_mvn(mean = numeric(d), covariance = sig_Z, num_samples = B * n)
    
    # sample X|Z, Y|Z from linear regression model
    X <- symcrt::generate_glm_response_data(Z, gamma, "gaussian")
    Y <- symcrt::generate_glm_response_data(Z, beta, "gaussian")
    
    # generate ground truth mean vectors for oracle methods and debugging
    E_X_given_Z_oracle <- Z %*% gamma
    E_Y_given_Z_oracle <- Z %*% beta
    Var_X_given_Z_oracle <- rep(1, n)
    
    # split the data into B lists
    data_list <- sapply(seq(1, B), function(i) {
      start <- (i - 1) * n + 1
      end <- i * n
      df <- list(
        X = X[start:end],
        Y = Y[start:end],
        Z = Z[start:end, ],
        Var_X_given_Z_oracle = Var_X_given_Z_oracle,
        E_X_given_Z_oracle = E_X_given_Z_oracle[start:end],
        E_Y_given_Z_oracle = E_Y_given_Z_oracle[start:end]
      )
      return(df)
    }, simplify = FALSE)
    
    # return the list of datasets
    return(data_list)
  }
  
  # translate data generation function into a simulatr function
  generate_data_spec_f <- simulatr::simulatr_function(
    f = generate_data_f,
    arg_names = formalArgs(generate_data_f),
    loop = FALSE
  )
  
  
  methods_df <- tibble::tribble(
    ~test_type,
    ~X_on_Z_reg,
    ~Y_on_Z_reg,
    ~test_hyperparams,
    ##############################
    "GCM",
    list(mean_method_type = "LASSO",
         mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
         var_method_type = "squared_residual"),
    list(mean_method_type = "LASSO",
         mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
         var_method_type = "squared_residual"),
    NULL,
    ##############################
    "MX2_F_test",
    list(mean_method_type = "MLE",
         mean_method_hyperparams = list(family = "gaussian"),
         var_method_type = "homoskedastic"),
    list(mean_method_type = "MLE",
         mean_method_hyperparams = list(family = "gaussian"),
         var_method_type = "homoskedastic"),
    NULL
  )
  
  
  # translate method list into a list of simulatr functions
  method_list <- symcrt::generate_method_list(methods_df)
  
  # create simulatr specifier object
  sim_spec <- simulatr::simulatr_specifier(
    parameter_grid = p_grid,
    fixed_parameters = fixed_params,
    generate_data_function = generate_data_spec_f,
    run_method_functions = method_list
  )
  
  output <- simulatr::check_simulatr_specifier_object(
    simulatr_spec = sim_spec,
    B = B,
    parallel = FALSE
  )
  
  data <- generate_data_f(n = n, 
                          d = d, 
                          B = B,  
                          sig_Z = sig_Z)
  
  # specify the argument of X_on_Z_reg/Y_on_Z_reg for GCM
  X_on_Z_reg_GCM <- list(mean_method_type = "LASSO",
                         mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
                         var_method_type = "squared_residual")
  Y_on_Z_reg_GCM <- list(mean_method_type = "LASSO",
                         mean_method_hyperparams = list(s = "lambda.min", nfolds = 10, family = "gaussian"),
                         var_method_type = "squared_residual")
  
  # specify the argument of X_on_Z_reg/Y_on_Z_reg for MX2
  X_on_Z_reg_MX2 <- list(mean_method_type = "MLE",
                         mean_method_hyperparams = list(family = "gaussian"),
                         var_method_type = "homoskedastic")
  Y_on_Z_reg_MX2 <- list(mean_method_type = "MLE",
                         mean_method_hyperparams = list(family = "gaussian"),
                         var_method_type = "homoskedastic")
  
  true_output_GCM <- symcrt::GCM(data = data[[1]],
                                 X_on_Z_reg = X_on_Z_reg_GCM,
                                 Y_on_Z_reg = Y_on_Z_reg_GCM,
                                 test_hyperparams = NULL)
  true_output_MX2 <- symcrt::GCM(data = data[[1]],
                                 X_on_Z_reg = X_on_Z_reg_MX2,
                                 Y_on_Z_reg = Y_on_Z_reg_MX2,
                                 test_hyperparams = NULL)
  true_output <- rbind(true_output_GCM, true_output_MX2)
  expect_lt(
    max(abs(true_output$value - output$results$value)),
    0.1
  )
})

