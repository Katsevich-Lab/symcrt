######################################################################
#
# Helpers functions for numerical simulations.
#
######################################################################

#' Helper function to translate the data frame of method specifications into a
#' list of simulatr functions.
#'
#' @param methods_df A data frame of the kind we define in sim_setting_x.R
#'
#' @return A list of simulatr functions whose length is the number of rows of \code{methods_df}.
#' @export
generate_method_list <- function(methods_df) {
  num_methods <- nrow(methods_df)
  simulatr_functs <- list()
  # loop through methods, which are rows of methods_df
  for (method_idx in 1:num_methods) {
    # extract the four pieces of information about each method
    test_type <- methods_df$test_type[method_idx]
    X_on_Z_reg <- methods_df$X_on_Z_reg[[method_idx]]
    Y_on_Z_reg <- methods_df$Y_on_Z_reg[[method_idx]]
    test_hyperparams <- methods_df$test_hyperparams[[method_idx]]
    test_hyperparams <- set_default_test_hyperparams(test_type, test_hyperparams)

    # test_type is a function in the simulatr package (e.g., GCM), and we want to
    # prepend "simulatr::" to this function, so we need to do it in this fancy way
    test_type_package <- eval(parse(text = sprintf("symcrt::%s", test_type)))
    # the method function takes just data as an argument, and call the function
    # test_type (e.g. GCM), after prepending "simulatr::" on data as well as the
    # three pieces of information about the method (X|Z reg, Y|Z reg, hyperparams)
    method_f <-
      local({
        test_type_package <- test_type_package
        X_on_Z_reg <- X_on_Z_reg
        Y_on_Z_reg <- Y_on_Z_reg
        test_hyperparams <- test_hyperparams
        function(data) {
          do.call(
            test_type_package,
            list(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams)
          )
        }
      })
    # we also create a name for each method, based on the test type, X|Z method,
    # Y|Z method, and then the row of the methods data frame, the latter to
    # distinguish between the same method with different hyperparameters
    method_name <- sprintf(
      "%s_%s_%s_%s_%d",
      test_type,
      test_hyperparams$way_to_learn,
      X_on_Z_reg$mean_method_type,
      Y_on_Z_reg$mean_method_type,
      method_idx
    )

    # define the simulatr function, specifying arg_names = NA_character_, which
    # means that there are no arguments to the method function besides data
    simulatr_functs[[method_name]] <- simulatr::simulatr_function(
      f = method_f,
      arg_names = NA_character_,
      loop = TRUE
    )
  }
  simulatr_functs
}

#' A line search algorithm for detecting the magnitude of coefficient reaching the confounding level
#'
#' @param n The number of data 
#' @param data A data list containing the noise vector for X|Z and Y|Z and data matrix Z
#' @param c The target confounding level
#' @param alpha The step size when doing line search
#' @param beta The coefficient vector for Y|Z
#' @param gamma The coefficient vector for X|Z
#' @param eps The precision number
#' @param type The type of response (either Gaussian or Binary)
#'
#' @return The magnitude of coefficient reaches the confounding level
#' @export
magnitude_detect <- function(n, data, c, alpha, beta, gamma, eps = 0.0001, type = "Gaussian") {
  Z <- data$Z
  B <- nrow(Z)
  confoun_level <- 0
  response_type <- type
  predictor.X <- Z %*% gamma
  predictor.Y <- Z %*% beta
  switch(response_type,
         Gaussian = {
           base_confoun <- simulate_confounding(n, 
                                                X = predictor.X + stats::rnorm(B), 
                                                Y = predictor.Y + stats::rnorm(B))
           if(base_confoun*c<0){
             stop("The sign of target confounding does not match with that of base line!")
           }
           i <- 1
           while (abs(confoun_level-c) > eps) {
             kappa <- alpha*i
             X <- kappa*predictor.X + stats::rnorm(B)
             Y <- kappa*predictor.Y + stats::rnorm(B)
             confoun_level <- simulate_confounding(n, X, Y)
             i <- i + 1
             if (i > 1E10){
               stop("Exceed the maximum iteration!")
             }
           }
         },
         Binary = {
           X_base <- stats::rbinom(B, 1, exp(predictor.X)/(1+exp(predictor.X)))
           Y_base <- stats::rbinom(B, 1, exp(predictor.Y)/(1+exp(predictor.Y)))
           base_confoun <- simulate_confounding(n, 
                                                X = X_base, 
                                                Y = Y_base)
           if(base_confoun*c<0){
             stop("The sign of target confounding does not match with that of base line!")
           }
           i <- 1
           while (abs(confoun_level-c) > eps) {
             kappa <- alpha*i
             X <- stats::rbinom(B, 1, exp(kappa*predictor.X)/(1+exp(kappa*predictor.X)))
             Y <- stats::rbinom(B, 1, exp(kappa*predictor.Y)/(1+exp(kappa*predictor.Y)))
             confoun_level <- simulate_confounding(n, X, Y)
             i <- i + 1
             if (i > 1E10){
               stop("Exceed the maximum iteration!")
             }
           }
         },
         {
           stop("Invalid specification of response type.")
         }
  )
  kappa
}

#' A function for calculating the confoudning level via simulation
#'
#' @param n number of samples
#' @param X A vector of treatment (length may be much larger than n)
#' @param Y A vector of outcome (length may be much larger than n)
#'
#' @return The confounding level
#' @export
simulate_confounding <- function(n, X, Y){
  sqrt(n)*mean((X-mean(X))*(Y-mean(Y)))/stats::sd((X-mean(X))*(Y-mean(Y)))
}

#' A function for calculating the approximation power for oracle GCM
#'
#' @param n number of samples
#' @param X A vector of treatment (length may be much larger than n)
#' @param Y A vector of outcome (length may be much larger than n)
#' @param E_X_given_Z A vector of conditional expectation of X|Z
#' @param E_Y_given_Z A vector of conditional expectation of Y|Z
#'
#' @return A number representing the power
#' @export

simulate_power <- function(n, X, Y, E_X_given_Z, E_Y_given_Z){
  sqrt(n)*mean((X - E_X_given_Z)*(Y - E_Y_given_Z))/stats::sd((X - E_X_given_Z)*(Y - E_Y_given_Z))
}


#' Title
#'
#' @param grid A data frame containing different combinations of parameters
#' @param c A vector containing a sequence of confounding level. Default is 1:4
#' @param B A large number used for simulating confounding level
#' @param no_nu_grid The number of grid that spread from 0 to the maximum confounding
#' level which is searched by the algorithm
#' @param response_type Type of response, Gaussian or Binary
#'
#' @return A data frame containing extra two parameters: confounding level and magnitude.
#' @export

compute_nu <- function(grid, c, B, no_nu_grid, response_type){
  no_grid <- nrow(grid)
  nu_mat <- matrix(0, nrow = length(c), ncol = no_grid)
  grid$coef_pos <- 0
  grid$coef_neg <- 0
  for (r in 1:no_grid) {
    # extract key variables
    n <- grid$n[r]
    d <- grid$d[r]
    s <- grid$s[r]
    rho <- grid$rho[r]
    
    # generate covariance matrix
    sig <- katlabutils::generate_cov_ar1(rho, d)
    
    # base beta and gamma
    base.beta <- numeric(d)
    set.seed(s)
    base.beta[1:s] <- 2*stats::rbinom(s, 1, 0.5) - 1
    base.gamma <- base.beta
    
    # generate Z data with B rows
    Z <- katlabutils::fast_generate_mvn(mean = numeric(d), 
                                        covariance = sig, 
                                        num_samples = B)
    # roughly search the nu corresponding to c_max
    nu_search <- 0
    predictor.X <- Z %*% base.gamma
    predictor.Y <- Z %*% base.beta 
    c_max <- max(c)
    switch(response_type,
           Gaussian = {
             base_confoun <- symcrt::simulate_confounding(n, 
                                                  X = predictor.X + stats::rnorm(B), 
                                                  Y = predictor.Y + stats::rnorm(B))
             if(base_confoun*c_max<0){
               stop("The sign of target confounding does not match with that of base line!")
             }
             c_search <- 0
             while (c_search-c_max < 0) {
               nu_search <- nu_search + 0.1
               
               X <- nu_search*predictor.X + stats::rnorm(B)
               Y <- nu_search*predictor.Y + stats::rnorm(B)
               
               # compute the actual level c based on X, Y
               c_search <- symcrt::simulate_confounding(n, X, Y)
             }
           },
           Binary = {
             X_base <- stats::rbinom(B, 1, 1/(1+exp(-predictor.X)))
             Y_base <- stats::rbinom(B, 1, 1/(1+exp(-predictor.Y)))
             base_confoun <- simulate_confounding(n, 
                                                  X = X_base, 
                                                  Y = Y_base)
             if(base_confoun*c_max<0){
               stop("The sign of target confounding does not match with that of base line!")
             }
             c_search <- 0
             while (c_search-c_max < 0) {
               nu_search <- nu_search + 0.1
               X <- stats::rbinom(n = B, size = 1, prob = 1/(1+exp(-nu_search*predictor.X)))
               Y <- stats::rbinom(n = B, size = 1, prob = 1/(1+exp(-nu_search*predictor.Y)))
               
               # compute the actual level c based on X, Y
               c_search <- symcrt::simulate_confounding(n, X, Y)
             }
           },
           {
             stop("Invalid specification of response type.")
           }
    )
    
    # generate a sequence of confounding level
    nu_max <- nu_search
    nu_seq <- seq(0, nu_max, length.out = no_nu_grid)
    c_seq <- numeric(length(nu_seq))
    switch(response_type,
           Gaussian = {
             for (k in 1:length(nu_seq)) {
               # generate X and Y
               X <- nu_seq[k]*predictor.X + stats::rnorm(B)
               Y <- nu_seq[k]*predictor.Y + stats::rnorm(B)
               # generate corresponding confounding level
               c_seq[k] <- symcrt::simulate_confounding(n, X, Y)
             }
           },
           Binary = {
             for (k in 1:length(nu_seq)) {
               # generate X and Y
               X <- stats::rbinom(n = B, size = 1, prob = 1/(1+exp(-nu_seq[k]*predictor.X)))
               Y <- stats::rbinom(n = B, size = 1, prob = 1/(1+exp(-nu_seq[k]*predictor.Y)))
               # generate corresponding confounding level
               c_seq[k] <- symcrt::simulate_confounding(n, X, Y)
             }
           },
           {
             stop("Invalid specification of response type.")
           }
    )
    # fit the (c_seq, nu_seq) curve with loess
    nu_seq <- c(numeric(no_nu_grid), nu_seq)
    c_seq <- c(numeric(no_nu_grid), c_seq)
    c_nu <- data.frame(nu = nu_seq, c = c_seq)
    poly_fit <- stats::loess(nu ~ c, c_nu)
    nu_mat[,r] <- as.vector(stats::predict(poly_fit, data.frame(c = c), se = FALSE))
    # force nu=0 when c =0
    if(0 %in% c){
      pos <- which(c == 0)
      nu_mat[pos, r] <- 0
    }
    grid$coef_pos[r] <- list(which(base.beta > 0))
    grid$coef_neg[r] <- list(setdiff(1:s, which(base.beta > 0)))
  }
  
  # replicate the grid to accomodate c and nu
  grid <- grid[rep(1:no_grid, each = length(c)),] 
  
  # fill the last two columns with c and nu
  grid$c <- rep(c, no_grid)
  grid$nu <- c(nu_mat)
  grid$grid_id <- 1:nrow(grid)
  grid
}

#' Compute the sequence of nu according to presepcfied level of type_I error
#'
#' @param maxType_I A maximum value of type_I error we want to reach
#' @param alpha Level of test
#' @param test_type One side or Two side
#' @param grid The original parameter grid without nu
#' @param no_nu Number of nu we want
#' @param B Monte Carlo times
#' @param response_type Gaussian or Binary
#'
#' @return A new parameter grid with extra column of nu
#' @export

compute_nu_via_Type_I <- function(maxType_I, alpha, test_type, grid, no_nu, B, response_type){
  no_grid <- nrow(grid)
  nu_mat <- matrix(0, nrow = no_nu, ncol = no_grid)
  grid$coef_pos <- 0
  grid$coef_neg <- 0
  for (r in 1:no_grid) {
    # extract key variables
    n <- grid$n[r]
    d <- grid$d[r]
    s <- grid$s[r]
    rho <- grid$rho[r]
    
    # generate covariance matrix
    sig <- katlabutils::generate_cov_ar1(rho, d)
    
    # base beta and gamma
    base.beta <- numeric(d)
    set.seed(s)
    base.beta[1:s] <- 2*stats::rbinom(s, 1, 0.5) - 1
    base.gamma <- base.beta
    
    # generate Z data with B rows
    Z <- katlabutils::fast_generate_mvn(mean = numeric(d), 
                                        covariance = sig, 
                                        num_samples = B)
    # roughly search the nu corresponding to maxType_I
    nu_search <- 0
    predictor.X <- Z %*% base.gamma
    predictor.Y <- Z %*% base.beta 
    switch(response_type,
           Gaussian = {
             base_confoun <- symcrt::simulate_confounding(n, 
                                                          X = predictor.X + stats::rnorm(B), 
                                                          Y = predictor.Y + stats::rnorm(B))
             type_I_search <- 0
             while (type_I_search-maxType_I < 0) {
               nu_search <- nu_search + sign(base_confoun)*0.02
               
               X <- nu_search*predictor.X + stats::rnorm(B)
               Y <- nu_search*predictor.Y + stats::rnorm(B)
               
               # compute the actual level c and type_I error based on X, Y
               c_search <- symcrt::simulate_confounding(n, X, Y)
               type_I_search <- symcrt::pval_shift(alpha = alpha, 
                                                   c = c_search,
                                                   type = test_type)
             }
           },
           Binary = {
             X_base <- stats::rbinom(B, 1, 1/(1+exp(-predictor.X)))
             Y_base <- stats::rbinom(B, 1, 1/(1+exp(-predictor.Y)))
             base_confoun <- simulate_confounding(n, 
                                                  X = X_base, 
                                                  Y = Y_base)
             type_I_search <- 0
             while (type_I_search-maxType_I < 0) {
               nu_search <- nu_search + sign(base_confoun)*0.02
               X <- stats::rbinom(n = B, size = 1, prob = 1/(1+exp(-nu_search*predictor.X)))
               Y <- stats::rbinom(n = B, size = 1, prob = 1/(1+exp(-nu_search*predictor.Y)))
               
               # compute the actual level c and type_I error based on X, Y
               c_search <- symcrt::simulate_confounding(n, X, Y)
               type_I_search <- symcrt::pval_shift(alpha = alpha, 
                                                   c = c_search,
                                                   type = test_type)
             }
           },
           {
             stop("Invalid specification of response type.")
           }
    )
    
    # generate a sequence of confounding level
    nu_max <- nu_search
    nu_mat[,r] <- seq(0, nu_max, nu_max/(no_nu-1))
    
    # store the sign of coefficient
    grid$coef_pos[r] <- list(which(base.beta > 0))
    grid$coef_neg[r] <- list(setdiff(1:s, which(base.beta > 0)))
  }
  
  # replicate the grid to accommodate nu
  grid <- grid[rep(1:no_grid, each = no_nu),] 
  
  # fill the last two columns with nu
  grid$nu <- c(nu_mat)
  grid$grid_id <- 1:nrow(grid)
  grid
}


#' A function to compute a sequence of theta
#'
#' @param maxPower The maximum power prespecified to detect maximum theta
#' @param alpha Level of test
#' @param test_type One-side or Two-side test
#' @param grid A parameter grid
#' @param no_theta Number of theta to obtain
#' @param B Number of MC replications
#' @param response_type Gaussian or Binary
#'
#' @return A new parameter grid
#' @export
#' 
compute_theta_via_power <- function(maxPower, alpha, test_type, grid, no_theta, B, response_type){
  no_grid <- nrow(grid)
  theta_mat <- matrix(0, nrow = no_theta, ncol = no_grid)
  for (r in 1:no_grid) {
    # extract key variables
    n <- grid$n[r]
    d <- grid$d[r]
    s <- grid$s[r]
    rho <- grid$rho[r]
    beta <- numeric(d)
    beta[grid$coef_neg[[r]]] <- (-1)*grid$nu[r]
    beta[grid$coef_pos[[r]]] <- grid$nu[r]
    gamma <- beta
    
    # generate covariance matrix
    sig <- katlabutils::generate_cov_ar1(rho, d)
    
    # base theta
    base.theta <- 1
    
    # generate Z data with B rows
    Z <- katlabutils::fast_generate_mvn(mean = numeric(d), 
                                        covariance = sig, 
                                        num_samples = B)
    # roughly search the nu corresponding to maxPower
    theta_search <- 0
    predictor.X <- Z %*% gamma
    predictor.Y <- predictor.X*base.theta + Z %*% beta 
    switch(response_type,
           Gaussian = {
             X <- predictor.X + stats::rnorm(B)
             Y <- X*base.theta + Z %*% beta + stats::rnorm(B)
             E_X_given_Z <- Z %*% gamma
             E_Y_given_Z <- E_X_given_Z*base.theta + Z %*% beta
             base_power <- symcrt::simulate_power(n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
             power_search <- 0
             while (power_search-maxPower < 0) {
               theta_search <- theta_search + sign(base_power)*0.02
               
               # update the conditional expectation of E_Y_given_Z
               X <- predictor.X + stats::rnorm(B)
               Y <- X*theta_search + Z %*% beta + stats::rnorm(B)
               E_Y_given_Z <- E_X_given_Z*theta_search + Z %*% beta
               
               # compute the actual level c and power based on X, Y
               c_search <- symcrt::simulate_power(n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
               power_search <- symcrt::pval_shift(alpha = alpha, 
                                                  c = c_search,
                                                  type = test_type)
             }
           },
           Binary = {
             E_X_given_Z <- symcrt::glogit(Z %*% gamma)
             X <- stats::rbinom(B, 1, E_X_given_Z)
             E_Y_given_Z <- glogit(base.theta + Z %*% beta)*glogit(Z %*% gamma) + 
               glogit(Z %*% beta)*(1 - glogit(Z %*% gamma))
             Y <- stats::rbinom(B, 1, glogit(X*base.theta + Z %*% beta ))
             base_power <- symcrt::simulate_power(n = n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
             power_search <- 0
             while (power_search-maxPower < 0) {
               theta_search <- theta_search + sign(base_power)*0.01
               
               # update the conditional expectation
               X <- stats::rbinom(n = B, size = 1, prob = E_X_given_Z)
               E_Y_given_Z <- glogit(theta_search + Z %*% beta)*glogit(Z %*% beta) + 
                 glogit(Z %*% beta)*(1 - glogit(Z %*% beta))
               Y <- stats::rbinom(n = B, size = 1, prob = glogit(X*theta_search + Z %*% beta))
               
               # compute the actual level c and power based on X, Y
               c_search <- symcrt::simulate_power(n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
               power_search <- symcrt::pval_shift(alpha = alpha, 
                                                  c = c_search,
                                                  type = test_type)
             }
           },
           {
             stop("Invalid specification of response type.")
           }
    )
    
    # generate a sequence of theta
    theta_max <- theta_search
    theta_mat[,r] <- seq(0, theta_max, theta_max/(no_theta-1))
  }
  
  # replicate the grid to accommodate theta
  grid <- grid[rep(1:no_grid, each = no_theta),] 
  
  # fill the last two columns with theta
  grid$theta <- c(theta_mat)
  grid$grid_id <- 1:nrow(grid)
  grid
}


#' Compute a sequence of theta and nu for a given level of confounding level
#'
#' @param maxPower The maximum power prespecified to detect maximum theta
#' @param alpha Level of test
#' @param test_type One-side or Two-side test
#' @param grid A parameter grid
#' @param no_theta Number of theta to obtain
#' @param B Number of MC replications
#' @param response_type Gaussian or Binary
#'
#' @return A new parameter grid
#' @export
#' 
compute_theta_via_power_fixed_c <- function(maxPower, alpha, test_type, grid, no_theta, B, response_type){
  no_grid <- nrow(grid)
  theta_mat <- matrix(0, nrow = no_theta, ncol = no_grid)
  nu_mat <- matrix(0, nrow = no_theta, ncol = no_grid)
  for (r in 1:no_grid) {
    # extract key variables
    n <- grid$n[r]
    d <- grid$d[r]
    s <- grid$s[r]
    rho <- grid$rho[r]
    beta <- numeric(d)
    beta[grid$coef_neg[[r]]] <- (-1)*grid$nu[r]
    beta[grid$coef_pos[[r]]] <- grid$nu[r]
    gamma <- beta
    
    # generate covariance matrix
    sig <- katlabutils::generate_cov_ar1(rho, d)
    
    # base theta
    base.theta <- 0.7
    
    # generate Z data with B rows
    Z <- katlabutils::fast_generate_mvn(mean = numeric(d), 
                                        covariance = sig, 
                                        num_samples = B)
    # roughly search the nu corresponding to maxPower
    theta_search <- 0
    switch(response_type,
           Gaussian = {
             X <- Z %*% gamma + stats::rnorm(B)
             # set base.theta = 1 as before
             beta_star <- beta - base.theta*gamma
             Y <- X*base.theta + Z %*% beta_star + stats::rnorm(B)
             E_X_given_Z <- Z %*% gamma
             E_Y_given_Z <- Z %*% beta
             base_power <- symcrt::simulate_power(n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
             power_search <- 0
             while (power_search-maxPower < 0) {
               theta_search <- theta_search + sign(base_power)*0.02
               
               # update the conditional expectation of E_Y_given_Z
               X <- Z %*% gamma + stats::rnorm(B)
               beta_star <- beta - theta_search*gamma
               Y <- X*theta_search + Z %*% beta_star + stats::rnorm(B)
               E_X_given_Z <- Z %*% gamma
               E_Y_given_Z <- Z %*% beta
               # compute the actual level c and power based on X, Y
               c_search <- symcrt::simulate_power(n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
               power_search <- symcrt::pval_shift(alpha = alpha, 
                                                  c = c_search,
                                                  type = test_type)
             }
           },
           Binary = {
             E_X_given_Z <- symcrt::glogit(Z %*% gamma)
             X <- stats::rbinom(B, 1, E_X_given_Z)
             E_Y_given_Z <- glogit(base.theta + Z %*% beta)*glogit(Z %*% gamma) + 
               glogit(Z %*% beta)*(1 - glogit(Z %*% gamma))
             Y <- stats::rbinom(B, 1, glogit(X*base.theta + Z %*% beta ))
             base_power <- symcrt::simulate_power(n = n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
             power_search <- 0
             while (power_search-maxPower < 0) {
               theta_search <- theta_search + sign(base_power)*0.01
               
               # update the conditional expectation
               X <- stats::rbinom(n = B, size = 1, prob = E_X_given_Z)
               E_Y_given_Z <- glogit(theta_search + Z %*% beta)*glogit(Z %*% beta) + 
                 glogit(Z %*% beta)*(1 - glogit(Z %*% beta))
               Y <- stats::rbinom(n = B, size = 1, prob = glogit(X*theta_search + Z %*% beta))
               
               # compute the actual level c and power based on X, Y
               c_search <- symcrt::simulate_power(n,
                                                  X = X, 
                                                  Y = Y,
                                                  E_X_given_Z = E_X_given_Z,
                                                  E_Y_given_Z = E_Y_given_Z)
               power_search <- symcrt::pval_shift(alpha = alpha, 
                                                  c = c_search,
                                                  type = test_type)
             }
           },
           {
             stop("Invalid specification of response type.")
           }
    )
    
    # generate a sequence of theta
    theta_max <- theta_search
    theta_mat[,r] <- seq(0, theta_max, theta_max/(no_theta-1))
  }
  
  # replicate the grid to accommodate theta
  grid <- grid[rep(1:no_grid, each = no_theta),] 
  
  # fill the last two columns with theta
  grid$theta <- c(theta_mat)
  grid$grid_id <- 1:nrow(grid)
  grid
}


#' Compute the shifted p_value
#'
#' @param alpha Level of test
#' @param c Confounding level
#' @param type Type of test: one-sided or two-sided
#'
#' @return Shifted p_value
#' @export

pval_shift <- function(alpha, c, type = "two_side"){
  if(type == "one_side" & c > 0){
    1 - stats::pnorm(stats::qnorm(1-alpha) - c)
  }else if(type == "one_side" & c < 0){
    stats::pnorm(stats::qnorm(alpha) - c)
  }else if(type == "two_side" & c > 0){
    1 - stats::pnorm(stats::qnorm(1-alpha/2) - c)
  }else if (type == "two_side" & c < 0){
    stats::pnorm(stats::qnorm(alpha/2) - c)
  }
}
