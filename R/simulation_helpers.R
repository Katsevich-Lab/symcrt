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
      "%s_%s_%s_%d",
      test_type,
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
#' @param data A data list containing the noise vector for X|Z and Y|Z and data matrix Z
#' @param c The target confounding level
#' @param alpha The step size when doing line search
#' @param beta The coefficient vector for Y|Z
#' @param gamma The coefficient vector for X|Z
#' @param eps The precision number
#'
#' @return The magnitude of coefficient reaches the confounding level
#' @export
magnitude_detect <- function(data, c, alpha, beta, gamma, eps = 0.0001) {
  res_X_Z <- data$res_X_Z
  res_Y_Z <- data$res_Y_Z
  Z <- data$Z
  base_confoun <- simulate_confounding(X = Z %*% gamma, Y = Z %*% beta)
  if(base_confoun*c<0){
    stop("The sign of target confounding does not match with that of base line!")
  }
  i <- 1
  confoun_level <- 0
  while (abs(confoun_level-c) > eps) {
    kappa <- alpha*i
    X <- kappa*Z %*% gamma + res_X_Z
    Y <- kappa*Z %*% beta + res_Y_Z
    confoun_level <- simulate_confounding(X, Y)
    i <- i + 1
    if (i > 1E10){
      stop("Exceed the maximum iteration!")
    }
  }
  kappa
}

#' A function for calculating the confoudning level via simulation
#'
#' @param X A vector of treatment
#' @param Y A vector of outcome
#'
#' @return The confounding level
#' @export
simulate_confounding <- function(X, Y){
  n <- length(X)
  sqrt(n)*mean(X*Y)/stats::sd(X*Y)
}

