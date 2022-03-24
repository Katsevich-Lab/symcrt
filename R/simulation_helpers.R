######################################################################
#
# Helpers functions for numerical simulations.
#
######################################################################

# generate dxd AR(1) covariance matrix with parameter rho

#' Generate dxd AR(1) covariance matrix with parameter rho
#'
#' @param rho
#' @param d
#'
#' @return A dxd covariance matrix.
#' @export
generate_cov_ar1 = function(rho, d){
  sig <- matrix(0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (j in 1:d) {
      sig[i,j] <- rho**(abs(i-j))
    }
  }
  sig
}
