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
    c(d,d)
  )
})

test_that("fast_generate_mvn works",{
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
  expect_lt(max(abs(sample_cov - true_cov)),
                   0.1)
  expect_equal(dim(mvn_data), c(num_samples, d))
})
