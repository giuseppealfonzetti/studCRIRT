test_that("pGreaterGrades() output", {

  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  det(S)
  solve(S)

  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    val <- latent_distr(lat[1], lat[2], rho, sig, FALSE)
    Rval <- mvtnorm::dmvnorm(lat, sigma = S, log = F)
    expect_true(abs(val-Rval)<1e-6)
  }
})
