test_that("pGreaterGrades() output", {

  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  det(S)
  solve(S)

  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    val <- latent_distr(lat[1], lat[2], atanh(rho), log(sig), FALSE)
    Rval <- mvtnorm::dmvnorm(lat, sigma = S, log = F)
    expect_true(abs(val-Rval)<1e-6)
  }
})

test_that("check lat_distr() log-likelihood", {

  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  L <- t(chol(S))
  pars <- c(L[2,1], L[2,2])

  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    ll <- lat_distr(lat[1], lat[2], pars)$ll
    Rval <- mvtnorm::dmvnorm(lat, sigma = S, log = TRUE)
    expect_true(abs(ll-Rval)<1e-6)
  }
})

test_that("check lat_distr() gradient", {

  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  L <- t(chol(S))
  pars <- c(L[2,1], L[2,2])

  rfun <- function(PAR, LAT){
    lat_distr(LAT[1], LAT[2], PAR)$ll
  }
  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    grll <- lat_distr(lat[1], lat[2], pars)$gr
    Rval <- numDeriv::grad(func = rfun, x = pars, LAT = lat)
    expect_equal(grll, Rval)
  }
})
