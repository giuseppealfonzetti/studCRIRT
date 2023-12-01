test_that("hazard() output", {
  ## External covariates CR
  dim_ext_cr <- 3L;

  set.seed(1)
  beta_d <- rnorm(2 + dim_ext_cr, 0, 1)
  beta_t <- rnorm(2 + dim_ext_cr, 0, 1)
  beta_g <- rnorm(2 + dim_ext_cr, 0, 1)

  beta_cr <- list('d' = beta_d, 't' = beta_t, 'g' = beta_g)
  ## Number of years
  number_years_before <- 6L
  number_years_after <- 2L

  ## Year-related intercept
  set.seed(2)

  ## year-specific intercepts for `d`, `t` and `g`
  beta0_d <- rnorm(number_years_before, 0, 1)
  beta0_t <- rnorm(number_years_before, 0, 1)
  beta0_g <- rnorm(number_years_after, 0, 1)
  beta0_cr <- list('d' = beta0_d, 't' = beta0_t, 'g' = beta0_g)

  theta_cr <- c(beta_d, beta_t, beta_g, beta0_d, beta0_t, beta0_g)

  ## Exam-related params
  grades <- 4L
  exams <- 3L

  set.seed(3)
  alpha <- rnorm(exams, 0, 1)
  beta_exams <- list()
  for(ex in 1:exams){
    beta_exams[[ex]] <- rnorm(grades, 0, 1)
  }

  theta_irt <- c(alpha, unlist(beta_exams))

  theta <- c(theta_cr, theta_irt)

  set.seed(4)
  x <- rnorm(dim_ext_cr+2, 0, 1)
  for (year in 1:number_years_before) {
    haz_d <- hazard(
      OUTCOME = 1,
      YEAR = year,
      THETA_CR = theta_cr,
      COVARIATES = x,
      NYB = number_years_before,
      NYA = number_years_after
    )

    haz_t <- hazard(
      OUTCOME = 2,
      YEAR = year,
      THETA_CR = theta_cr,
      COVARIATES = x,
      NYB = number_years_before,
      NYA = number_years_after
    )

    expeta_d <- as.numeric(exp(beta0_d[year] + t(beta_d)%*%x))
    expeta_t <- as.numeric(exp(beta0_t[year] + t(beta_t)%*%x))

    expect_equal(haz_d, expeta_d/(1+expeta_d+expeta_t))
    expect_equal(haz_t, expeta_t/(1+expeta_d+expeta_t))
  }

  for (year in 1:number_years_after) {
    haz_g <- hazard(
      OUTCOME = 3,
      YEAR = year,
      THETA_CR = theta_cr,
      COVARIATES = x,
      NYB = number_years_before,
      NYA = number_years_after
    )

    expeta_g <- as.numeric(exp(beta0_g[year] + t(beta_g)%*%x))
    expect_equal(haz_g, expeta_g/(1+expeta_g))
  }
})

test_that("survival() output", {
  ## External covariates CR
  dim_ext_cr <- 3L;

  seed <- 123
  set.seed(seed)
  beta_d <- rnorm(2 + dim_ext_cr, 0, 1)
  beta_t <- rnorm(2 + dim_ext_cr, 0, 1)
  beta_g <- rnorm(2 + dim_ext_cr, 0, 1)

  beta_cr <- list('d' = beta_d, 't' = beta_t, 'g' = beta_g)
  ## Number of years
  number_years_before <- 6L
  number_years_after <- 2L

  ## year-specific intercepts for `d`, `t` and `g`
  set.seed(seed+1)
  beta0_d <- rnorm(number_years_before, 0, 1)
  beta0_t <- rnorm(number_years_before, 0, 1)
  beta0_g <- rnorm(number_years_after, 0, 1)

  theta_cr <- c(beta_d, beta_t, beta_g, beta0_d, beta0_t, beta0_g)

  set.seed(seed+2)
  x <- rnorm(dim_ext_cr+2, 0, 1)

  # (YEAR_FIRST, YEAR_LAST, YEAR_LAST_EXAM)
  triplets <- list(
    c(2, 5, 6), # regime before
    c(3, 4, 3), # regime after
    c(3, 5, 5), # mixed
    c(3, 5, 4), # mixed
    c(3, 4, 3), # regime after
    c(1, 4, 5) # regime before
  )

  for (trp_idx in 1:length(triplets)) {
    val <- survival(
      YEAR_FIRST = triplets[[trp_idx]][1],
      YEAR_LAST = triplets[[trp_idx]][2],
      THETA_CR = theta_cr,
      COVARIATES = x,
      NYB = number_years_before,
      NYA = number_years_after,
      YEAR_LAST_EXAM = triplets[[trp_idx]][3]
    )

    if(triplets[[trp_idx]][3] > triplets[[trp_idx]][2]){

      Rval <- 1
      for (year in triplets[[trp_idx]][1]:triplets[[trp_idx]][2]) {
        Rval <- Rval * (1 - hazard(
          OUTCOME = 1,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = x,
          NYB = number_years_before,
          NYA = number_years_after
        ) - hazard(
          OUTCOME = 2,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = x,
          NYB = number_years_before,
          NYA = number_years_after
        ))
      }
      expect_equal(val, Rval)

    }else if(triplets[[trp_idx]][3] <= triplets[[trp_idx]][1]){

      Rval <- 1
      for (year in triplets[[trp_idx]][1]:triplets[[trp_idx]][2]) {
        Rval <- Rval * (1 - hazard(
          OUTCOME = 3,
          YEAR = year - triplets[[trp_idx]][3] + 1,
          THETA_CR = theta_cr,
          COVARIATES = x,
          NYB = number_years_before,
          NYA = number_years_after
        ))
      }
      expect_equal(val, Rval)
    }else if((triplets[[trp_idx]][3] > triplets[[trp_idx]][1]) & (triplets[[trp_idx]][3] <= triplets[[trp_idx]][2])){

      Rval <- 1
      for (year in triplets[[trp_idx]][1]:(triplets[[trp_idx]][3]-1)) {
        Rval <- Rval * (1 - hazard(
          OUTCOME = 1,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = x,
          NYB = number_years_before,
          NYA = number_years_after
        ) - hazard(
          OUTCOME = 2,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = x,
          NYB = number_years_before,
          NYA = number_years_after
        ))
      }
      for (year in triplets[[trp_idx]][3]:triplets[[trp_idx]][2]) {
        Rval <- Rval * (1 - hazard(
          OUTCOME = 3,
          YEAR = year - triplets[[trp_idx]][3] + 1,
          THETA_CR = theta_cr,
          COVARIATES = x,
          NYB = number_years_before,
          NYA = number_years_after
        ))
      }
      expect_equal(val, Rval)

    }
  }

})
