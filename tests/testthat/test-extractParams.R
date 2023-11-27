test_that("extract_params() output", {
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

  set.seed(4)
  gamma <- rnorm(exams, 0, 1)

  set.seed(5)
  lambda <- rnorm(exams, 0, 1)

  theta_irt <- c(alpha, unlist(beta_exams), gamma, lambda)

  theta <- c(theta_cr, theta_irt)

  ## Test CR parameters extraction
  for (outcome in c('d', 't', 'g')) {
    coeff <- extract_params(
      THETA = theta,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      PAR_TYPE = 1,
      OUTCOME = outcome,
      EXAM = NA
    )
    expect_identical(coeff, beta_cr[[outcome]])

    inter <- extract_params(
      THETA = theta,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      PAR_TYPE = 2,
      OUTCOME = outcome,
      EXAM = NA
    )
    expect_identical(inter, beta0_cr[[outcome]])

  }

  ## Test IRT parameters extraction
  for (exam in 1:exams) {
    coeff <- extract_params(
      THETA = theta,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      PAR_TYPE = 1,
      OUTCOME = NA,
      EXAM = exam
    )
    expect_identical(coeff, alpha[exam])

    inter <- extract_params(
      THETA = theta,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      PAR_TYPE = 2,
      OUTCOME = NA,
      EXAM = exam
    )
    expect_identical(inter,  beta_exams[[exam]])

    speed <- extract_params(
      THETA = theta,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      PAR_TYPE = 3,
      OUTCOME = NA,
      EXAM = exam
    )
    expect_identical(speed,  c(gamma[exam], lambda[exam]))

  }


})
