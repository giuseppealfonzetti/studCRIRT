test_that("paramsList2vec() and paramsVec2list() output", {

  ## External covariates CR
  dim_ext_cr <- 3L;

  seed <- 1234
  set.seed(seed+1)
  beta_d <- rnorm(2 + dim_ext_cr, 0, 1)
  beta_t <- rnorm(2 + dim_ext_cr, 0, 1)
  beta_g <- rnorm(2 + dim_ext_cr, 0, 1)

  beta_cr <- list('d' = beta_d, 't' = beta_t, 'g' = beta_g)
  ## Number of years
  number_years_before <- 6L
  number_years_after <- 2L

  ## Year-related intercept
  set.seed(seed+2)

  ## year-specific intercepts for `d`, `t` and `g`
  beta0_d <- rnorm(number_years_before, 0, 1)
  beta0_t <- rnorm(number_years_before, 0, 1)
  beta0_g <- rnorm(number_years_after, 0, 1)
  beta0_cr <- list('d' = beta0_d, 't' = beta0_t, 'g' = beta0_g)

  theta_cr <- c(beta_d, beta_t, beta_g, beta0_d, beta0_t, beta0_g)

  ## Exam-related params
  grades <- 3L
  exams <- 3L

  set.seed(seed+3)
  alpha <- runif(exams, 0, 5)
  beta_exams <- list()
  for(ex in 1:exams){
    beta_exams[[ex]] <- sort(rnorm(grades, 0, 5), decreasing = T)
  }

  set.seed(seed+4)
  gamma <- runif(exams, 0, 5)

  set.seed(seed+5)
  lambda <- runif(exams, .1, 3)

  theta_irt <- c(alpha, unlist(beta_exams), gamma, lambda)


  ### theta latent distr ####
  theta_lat <- c(.3, 1.5)


  theta <- c(theta_cr, theta_irt, theta_lat)

  # construct list
  params_list <- list()
  params_list[['CR']] <- list()
  params_list[['IRT']] <- list()
  params_list[['LAT']] <- list()
  params_list[['CR']][['Dropout']] <- list(
    'Slope_ability' = beta_d[1],
    'Slope_speed' = beta_d[2],
    'Slope_covariates' = beta_d[3:(2+dim_ext_cr)],
    'Intercepts_year' = beta0_d
  )
  params_list[['CR']][['Transfer']] <- list(
    'Slope_ability' = beta_t[1],
    'Slope_speed' = beta_t[2],
    'Slope_covariates' = beta_t[3:(2+dim_ext_cr)],
    'Intercepts_year' = beta0_t
  )
  params_list[['CR']][['Graduation']] <- list(
    'Slope_ability' = beta_g[1],
    'Slope_speed' = beta_g[2],
    'Slope_covariates' = beta_g[3:(2+dim_ext_cr)],
    'Intercepts_year' = beta0_g
  )
  params_list[['IRT']][['Exams_slopes']] <- alpha
  params_list[['IRT']][['Exams_grades_intercepts']] <- unlist(beta_exams)
  params_list[['IRT']][['Exams_average_time']] <- gamma
  params_list[['IRT']][['Exams_variability_time']] <- lambda
  params_list[['LAT']][['Corr']] <- theta_lat[1]
  params_list[['LAT']][['Speed_variability']] <- theta_lat[2]

  expect_identical(
    paramsList2vec(
      PARAMS_LIST = params_list,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      CONSTRFLAG = T), theta
  )

  expect_identical(
    paramsList2vec(
      PARAMS_LIST = paramsVec2list(
        THETA = theta,
        DIM_EXT = dim_ext_cr,
        NYB = number_years_before,
        NYA = number_years_after,
        N_GRADES = grades,
        N_EXAMS = exams,
        CONSTRFLAG = T),
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams,
      CONSTRFLAG = T), theta
  )
})
