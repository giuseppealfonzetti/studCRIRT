test_that("reparInt() errors",{

  set.seed(123)
  param <- rnorm(5)
  thr <- reparInt(param, CON2UN = F)

  expect_error(reparInt(sample(thr,length(thr)), CON2UN = T))

})


test_that("reparInt() output",{

  set.seed(123)
  for (trial in 1:10) {
    param <- rnorm(5)
    thr <- reparInt(param, CON2UN = F)
    expect_true(sum(sort(thr, decreasing = T) == thr)==length(thr))

    expect_equal(reparInt(thr, CON2UN = T), param)

  }
})


test_that("paramsList2vec() and paramsVec2list() errors", {

  seed <- 123

  ## Setup dimensions
  external_covariates <- 3L     # number of external covariates
  years_before <- 6L            # number of possible years in the regime without graduation
  years_after <- 2L             # number of possible years in the regime with graduation
  grades <- 4L                  # number of grades
  exams <- 10L                  # number of exams


  ## Construct the list
  params_list <- list()

  # Three main blocks
  params_list[['CR']] <- list()
  params_list[['IRT']] <- list()
  params_list[['LAT']] <- list()

  # Parameter-blocks related to competing risks
  set.seed(seed)
  params_list[['CR']][['Dropout']] <- list(
   'Slope_ability' = rnorm(1),
   'Slope_speed' = rnorm(1),
   'Slope_covariates' = rnorm(external_covariates),
   'Intercepts_year' = rnorm(years_before)
  )

  set.seed(seed+1)
  params_list[['CR']][['Transfer']] <- list(
   'Slope_ability' = rnorm(1),
   'Slope_speed' = rnorm(1),
   'Slope_covariates' = rnorm(external_covariates),
   'Intercepts_year' = rnorm(years_before)
  )

  set.seed(seed+2)
  params_list[['CR']][['Graduation']] <- list(
   'Slope_ability' = rnorm(1),
   'Slope_speed' = rnorm(1),
   'Slope_covariates' = rnorm(external_covariates),
   'Intercepts_year' = rnorm(years_after)
  )

  ## Parameter-blocks related to exams IRT-modelling
  set.seed(seed+3)
  params_list[['IRT']][['Exams_slopes']] <- runif(exams, 0, 5)

  set.seed(seed+4)
  examsInt <- matrix(0, exams, grades)
  for(ex in 1:exams){
  examsInt[ex,] <- sort(rnorm(grades, 0, 5), decreasing = TRUE)
  }
  params_list[['IRT']][['Exams_grades_intercepts']] <- examsInt

  set.seed(seed+5)
  params_list[['IRT']][['Exams_average_time']] <- rnorm(exams)

  set.seed(seed+6)
  params_list[['IRT']][['Exams_variability_time']] <- runif(exams)

  params_list[['LAT']][['Corr']] <- .6
  params_list[['LAT']][['Speed_variability']] <- 1.5


  ## TESTS

  tmp <- params_list
  for (tst in c('a', 3.5, NA, NaN, TRUE)) {
    tmp[['CR']][['Dropout']][['Slope_ability']] <- tst
    expect_error(paramsList2vec(PARAMS_LIST = tmp,
                                DIM_EXT = external_covariates,
                                NYB = years_before,
                                NYA = years_after,
                                N_GRADES = grades,
                                N_EXAMS = exams))
  }

  tmp <- params_list
  for (tst in c('a', 3.5, NA, NaN, TRUE)) {
    expect_error(paramsList2vec(PARAMS_LIST = tmp,
                                DIM_EXT = external_covariates,
                                NYB = tst,
                                NYA = years_after,
                                N_GRADES = grades,
                                N_EXAMS = exams))
  }

  for (tst in c('Exams_slopes', 'Exams_average_time', 'Exams_variability_time', 'Exams_grades_intercepts')) {
    tmp <- params_list
    tmp[['IRT']][[tst]] <- rnorm(20)
    expect_error(paramsList2vec(PARAMS_LIST = tmp,
                                DIM_EXT = external_covariates,
                                NYB = years_before,
                                NYA = years_after,
                                N_GRADES = grades,
                                N_EXAMS = exams))
  }


  #### CHECK CONSTRAINTS ####

  #### Latent correlation ####
  tmp <- params_list
  tmp[['LAT']][['Corr']] <- 1
  expect_error(paramsList2vec(PARAMS_LIST = tmp,
                              DIM_EXT = external_covariates,
                              NYB = years_before,
                              NYA = years_after,
                              N_GRADES = grades,
                              N_EXAMS = exams))

  #### latent var #####
  tmp <- params_list
  tmp[['LAT']][['Speed_variability']] <- -1
  expect_error(paramsList2vec(PARAMS_LIST = tmp,
                              DIM_EXT = external_covariates,
                              NYB = years_before,
                              NYA = years_after,
                              N_GRADES = grades,
                              N_EXAMS = exams))

  #### grades intercepts ####
  tmp <- params_list
  mod <- tmp[['IRT']][['Exams_grades_intercepts']]
  mod[1,] <- sample(mod[1,],grades)
  tmp[['IRT']][['Exams_grades_intercepts']] <- mod
  expect_error(paramsList2vec(PARAMS_LIST = tmp,
                              DIM_EXT = external_covariates,
                              NYB = years_before,
                              NYA = years_after,
                              N_GRADES = grades,
                              N_EXAMS = exams))

  #### exams var ####
  tmp <- params_list
  mod <- tmp[['IRT']][['Exams_variability_time']]
  mod[1] <- -1
  tmp[['IRT']][['Exams_variability_time']] <- mod
  expect_error(paramsList2vec(PARAMS_LIST = tmp,
                              DIM_EXT = external_covariates,
                              NYB = years_before,
                              NYA = years_after,
                              N_GRADES = grades,
                              N_EXAMS = exams))



})
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
  grades <- 4L
  exams <- 5L

  set.seed(seed+3)
  alpha <- runif(exams, 0, 5)
  beta_exams <- matrix(0, exams, grades)
  for(ex in 1:exams){
    beta_exams[ex,] <- sort(rnorm(grades, 0, 5), decreasing = T)
  }

  set.seed(seed+4)
  gamma <- runif(exams, 0, 5)

  set.seed(seed+5)
  lambda <- runif(exams, .1, 3)

  theta_irt <- c(alpha, as.numeric(unlist(lapply(as.list(as.data.frame(t(beta_exams))), reparInt))),
                 gamma, log(lambda))



  ### theta latent distr ####
  theta_lat <- c(atanh(.3), log(1.5))


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
  params_list[['IRT']][['Exams_grades_intercepts']] <- beta_exams
  params_list[['IRT']][['Exams_average_time']] <- gamma
  params_list[['IRT']][['Exams_variability_time']] <- lambda
  params_list[['LAT']][['Corr']] <- tanh(theta_lat[1])
  params_list[['LAT']][['Speed_variability']] <- exp(theta_lat[2])

  expect_equal(
    paramsList2vec(
      PARAMS_LIST = params_list,
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams), theta
  )

  expect_equal(
    paramsList2vec(
      PARAMS_LIST = paramsVec2list(
        THETA = theta,
        DIM_EXT = dim_ext_cr,
        NYB = number_years_before,
        NYA = number_years_after,
        N_GRADES = grades,
        N_EXAMS = exams),
      DIM_EXT = dim_ext_cr,
      NYB = number_years_before,
      NYA = number_years_after,
      N_GRADES = grades,
      N_EXAMS = exams), theta
  )
})
