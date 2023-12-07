seed <- 1234
{

  ## Setup dimensions
  external_covariates <- 3L     # number of external covariates
  years_before <- 6L            # number of possible years in the regime without graduation
  years_after <- 2L             # number of possible years in the regime with graduation
  grades <- 3L                  # number of grades
  exams <- 4L                  # number of exams


  ## Construct the list
  params_list <- list()

  {
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
  }

  theta <- paramsList2vec(PARAMS_LIST = params_list,
                          DIM_EXT = external_covariates,
                          NYB = years_before,
                          NYA = years_after,
                          N_GRADES = grades,
                          N_EXAMS = exams)

  dim_cr <- 3*(external_covariates+2) + 2*(years_before) + years_after
  theta_cr <- theta[1:dim_cr]
  covariates <- rnorm(2+external_covariates)
}


test_that("hazard() output", {

  for (year in 1:years_before) {
    haz_d <- hazard(
      OUTCOME = 1,
      YEAR = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after
    )

    haz_t <- hazard(
      OUTCOME = 2,
      YEAR = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after
    )

    parsBlock <- params_list[['CR']][['Dropout']]
    expeta_d <- as.numeric(
      exp(
        parsBlock[['Intercepts_year']][year] +
          t(c(parsBlock[['Slope_ability']], parsBlock[['Slope_speed']], parsBlock[['Slope_covariates']]))%*%covariates))
    parsBlock <- params_list[['CR']][['Transfer']]
    expeta_t <- as.numeric(
      exp(
        parsBlock[['Intercepts_year']][year] +
          t(c(parsBlock[['Slope_ability']], parsBlock[['Slope_speed']], parsBlock[['Slope_covariates']]))%*%covariates))
    expect_equal(haz_d, expeta_d/(1+expeta_d+expeta_t))
    expect_equal(haz_t, expeta_t/(1+expeta_d+expeta_t))
  }

  for (year in 1:years_after) {
    haz_g <- hazard(
      OUTCOME = 3,
      YEAR = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after
    )

    parsBlock <- params_list[['CR']][['Graduation']]
    expeta_g <- as.numeric(
      exp(
        parsBlock[['Intercepts_year']][year] +
          t(c(parsBlock[['Slope_ability']], parsBlock[['Slope_speed']], parsBlock[['Slope_covariates']]))%*%covariates))
    expect_equal(haz_g, expeta_g/(1+expeta_g))
  }
})

test_that("survival() output", {

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
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after,
      YEAR_LAST_EXAM = triplets[[trp_idx]][3]
    )

    if(triplets[[trp_idx]][3] > triplets[[trp_idx]][2]){

      Rval <- 1
      for (year in triplets[[trp_idx]][1]:triplets[[trp_idx]][2]) {
        Rval <- Rval * (1 - hazard(
          OUTCOME = 1,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = covariates,
          NYB = years_before,
          NYA = years_after
        ) - hazard(
          OUTCOME = 2,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = covariates,
          NYB = years_before,
          NYA = years_after
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
          COVARIATES = covariates,
          NYB = years_before,
          NYA = years_after
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
          COVARIATES = covariates,
          NYB = years_before,
          NYA = years_after
        ) - hazard(
          OUTCOME = 2,
          YEAR = year,
          THETA_CR = theta_cr,
          COVARIATES = covariates,
          NYB = years_before,
          NYA = years_after
        ))
      }
      for (year in triplets[[trp_idx]][3]:triplets[[trp_idx]][2]) {
        Rval <- Rval * (1 - hazard(
          OUTCOME = 3,
          YEAR = year - triplets[[trp_idx]][3] + 1,
          THETA_CR = theta_cr,
          COVARIATES = covariates,
          NYB = years_before,
          NYA = years_after
        ))
      }
      expect_equal(val, Rval)

    }
  }

})

test_that("check outcomeLik() probability space",{


  probsMarg <- matrix(NA, years_before, 3)
  probs <- matrix(NA, years_before, 3)
  for (year in 1:years_before) {
    for (outcome in 0:2) {
      probs[year,outcome+1] <- outcomeLik(
        OUTCOME = outcome,
        YEAR_FIRST = 1,
        YEAR_LAST = year,
        THETA_CR = theta_cr,
        COVARIATES = covariates,
        NYB = years_before,
        NYA = years_after
      )
    }

    probsMarg[year, 1] <- survival(
      YEAR_FIRST = year,
      YEAR_LAST = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after,
    )
    for (outcome in 1:2) {
      probsMarg[year, outcome+1] <- hazard(
        OUTCOME = outcome,
        YEAR = year,
        THETA_CR = theta_cr,
        COVARIATES = covariates,
        NYB = years_before,
        NYA = years_after
      )
    }
  }
  expect_equal(rowSums(probsMarg), rep(1, years_before))
  expect_equal(sum(probs[,-1])+probs[years_before,1], 1)



  yle <- 3
  probsMarg <- matrix(NA, years_after, 2)
  probs <- matrix(NA, years_after, 2)
  for (year in yle:(yle+years_after-1)) {
    probs[year-yle+1,1] <- outcomeLik(
      OUTCOME = 0,
      YEAR_FIRST = yle,
      YEAR_LAST = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after,
      YEAR_LAST_EXAM = yle
    )
    probs[year-yle+1,2] <- outcomeLik(
      OUTCOME = 3,
      YEAR_FIRST = yle,
      YEAR_LAST = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after,
      YEAR_LAST_EXAM = yle
    )

    probsMarg[year-yle+1, 1] <- survival(
      YEAR_FIRST = year,
      YEAR_LAST = year,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after,
      YEAR_LAST_EXAM = yle
    )

    probsMarg[year-yle+1, 2] <- hazard(
      OUTCOME = 3,
      YEAR = year-yle+1,
      THETA_CR = theta_cr,
      COVARIATES = covariates,
      NYB = years_before,
      NYA = years_after
    )


  }
  rowSums(probsMarg)
  expect_equal(rowSums(probsMarg), rep(1, years_after))
  expect_equal(sum(probs[,-1])+probs[years_after,1], 1)


})
