test_that("pGreaterGrades() output", {

  seed <- 123

  ## Setup dimensions
  external_covariates <- 3L     # number of external covariates
  years_before <- 4L            # number of possible years in the regime without graduation
  years_after <- 2L             # number of possible years in the regime with graduation
  grades <- 3L                  # number of grades
  exams <- 5L                  # number of exams


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

  ## Parameter-blocks related to  <-  IRT-modelling
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

  theta <- paramsList2vec(PARAMS_LIST = params_list,
                          DIM_EXT = external_covariates,
                          NYB = years_before,
                          NYA = years_after,
                          N_GRADES = grades,
                          N_EXAMS = exams)

  theta_irt <- theta[
    (3*(external_covariates+2) + 2*(years_before) + years_after+1):
      (3*(external_covariates+2) + 2*(years_before) + years_after + 3*exams + exams*grades + 2)]

  set.seed(seed+4)
  abilities <- sort(rnorm(3, 0, 1), decreasing = T)
  for (abi_index in 1:length(abilities)) {
    for (grade in 1:grades) {
      for (exam in 1:exams) {
        prob <- pGreaterGrades(
          GRADE = grade,
          EXAM = exam,
          THETA_IRT = theta_irt,
          N_GRADES = grades,
          N_EXAMS = exams,
          ABILITY = abilities[abi_index]
        )

        cf <- params_list[['IRT']][['Exams_slopes']][exam]
        int <- params_list[['IRT']][['Exams_grades_intercepts']][exam, grade]

        # check probs matches
        expect_equal(
          prob, exp(int+cf*abilities[abi_index])/(1+exp(int+cf*abilities[abi_index]))
        )

        if(grade>1){
          probprev <- pGreaterGrades(
            GRADE = grade-1,
            EXAM = exam,
            THETA_IRT = theta_irt,
            N_GRADES = grades,
            N_EXAMS = exams,
            ABILITY = abilities[abi_index]
          )
          # check probs are decresing with higher grades
          expect_true(probprev>prob)
        }

        if(abi_index>1){
          probprev <- pGreaterGrades(
            GRADE = grade,
            EXAM = exam,
            THETA_IRT = theta_irt,
            N_GRADES = grades,
            N_EXAMS = exams,
            ABILITY = abilities[abi_index-1]
          )
          # check probs are decresing with lower ability
          expect_true(probprev>prob)
        }
      }
    }
  }


})

test_that("pGrade() output", {

  seed <- 123

  ## Setup dimensions
  external_covariates <- 3L     # number of external covariates
  years_before <- 4L            # number of possible years in the regime without graduation
  years_after <- 2L             # number of possible years in the regime with graduation
  grades <- 3L                  # number of grades
  exams <- 5L                  # number of exams


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

  theta <- paramsList2vec(PARAMS_LIST = params_list,
                          DIM_EXT = external_covariates,
                          NYB = years_before,
                          NYA = years_after,
                          N_GRADES = grades,
                          N_EXAMS = exams)

  theta_irt <- theta[
    (3*(external_covariates+2) + 2*(years_before) + years_after+1):
      (3*(external_covariates+2) + 2*(years_before) + years_after + 3*exams + exams*grades + 2)]

  set.seed(seed+4)
  abilities <- sort(rnorm(3, 0, 2), decreasing = T)
  for (abi_index in 1:length(abilities)) {
    for (grade in 0:grades) {
      for (exam in 1:exams) {

        if(grade == 0){
          Rval <- 1 - exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,1]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index])/
            (1+exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,1]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index]))
        }else if(grade < grades) {
          Rval <- exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,grade]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index])/
            (1+exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,grade]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index]))-
            exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,grade+1]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index])/
            (1+exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,grade+1]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index]))
        }else if(grade == grades){
          Rval <- exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,grade]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index])/
            (1+exp(params_list[['IRT']][['Exams_grades_intercepts']][exam,grade]+params_list[['IRT']][['Exams_slopes']][exam]*abilities[abi_index]))
        }
        val <- pGrade(
          GRADE = grade,
          EXAM = exam,
          THETA_IRT = theta_irt,
          N_GRADES = grades,
          N_EXAMS = exams,
          ABILITY = abilities[abi_index]
        )

        # check umeric output
        expect_equal(Rval, val)

      }
    }

  }

})

test_that("pTimeExam() output", {

  seed <- 123

  ## Setup dimensions
  external_covariates <- 3L     # number of external covariates
  years_before <- 4L            # number of possible years in the regime without graduation
  years_after <- 2L             # number of possible years in the regime with graduation
  grades <- 3L                  # number of grades
  exams <- 5L                  # number of exams


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

  theta <- paramsList2vec(PARAMS_LIST = params_list,
                          DIM_EXT = external_covariates,
                          NYB = years_before,
                          NYA = years_after,
                          N_GRADES = grades,
                          N_EXAMS = exams)

  theta_irt <- theta[
    (3*(external_covariates+2) + 2*(years_before) + years_after+1):
      (3*(external_covariates+2) + 2*(years_before) + years_after + 3*exams + exams*grades + 2)]

  set.seed(seed+6)
  speeds <- sort(rnorm(3,0,5), decreasing=T)

  set.seed(seed+7)
  for (speed_index in 1:length(speeds)) {
    for (day in runif(10, 100, 1000)) {
      for (exam in 1:exams) {

        Rval <- dlnorm(day,
                       params_list[['IRT']][['Exams_average_time']][exam]-speeds[speed_index],
                       1/params_list[['IRT']][['Exams_variability_time']][exam])
        val <- pTimeExam(
          EXAM = exam,
          DAY = day,
          THETA_IRT = theta_irt,
          N_GRADES = grades,
          N_EXAMS = exams,
          SPEED = speeds[speed_index],
          CDFFLAG = F
        )
        expect_equal(val, Rval)

        Rval <- plnorm(day,
                       params_list[['IRT']][['Exams_average_time']][exam]-speeds[speed_index],
                       1/params_list[['IRT']][['Exams_variability_time']][exam])
        val <- pTimeExam(
          EXAM = exam,
          DAY = day,
          THETA_IRT = theta_irt,
          N_GRADES = grades,
          N_EXAMS = exams,
          SPEED = speeds[speed_index],
          CDFFLAG = T
        )
        #check numeric output
        expect_equal(val, Rval)

        if(speed_index>1){
          val <- pTimeExam(
            EXAM = exam,
            DAY = day,
            THETA_IRT = theta_irt,
            N_GRADES = grades,
            N_EXAMS = exams,
            SPEED = speeds[speed_index],
            CDFFLAG = T
          )
          valprev <- pTimeExam(
            EXAM = exam,
            DAY = day,
            THETA_IRT = theta_irt,
            N_GRADES = grades,
            N_EXAMS = exams,
            SPEED = speeds[speed_index-1],
            CDFFLAG = T
          )
          #check that p(T<t) decreases with decreasing speeds
          expect_true(valprev>val)

        }


      }
    }

  }


})

test_that("examLik() output", {

  seed <- 123

  ## Setup dimensions
  external_covariates <- 3L     # number of external covariates
  years_before <- 4L            # number of possible years in the regime without graduation
  years_after <- 2L             # number of possible years in the regime with graduation
  grades <- 3L                  # number of grades
  exams <- 5L                  # number of exams


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

  theta <- paramsList2vec(PARAMS_LIST = params_list,
                          DIM_EXT = external_covariates,
                          NYB = years_before,
                          NYA = years_after,
                          N_GRADES = grades,
                          N_EXAMS = exams)

  theta_irt <- theta[
    (3*(external_covariates+2) + 2*(years_before) + years_after+1):
      (3*(external_covariates+2) + 2*(years_before) + years_after + 3*exams + exams*grades + 2)]

  # set.seed(seed+6)
  # speed <- rnorm(1,0,1)
  #
  # set.seed(seed+7)
  # ability <- rnorm(1,0,1)

  set.seed(seed+6)
  for (ability in rnorm(2,0,1)) {
    for (speed in rnorm(2,0,1)) {
      for (day in runif(2, 100, 1000)) {
        for (exam in 1:exams) {
          for (grade in 1:grades) {
            pG <- pGrade(GRADE = grade, EXAM = exam, THETA_IRT = theta_irt, N_GRADES = grades, N_EXAMS = exams, ABILITY = ability)
            pT <- dlnorm(day,
                         params_list[['IRT']][['Exams_average_time']][exam]-speed,
                         1/params_list[['IRT']][['Exams_variability_time']][exam])
            Rval <- pT*pG
            val <- examLik(
              EXAM = exam,
              GRADE = grade,
              DAY = day,
              OBSFLAG = T,
              THETA_IRT = theta_irt,
              N_GRADES = grades,
              N_EXAMS = exams,
              ABILITY = ability,
              SPEED = speed
            )
            expect_equal(val, Rval)

            pG <- pGreaterGrades(GRADE = 1, EXAM = exam, THETA_IRT = theta_irt, N_GRADES = grades, N_EXAMS = exams, ABILITY = ability)
            pT <- plnorm(day,
                         params_list[['IRT']][['Exams_average_time']][exam]-speed,
                         1/params_list[['IRT']][['Exams_variability_time']][exam])
            Rval <- 1-pT*pG
            val <- examLik(
              EXAM = exam,
              GRADE = grade,
              DAY = day,
              OBSFLAG = F,
              THETA_IRT = theta_irt,
              N_GRADES = grades,
              N_EXAMS = exams,
              ABILITY = ability,
              SPEED = speed
            )
            expect_equal(val, Rval)
          }
        }

      }
    }
  }



})
