test_that("pGreaterGrades() output", {

  ## External covariates CR
  dim_ext_cr <- 1L;

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

        cf <- alpha[exam]
        int <- beta_exams[[exam]][grade]

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

  ## External covariates CR
  dim_ext_cr <- 1L;

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

  set.seed(seed+4)
  abilities <- sort(rnorm(3, 0, 2), decreasing = T)
  for (abi_index in 1:length(abilities)) {
    for (grade in 0:grades) {
      for (exam in 1:exams) {

        if(grade == 0){
          Rval <- 1 - exp(beta_exams[[exam]][1]+alpha[exam]*abilities[abi_index])/(1+exp(beta_exams[[exam]][1]+alpha[exam]*abilities[abi_index]))
        }else if(grade < grades) {
          Rval <- exp(beta_exams[[exam]][grade]+alpha[exam]*abilities[abi_index])/(1+exp(beta_exams[[exam]][grade]+alpha[exam]*abilities[abi_index]))-
            exp(beta_exams[[exam]][grade+1]+alpha[exam]*abilities[abi_index])/(1+exp(beta_exams[[exam]][grade+1]+alpha[exam]*abilities[abi_index]))
        }else if(grade == grades){
          Rval <- exp(beta_exams[[exam]][grade]+alpha[exam]*abilities[abi_index])/(1+exp(beta_exams[[exam]][grade]+alpha[exam]*abilities[abi_index]))
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

  ## External covariates CR
  dim_ext_cr <- 1L;

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

  set.seed(seed+6)
  speeds <- sort(rnorm(3,0,5), decreasing=T)

  set.seed(seed+7)
  for (speed_index in 1:length(speeds)) {
    for (day in runif(10, 100, 1000)) {
      for (exam in 1:exams) {

        Rval <- dlnorm(day, gamma[exam]-speeds[speed_index], 1/lambda[exam])
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

        Rval <- plnorm(day, gamma[exam]-speeds[speed_index], 1/lambda[exam])
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

  ## External covariates CR
  dim_ext_cr <- 1L;

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
            pT <- dlnorm(day, gamma[exam]-speed, 1/lambda[exam])
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
            pT <- plnorm(day, gamma[exam]-speed, 1/lambda[exam])
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
