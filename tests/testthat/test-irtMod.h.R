test_that("pGreaterGrades() output", {

  ## External covariates CR
  dim_ext_cr <- 3L;

  seed <- 3
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
  exams <- 3L

  set.seed(seed+3)
  alpha <- rnorm(exams, 0, 1)
  beta_exams <- list()
  for(ex in 1:exams){
    beta_exams[[ex]] <- rnorm(grades, 0, 1)
  }

  theta_irt <- c(alpha, unlist(beta_exams))

  set.seed(seed+4)
  abi <- rnorm(1,0,1)
  for (grade in 1:grades) {
    for (exam in 1:exams) {
      prob <- pGreaterGrades(
        GRADE = grade,
        EXAM = exam,
        THETA_IRT = theta_irt,
        N_GRADES = grades,
        N_EXAMS = exams,
        ABILITY = abi
      )

      cf <- alpha[exam]
      int <- beta_exams[[exam]][grade]

      expect_equal(
        prob, exp(int+cf*abi)/(1+exp(int+cf*abi))
      )
    }
  }

})

test_that("pGrade() output", {

  ## External covariates CR
  dim_ext_cr <- 3L;

  seed <- 3
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
  exams <- 3L

  set.seed(seed+3)
  alpha <- rnorm(exams, 0, 1)
  beta_exams <- list()
  for(ex in 1:exams){
    beta_exams[[ex]] <- rnorm(grades, 0, 1)
  }

  theta_irt <- c(alpha, unlist(beta_exams))

  set.seed(seed+4)
  abi <- rnorm(1,0,1)
  for (grade in 0:grades) {
    for (exam in 1:exams) {

      if(grade == 0){
        Rval <- 1 - exp(beta_exams[[exam]][1]+alpha[exam]*abi)/(1+exp(beta_exams[[exam]][1]+alpha[exam]*abi))
      }else if(grade < grades) {
          Rval <- exp(beta_exams[[exam]][grade]+alpha[exam]*abi)/(1+exp(beta_exams[[exam]][grade]+alpha[exam]*abi))-
            exp(beta_exams[[exam]][grade+1]+alpha[exam]*abi)/(1+exp(beta_exams[[exam]][grade+1]+alpha[exam]*abi))
      }else if(grade == grades){
        Rval <- exp(beta_exams[[exam]][grade]+alpha[exam]*abi)/(1+exp(beta_exams[[exam]][grade]+alpha[exam]*abi))
      }
      val <- pGrade(
        GRADE = grade,
        EXAM = exam,
        THETA_IRT = theta_irt,
        N_GRADES = grades,
        N_EXAMS = exams,
        ABILITY = abi
      )
      expect_equal(Rval, val)
    }
  }

})

test_that("pTimeExam() output", {

  ## External covariates CR
  dim_ext_cr <- 3L;

  seed <- 3
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
  exams <- 3L

  set.seed(seed+3)
  alpha <- rnorm(exams, 0, 1)
  beta_exams <- list()
  for(ex in 1:exams){
    beta_exams[[ex]] <- rnorm(grades, 0, 1)
  }

  set.seed(seed+4)
  gamma <- rnorm(exams, 0, 1)

  set.seed(seed+5)
  lambda <- runif(exams, .1, 3)

  theta_irt <- c(alpha, unlist(beta_exams), gamma, lambda)

  set.seed(seed+6)
  speed <- rnorm(1,0,1)

  set.seed(seed+7)
  for (day in runif(10, 100, 1000)) {
    for (exam in 1:exams) {
      Rval <- plnorm(day, gamma[exam]-speed, 1/lambda[exam])
      val <- pTimeExam(
        EXAM = exam,
        DAY = day,
        THETA_IRT = theta_irt,
        N_GRADES = grades,
        N_EXAMS = exams,
        SPEED = speed
      )
      expect_equal(val, Rval)
    }
  }


})
