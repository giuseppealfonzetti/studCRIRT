seed <- 123

external_covariates <- 3L     # number of external covariates
years_before <- 6L            # number of possible years without graduation risk
years_after <- 2L             # number of possible years  with graduation risk
grades <- 4L                  # number of grades
exams <- 50L                  # number of exams

dim_cr <- 3*(external_covariates+2) + 2*(years_before) + years_after
dim_irt <- 3*exams + exams*grades

set.seed(seed)

test_that('rngGrades()', {
  for (i in 1:10) {
    set.seed(seed+i)
    theta <- rnorm(dim_cr+dim_irt+2)
    lat <- rnorm(2)
    grades_vec <- rngGrades(
      THETA_IRT = theta[(dim_cr+1):(dim_cr+dim_irt)],
      N_GRADES = grades,
      N_EXAMS = exams,
      ABILITY = lat[1],
      SEED = seed+i
    )$grades
    expect_true(sum(grades_vec>=0&grades_vec<=grades)==exams)
  }
  for (grade in 0:grades) {
    expect_true(grade %in% grades_vec)
  }
})

test_that('rngTimes()', {
  for (i in 1:10) {
    set.seed(seed+i)
    theta <- rnorm(dim_cr+dim_irt+2)
    lat <- rnorm(2)
    times_vec <- rngTimes(
      THETA_IRT = theta[(dim_cr+1):(dim_cr+dim_irt)],
      N_GRADES = grades,
      N_EXAMS = exams,
      SPEED = lat[2],
      SEED = seed
    )
    expect_true(sum(times_vec>=0)==exams)
  }
})
