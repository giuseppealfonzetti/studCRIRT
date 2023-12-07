#' Generate random grades
#'
#' @param THETA_IRT Portion of the parameter vector related to the IRT model.
#' @param N_GRADES Number of grades modelled.
#' @param N_EXAMS Number of exams.
#' @param ABILITY Ability value.
#' @param SEED Random seed.
#'
#' @returns It returns a vector of grades.
#'
#' @export
#' @importFrom stats rmultinom
rngGrades <- function(THETA_IRT, N_GRADES, N_EXAMS, ABILITY, SEED = 123){

  set.seed(123)
  gradesVec <- rep(NA, N_EXAMS)
  probsMat <- matrix(0, N_EXAMS, N_GRADES+1)
  for (exam in 1:N_EXAMS) {
    probs <- c()
    for(grade in 0:N_GRADES){
      probsMat[exam, grade+1] <- pGrade(GRADE = grade,
                             EXAM = exam,
                             THETA_IRT = THETA_IRT,
                             N_GRADES = N_GRADES,
                             N_EXAMS = N_EXAMS,
                             ABILITY = ABILITY)
    }

    gradesVec[exam] <- which(rmultinom(1, 1, probsMat[exam,])==1)-1
  }

  return(list(
    'probs' = probsMat,
    'grades' = gradesVec
  ))
}

#' Generate random Times
#'
#' @param THETA_IRT Portion of the parameter vector related to the IRT model.
#' @param N_GRADES Number of grades modelled.
#' @param N_EXAMS Number of exams.
#' @param SPEED Speed value.
#' @param SEED Random seed.
#'
#' @returns It returns a vector of times

#' @export
#' @importFrom stats rlnorm
rngTimes <- function(THETA_IRT, N_GRADES, N_EXAMS, SPEED, SEED = 123){

  timesVec <- rep(NA, N_EXAMS)
  for (exam in 1:N_EXAMS) {
    av <- extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 3, exam)[1] - SPEED
    va <- 1/extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 4, exam)[1]
    timesVec[exam] <- rlnorm(1, av, va)
  }

  return(timesVec)
}

