#' Approximate the marginal likelihood of a single observation
#'
#' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
#' @param EXTCOVARIATES  Vector of external covariates.
#' @param EXAMS_GRADES Vector of grades.
#' @param EXAMS_DAYS Vector of times.
#' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
#' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
#' @param YEAR Year of evaluation.
#' @param N_GRADES Number of grades modelled.
#' @param N_EXAMS Number of exams modelled
#' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
#' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
#' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
#' @param LOGFLAG Set TRUE to return log value.
#' @param CUBATURE_METHOD Argument to be passed to [cubature::cubintegrate].
#'
#' @export
marginal_likelihood <- function(
    THETA,
    EXTCOVARIATES,
    EXAMS_GRADES,
    EXAMS_DAYS,
    EXAMS_OBSFLAG,
    OUTCOME,
    YEAR,
    N_GRADES,
    N_EXAMS,
    NYB,
    NYA,
    YEAR_LAST_EXAM = 100,
    LOGFLAG = FALSE,
    CUBATURE_METHOD = 'hcubature'
){

  maxTime <- 365*YEAR
  obs_vec <- EXAMS_GRADES>0&EXAMS_DAYS<maxTime
  Rfun <- function(LAT){
    complete_likelihood(
      THETA = THETA,
      EXTCOVARIATES = EXTCOVARIATES,
      EXAMS_GRADES = EXAMS_GRADES,
      EXAMS_DAYS = EXAMS_DAYS,
      EXAMS_OBSFLAG = EXAMS_OBSFLAG,
      OUTCOME = OUTCOME,
      YEAR = YEAR,
      N_GRADES = N_GRADES,
      N_EXAMS = N_EXAMS,
      NYB = NYB,
      NYA = NYA,
      ABILITY = LAT[1],
      SPEED = LAT[2],
      YEAR_LAST_EXAM = YEAR_LAST_EXAM,
      LOGFLAG = FALSE
    )
  }

  out <- cubature::cubintegrate(f= Rfun, lower = c(-Inf, -Inf), upper = c(Inf, Inf), method = CUBATURE_METHOD)
  if(LOGFLAG){out <- exp(out)}
  return(out)
}


#' Approximate the marginal likelihood of the sample
#'
#' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
#' @param EXTCOVARIATES  Matrix of external covariates with `n` rows.
#' @param EXAMS_GRADES List of `n` vectors of grades.
#' @param EXAMS_DAYS List of `n` vectors of times.
#' @param EXAMS_OBSFLAG VList of `n` vectors of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
#' @param OUTCOME  Vector of `n` outcomes: `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
#' @param YEAR Vector of `n` years of evaluation.
#' @param N_GRADES Number of grades modelled.
#' @param N_EXAMS Number of exams modelled
#' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
#' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
#' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
#' @param LOGFLAG Set TRUE to return log value.
#' @param CUBATURE_METHOD Argument to be passed to [cubature::cubintegrate].
#'
#' @export
sample_marginal_likelihood <- function(
    THETA,
    EXTCOVARIATES,
    EXAMS_GRADES,
    EXAMS_DAYS,
    EXAMS_OBSFLAG,
    OUTCOME,
    YEAR,
    N_GRADES,
    N_EXAMS,
    NYB,
    NYA,
    YEAR_LAST_EXAM = 100,
    LOGFLAG = FALSE,
    CUBATURE_METHOD = 'hcubature'
){

  maxTime <- 365*YEAR
  obs_vec <- EXAMS_GRADES>0&EXAMS_DAYS<maxTime
  Rfun <- function(LAT){
    complete_likelihood(
      THETA = THETA,
      EXTCOVARIATES = EXTCOVARIATES,
      EXAMS_GRADES = EXAMS_GRADES,
      EXAMS_DAYS = EXAMS_DAYS,
      EXAMS_OBSFLAG = EXAMS_OBSFLAG,
      OUTCOME = OUTCOME,
      YEAR = YEAR,
      N_GRADES = N_GRADES,
      N_EXAMS = N_EXAMS,
      NYB = NYB,
      NYA = NYA,
      ABILITY = LAT[1],
      SPEED = LAT[2],
      YEAR_LAST_EXAM = YEAR_LAST_EXAM,
      LOGFLAG = FALSE
    )
  }

  out <- cubature::cubintegrate(f= Rfun, lower = c(-Inf, -Inf), upper = c(Inf, Inf), method = CUBATURE_METHOD)
  if(LOGFLAG){out <- exp(out)}
  return(out)
}
