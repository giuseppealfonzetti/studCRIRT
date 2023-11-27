
#' @export
extract_params_CR <- function(THETA_CR, DIM_EXT, NYB, NYA){
  beta_d <- THETA_CR[1:(DIM_EXT+2)]
  beta_t <- THETA_CR[(DIM_EXT+3):(2*DIM_EXT+4)]
  beta_g <- THETA_CR[(2*DIM_EXT+5):(3*DIM_EXT+6)]
  beta0_d <- THETA_CR[(3*DIM_EXT+7):(3*DIM_EXT+6+NYB)]
  beta0_t <- THETA_CR[(3*DIM_EXT+7+NYB):(3*DIM_EXT+6+2*NYB)]
  beta0_g <- THETA_CR[(3*DIM_EXT+7+2*NYB):length(THETA_CR)]



  out <- list(
    'beta_d' = beta_d,
    'beta_t' = beta_t,
    'beta_g' = beta_g,
    'beta0_d' = beta0_d,
    'beta0_t' = beta0_t,
    'beta0_g' = beta0_g
  )

  return(out)
}

#' Extract parameters of interest from theta vector
#'
#' @param THETA parameter vector
#' @param DIM_EXT number of external covariates in the competing risk model.
#' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param N_GRADES number of grades modelled.
#' @param N_EXAMS number of exams.
#' @param PAR_TYPE 1 for regression coeff, 2 for intercepts, 3 for speed-exam-related parameters
#' @param OUTCOME outcome of interest. Possible values: `d`, `t` and `g`.
#' @param EXAM exam of interest. Posible values in `1:N_EXAMS`.
#'
#' @returns It returns a vector containing the parameters of interests if either `OUTCOME` or `EXAM` have been specified.
#'  It returns a complete list of the parameters if both `OUTCOME` and `EXAM` are set to `NA`.
#' @export
extract_params <- function(THETA, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS, PAR_TYPE = NA, OUTCOME = NA, EXAM = NA){
  if(!is.numeric(THETA)) stop('THETA not accepted.')
  if(!is.integer(DIM_EXT)) stop('DIM_EXT not accepted.')
  if(!is.integer(NYB)) stop('NYB not accepted.')
  if(!is.integer(NYA)) stop('NYA not accepted.')
  if(!is.integer(N_GRADES)) stop('N_GRADES not accepted.')
  if(!is.integer(N_EXAMS)) stop('N_EXAMS not accepted.')

  if(!(OUTCOME %in% c('d', 't', 'g', NA) )) stop('OUTCOME not accepted.')
  if(!(EXAM %in% c(1:N_EXAMS, NA) )) stop('EXAM chosen not accepted.')
  if(!is.na(EXAM) & !is.na(OUTCOME)) stop('Both EXAM and OUTCOME selected.')
  out <- c()

  theta_cr <- THETA[1:(3*DIM_EXT + 6 + 2*NYB + NYA)]
  theta_irt <-THETA[(3*DIM_EXT + 7 + 2*NYB + NYA):(3*DIM_EXT + 6 + 2*NYB + NYA + N_EXAMS*(N_GRADES+1) + 2*N_EXAMS)]

  if(!is.na(EXAM)){
    if(!(PAR_TYPE %in% c(1,2,3))) stop('PAR_TYPE not accepted.')
    out <- extract_params_irt(
      THETA_IRT = theta_irt,
      N_GRADES = N_GRADES,
      N_EXAMS = N_EXAMS,
      OPTION = PAR_TYPE,
      EXAM = EXAM
    )
  }

  if(!is.na(OUTCOME)){
    if(!(PAR_TYPE %in% c(1,2))) stop('PAR_TYPE not accepted.')
    outcode <- switch(OUTCOME, 'd' = 1, 't' = 2, 'g' = 3)
    option <- (PAR_TYPE-1)*3 + outcode
    out <- extract_params_cr(
      THETA_CR = theta_cr,
      DIM_EXT = DIM_EXT,
      NYB = NYB,
      NYA = NYA,
      OPTION = option
    )

  }

  if(is.na(EXAM) & is.na(OUTCOME)){out <- list('cr' = theta_cr, 'irt' = theta_irt)}


  return(out)
}
