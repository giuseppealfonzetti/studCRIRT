#' Construct parameter vector from list
#'
#' @param PARAMS_LIST Parameter list. Structure to be documented.
#' @param DIM_EXT number of external covariates in the competing risk model.
#' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param N_GRADES number of grades modelled.
#' @param N_EXAMS number of exams.
#' @param CONSTRFLAG TRUE if input is defined in the constrained parameter space.
#'
#' @export
paramsList2vec <- function(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS, CONSTRFLAG = T){

  # Check if all parameter blocks are specified
  sapply(c('CR', 'IRT', 'LAT'), function(x) {if(!(x %in% names(PARAMS_LIST))) stop(paste0('Parameters block ', x, ' not detected'))})
  sapply(c('Dropout', 'Transfer', 'Graduation'), function(x) {if(!(x %in% names(PARAMS_LIST[['CR']]))) stop(paste0('Parameters block ', x, ' not detected'))})
  sapply(c('Exams_slopes', 'Exams_grades_intercepts', 'Exams_average_time', 'Exams_variability_time'), function(x) {if(!(x %in% names(PARAMS_LIST[['IRT']]))) stop(paste0('Parameters block ', x, ' not detected'))})
  sapply(c('Corr', 'Speed_variability'), function(x) {if(!(x %in% names(PARAMS_LIST[['LAT']]))) stop(paste0('Parameters block ', x, ' not detected'))})

  # Check input types
  if(!is.integer(DIM_EXT)) stop('DIM_EXT not accepted.')
  if(!is.integer(NYB)) stop('NYB not accepted.')
  if(!is.integer(NYA)) stop('NYA not accepted.')
  if(!is.integer(N_GRADES)) stop('N_GRADES not accepted.')
  if(!is.integer(N_EXAMS)) stop('N_EXAMS not accepted.')

  # Check dimensions CR inputs
  for(outcome in c('Dropout', 'Transfer', 'Graduation')){
    sapply(c('Slope_ability', 'Slope_speed'), function(x) {if(!(x %in% names(PARAMS_LIST[['CR']][[outcome]]))) stop(paste0('Parameters block ', x, ' not detected for ', outcome, '.'))})
    if(length(PARAMS_LIST[['CR']][[outcome]][['Slope_covariates']])!= DIM_EXT) {stop(paste0('Slope covariates dimension does not match DIM_EXT for ', outcome, '.'))}

    if(outcome %in% c('Dropout', 'Transfer')){
      if(length(PARAMS_LIST[['CR']][[outcome]][['Intercepts_year']])!= NYB) {stop(paste0('Intercepts_year does not match NYB for ', outcome, '.'))}
      } else if (outcome == 'Graduation'){
          if(length(PARAMS_LIST[['CR']][[outcome]][['Intercepts_year']])!= NYA) {stop(paste0('Intercepts_year does not match NYA for ', outcome, '.'))}
          }
    }

  # Check dimensions IRT inputs
  if(length(PARAMS_LIST[['IRT']][['Exams_slopes']])!= N_EXAMS) stop(paste0('The number of exams slopes does not mathch N_EXAMS.'))
  if(length(PARAMS_LIST[['IRT']][['Exams_average_time']])!= N_EXAMS) stop(paste0('The number of exams average times does not mathch N_EXAMS.'))
  if(length(PARAMS_LIST[['IRT']][['Exams_variability_time']])!= N_EXAMS) stop(paste0('The number of exams time-variance parameters does not mathch N_EXAMS.'))
  if(length(PARAMS_LIST[['IRT']][['Exams_grades_intercepts']])!= N_EXAMS*N_GRADES) stop(paste0('The number of exams/grades intercepts does not mathch N_EXAMS*N_GRADES.'))

  # Check LAT inputs
  if((length(PARAMS_LIST[['LAT']][['Corr']])>1) | !is.numeric(PARAMS_LIST[['LAT']][['Corr']])) stop('Latent correlation parameter must be scalar')
  if((length(PARAMS_LIST[['LAT']][['Speed_variability']])>1) | !is.numeric(PARAMS_LIST[['LAT']][['Speed_variability']])) stop('Speed_variability must be scalar')


  # Construct the parameter vector
  theta_cr <- c(PARAMS_LIST[['CR']][['Dropout']][['Slope_ability']],
                PARAMS_LIST[['CR']][['Dropout']][['Slope_speed']],
                PARAMS_LIST[['CR']][['Dropout']][['Slope_covariates']],
                PARAMS_LIST[['CR']][['Transfer']][['Slope_ability']],
                PARAMS_LIST[['CR']][['Transfer']][['Slope_speed']],
                PARAMS_LIST[['CR']][['Transfer']][['Slope_covariates']],
                PARAMS_LIST[['CR']][['Graduation']][['Slope_ability']],
                PARAMS_LIST[['CR']][['Graduation']][['Slope_speed']],
                PARAMS_LIST[['CR']][['Graduation']][['Slope_covariates']],
                PARAMS_LIST[['CR']][['Dropout']][['Intercepts_year']],
                PARAMS_LIST[['CR']][['Transfer']][['Intercepts_year']],
                PARAMS_LIST[['CR']][['Graduation']][['Intercepts_year']])
  theta_irt <- c(PARAMS_LIST[['IRT']][['Exams_slopes']],
                 PARAMS_LIST[['IRT']][['Exams_grades_intercepts']],
                 PARAMS_LIST[['IRT']][['Exams_average_time']],
                 PARAMS_LIST[['IRT']][['Exams_variability_time']])
  theta_lat <- c(PARAMS_LIST[['LAT']][['Corr']], PARAMS_LIST[['LAT']][['Speed_variability']])
  theta <- c(theta_cr, theta_irt, theta_lat)

  return(theta)
}

#' Construct named list from parameter vector
#'
#' @param THETA Parameter vector. Structure to be documented.
#' @param DIM_EXT number of external covariates in the competing risk model.
#' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param N_GRADES number of grades modelled.
#' @param N_EXAMS number of exams.
#' @param CONSTRFLAG TRUE if input is defined in the constrained parameter space.
#'
#' @export
paramsVec2list <- function(THETA, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS, CONSTRFLAG = T){

  # Check input type
  if(!is.numeric(THETA)) stop('`THETA` not accepted.')
  if(!is.integer(DIM_EXT)) stop('`DIM_EXT` not accepted.')
  if(!is.integer(NYB)) stop('`NYB` not accepted.')
  if(!is.integer(NYA)) stop('`NYA` not accepted.')
  if(!is.integer(N_GRADES)) stop('`N_GRADES` not accepted.')
  if(!is.integer(N_EXAMS)) stop('`N_EXAMS` not accepted.')

  # Check dimension
  dim <- (3*(DIM_EXT+2) + 2*(NYB) + NYA) +
    3*N_EXAMS + N_EXAMS*N_GRADES + 2
  if(dim!=length(THETA)) stop('THETA dimension does not match with info provided by `DIM_EXT`, `NYB`, `NYA`, `N_GRADES`, `N_EXAMS`')

  # split params
  theta_cr <- THETA[1:(3*(DIM_EXT+2) + 2*(NYB) + NYA)]
  theta_irt <- THETA[(3*(DIM_EXT+2) + 2*(NYB) + NYA+1):(3*(DIM_EXT+2) + 2*(NYB) + NYA + 3*N_EXAMS + N_EXAMS*N_GRADES + 2)]
  theta_lat <- THETA[(3*(DIM_EXT+2) + 2*(NYB) + NYA + 3*N_EXAMS + N_EXAMS*N_GRADES + 1):(3*(DIM_EXT+2) + 2*(NYB) + NYA +
                                                            3*N_EXAMS + N_EXAMS*N_GRADES + 2)]

  # Construct list
  params_list <- list()
  params_list[['CR']] <- list()
  params_list[['IRT']] <- list()
  params_list[['LAT']] <- list()

  # CR parameters
  beta_d <- extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 1)
  beta_t <- extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 2)
  beta_g <- extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 3)

  params_list[['CR']][['Dropout']] <- list(
    'Slope_ability' = beta_d[1],
    'Slope_speed' = beta_d[2],
    'Slope_covariates' = beta_d[3:(2+DIM_EXT)],
    'Intercepts_year' = extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 4)
  )
  params_list[['CR']][['Transfer']] <- list(
    'Slope_ability' = beta_t[1],
    'Slope_speed' = beta_t[2],
    'Slope_covariates' = beta_t[3:(2+DIM_EXT)],
    'Intercepts_year' = extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 5)
  )
  params_list[['CR']][['Graduation']] <- list(
    'Slope_ability' = beta_g[1],
    'Slope_speed' = beta_g[2],
    'Slope_covariates' = beta_g[3:(2+DIM_EXT)],
    'Intercepts_year' = extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 6)
  )

  # IRT parameters
  params_list[['IRT']][['Exams_slopes']] <- unlist(lapply(1:N_EXAMS,
                                                   function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                  N_GRADES = N_GRADES,
                                                                                  N_EXAMS = N_EXAMS,
                                                                                  OPTION = 1, EXAM = x)))
  params_list[['IRT']][['Exams_grades_intercepts']] <- unlist(lapply(1:N_EXAMS,
                                                                     function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                                    N_GRADES = N_GRADES,
                                                                                                    N_EXAMS = N_EXAMS,
                                                                                                    OPTION = 2, EXAM = x)))
  params_list[['IRT']][['Exams_average_time']] <- unlist(lapply(1:N_EXAMS,
                                                         function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                        N_GRADES = N_GRADES,
                                                                                        N_EXAMS = N_EXAMS,
                                                                                        OPTION = 3, EXAM = x)))
  params_list[['IRT']][['Exams_variability_time']] <- unlist(lapply(1:N_EXAMS,
                                                             function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                            N_GRADES = N_GRADES,
                                                                                            N_EXAMS = N_EXAMS,
                                                                                            OPTION = 4, EXAM = x)))
  params_list[['LAT']][['Corr']] <- theta_lat[1]
  params_list[['LAT']][['Speed_variability']] <- theta_lat[2]
  params_list
}
