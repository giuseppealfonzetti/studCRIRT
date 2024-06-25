#' Intercepts reparameterisation
#'
#' @description
#' It allows to go back and forth from constrained intercepts
#' to unconstrained parameters.
#'
#' @param X Intercepts or unconstrained parameter values for grades low-to-high.
#' @param CON2UN TRUE if going from the constrained space to the unconstrained one.
#'
#' @returns The numeric vector containing the constrained intercepts
#' or the unconstrained parameters.
#'
#' @examples
#' set.seed(123)
#'
#' # random unconstrained parameters
#' params <- rnorm(4)
#'
#' # from unconstrained parameters to constrained intercepts
#' intercepts <- reparInt(params, CON2UN = FALSE)
#'
#' # retrieve params
#' reparInt(intercepts, CON2UN = TRUE)
#'
#' @export
reparInt <- function(X, CON2UN = T){
  out <- c()
  if(!CON2UN){
    # From unconstrained to constrained
    # reverse <- rev(X)
    # param <- reverse[1]
    # tmp <- reverse[-1]
    # tmp <- exp(tmp)
    # param <- c(param, tmp)
    # out <- rev(cumsum(param))
    out <- cppReparInt(X, CON2UN)
  }else{
    if(!(sum(sort(X, decreasing = T) == X)==length(X)))stop('Constrained INTERCEPTS must be decreasing')
    # From constrained to unconstrained
    # out <- c(rev(log(diff(rev(X)))), X[length(X)])
    out <- cppReparInt(X, CON2UN)

  }

  return(out)
}

#' Check if the provided parameter list is genuine
#'
#' @inheritParams paramsList2vec
#' @export
checkerList2vec <- function(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS){

  if(!is.numeric(unlist(PARAMS_LIST)))stop('PARAMS_LIST contains non-numeric values.')
  # Check if all parameter blocks are specified
  sapply(c('CR', 'IRT', 'LAT'), function(x) {if(!(x %in% names(PARAMS_LIST))) stop(paste0('Parameters block ', x, ' not detected'))})
  sapply(c('Dropout', 'Transfer', 'Graduation'), function(x) {if(!(x %in% names(PARAMS_LIST[['CR']]))) stop(paste0('Parameters block ', x, ' not detected'))})
  sapply(c('Exams_slopes', 'Exams_grades_intercepts', 'Exams_average_time', 'Exams_variability_time'), function(x) {if(!(x %in% names(PARAMS_LIST[['IRT']]))) stop(paste0('Parameters block ', x, ' not detected'))})
  sapply(c('Corr', 'Speed_variability'), function(x) {if(!(x %in% names(PARAMS_LIST[['LAT']]))) stop(paste0('Parameters block ', x, ' not detected'))})

  # Check input types
  if(!is.integer(c(DIM_EXT,NYB,NYA, N_GRADES, N_EXAMS))) stop('Some of DIM_EXT,NYB,NYA, N_GRADES, N_EXAMS is not an integer.')

  # Check dimensions CR inputs
  for(outcome in c('Dropout', 'Transfer', 'Graduation')){
    sapply(c('Slope_ability', 'Slope_speed'), function(x) {if(!(x %in% names(PARAMS_LIST[['CR']][[outcome]]))) stop(paste0('Parameters block ', x, ' not detected for ', outcome, '.'))})
    sapply(c('Slope_ability', 'Slope_speed'), function(x) {if((length(PARAMS_LIST[['CR']][[outcome]][[x]])>1)) stop(paste0('Parameters block ', x, ' not valid for ', outcome, '.'))})
    if(length(PARAMS_LIST[['CR']][[outcome]][['Slope_covariates']])!= DIM_EXT) {stop(paste0('Slope covariates dimension does not match DIM_EXT for ', outcome, '.'))}

    if(outcome %in% c('Dropout', 'Transfer')){
      if(length(PARAMS_LIST[['CR']][[outcome]][['Intercepts_year']])!= NYB) {stop(paste0('Intercepts_year does not match NYB for ', outcome, '.'))}
    } else if (outcome == 'Graduation'){
      if(length(PARAMS_LIST[['CR']][[outcome]][['Intercepts_year']])!= NYA) {stop(paste0('Intercepts_year does not match NYA for ', outcome, '.'))}
    }
  }

  # Check dimensions IRT inputs
  if(length(PARAMS_LIST[['IRT']][['Exams_slopes']])!= N_EXAMS) stop(paste0('The number of exams slopes does not match N_EXAMS.'))
  if(length(PARAMS_LIST[['IRT']][['Exams_average_time']])!= N_EXAMS) stop(paste0('The number of exams average times does not mathch N_EXAMS.'))
  if(length(PARAMS_LIST[['IRT']][['Exams_variability_time']])!= N_EXAMS) stop(paste0('The number of exams time-variance parameters does not mathch N_EXAMS.'))
  if(ncol(PARAMS_LIST[['IRT']][['Exams_grades_intercepts']])!= N_GRADES) stop(paste0('The number of columns for exams-grades intercepts does not mathch N_GRADES.'))
  if(nrow(PARAMS_LIST[['IRT']][['Exams_grades_intercepts']])!= N_EXAMS) stop(paste0('The number of rows for exams-grades intercepts does not mathch N_EXAMS'))

  # Check LAT inputs
  if((length(PARAMS_LIST[['LAT']][['Corr']])>1) | !is.numeric(PARAMS_LIST[['LAT']][['Corr']])) stop('Latent correlation parameter not valid')
  if((length(PARAMS_LIST[['LAT']][['Speed_variability']])>1) | !is.numeric(PARAMS_LIST[['LAT']][['Speed_variability']])) stop('Speed_variability not valid')
  if((abs(PARAMS_LIST[['LAT']][['Corr']])>=1) ) stop('Latent correlation parameter not valid')
  if((PARAMS_LIST[['LAT']][['Speed_variability']]<=0)) stop('Speed_variability not valid')
  if(sum(PARAMS_LIST[['IRT']][['Exams_variability_time']]>0)!= N_EXAMS) stop(paste0('The number of exams time-variance parameters does not mathch N_EXAMS.'))

  return(TRUE)
}


#' Construct parameter vector from list
#'
#' @description
#' Given a properly-structered list of parameters,
#' it returns the implied unconstrained parameter vector.
#'
#' @section Internal reparameterisations:
#'
#' - Exam-specific intercepts vector. Decreasing from lower grade to higher.
#' Reparametrised using \link{reparInt}(),
#'
#' - Latent correlation. Reparameterised via Fisher's transformation.
#'
#' - Standard deviations. Reparameterised via exponential transformation.
#'
#'
#' @param PARAMS_LIST Parameter list. Structure to be documented. See the example.
#' @param DIM_EXT number of external covariates in the competing risk model.
#' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts in the competing risk model.
#' @param N_GRADES number of grades modelled.
#' @param N_EXAMS number of exams.
#'
#' @returns It returns the unconstrained parameter vector.
#'
#' @examples
#' seed <- 123
#'
#' ## Setup dimensions
#' external_covariates <- 3L     # number of external covariates
#' years_before <- 6L            # number of possible years in the regime without graduation
#' years_after <- 2L             # number of possible years in the regime with graduation
#' grades <- 4L                  # number of grades
#' exams <- 10L                  # number of exams
#'
#'
#' ## Construct the list
#' params_list <- list()
#'
#' # Three main blocks
#' params_list[['CR']] <- list()
#' params_list[['IRT']] <- list()
#' params_list[['LAT']] <- list()
#'
#' # Parameter-blocks related to competing risks
#' set.seed(seed)
#' params_list[['CR']][['Dropout']] <- list(
#'  'Slope_ability' = rnorm(1),
#'  'Slope_speed' = rnorm(1),
#'  'Slope_covariates' = rnorm(external_covariates),
#'  'Intercepts_year' = rnorm(years_before)
#')
#'
#' set.seed(seed+1)
#' params_list[['CR']][['Transfer']] <- list(
#'  'Slope_ability' = rnorm(1),
#'  'Slope_speed' = rnorm(1),
#'  'Slope_covariates' = rnorm(external_covariates),
#'  'Intercepts_year' = rnorm(years_before)
#' )
#'
#' set.seed(seed+2)
#' params_list[['CR']][['Graduation']] <- list(
#'  'Slope_ability' = rnorm(1),
#'  'Slope_speed' = rnorm(1),
#'  'Slope_covariates' = rnorm(external_covariates),
#'  'Intercepts_year' = rnorm(years_after)
#' )
#'
#' ## Parameter-blocks related to exams IRT-modelling
#' set.seed(seed+3)
#' params_list[['IRT']][['Exams_slopes']] <- runif(exams, 0, 5)
#'
#' set.seed(seed+4)
#' examsInt <- matrix(0, exams, grades)
#' for(ex in 1:exams){
#' examsInt[ex,] <- sort(rnorm(grades, 0, 5), decreasing = TRUE)
#' }
#' params_list[['IRT']][['Exams_grades_intercepts']] <- examsInt
#'
#' set.seed(seed+5)
#' params_list[['IRT']][['Exams_average_time']] <- rnorm(exams)
#'
#' set.seed(seed+6)
#' params_list[['IRT']][['Exams_variability_time']] <- runif(exams)
#'
#' params_list[['LAT']][['Corr']] <- .6
#' params_list[['LAT']][['Speed_variability']] <- 1.5
#'
#' theta <- paramsList2vec(
#'   PARAMS_LIST = params_list,
#'   DIM_EXT = external_covariates,
#'   NYB = years_before,
#'   NYA = years_after,
#'   N_GRADES = grades,
#'   N_EXAMS = exams)
#'
#'
#' @export
paramsList2vec <- function(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS){

  # handle errors
  checkerList2vec(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS)

  # Intercepts reparameterisation
  repExamsGradesInt <-
    as.numeric(unlist(
      lapply(as.list(as.data.frame(t(PARAMS_LIST[['IRT']][['Exams_grades_intercepts']]))),
             reparInt)))

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
                 repExamsGradesInt,
                 # as.vector(t(PARAMS_LIST[['IRT']][['Exams_grades_intercepts']])),
                 PARAMS_LIST[['IRT']][['Exams_average_time']],
                 log(PARAMS_LIST[['IRT']][['Exams_variability_time']]))

  rho <- PARAMS_LIST[['LAT']][['Corr']]
  sig <- PARAMS_LIST[['LAT']][['Speed_variability']]
  Slat <- matrix(c(1,rho*sig,rho*sig,sig^2),2,2)
  L <- t(chol(Slat))
  theta_lat <- c(L[2,1], L[2,2])
  # theta_lat <- c(atanh(PARAMS_LIST[['LAT']][['Corr']]), log(PARAMS_LIST[['LAT']][['Speed_variability']]))
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
#' @param LABS_EXAMS optional label for exams
#' @param LABS_GRADES optional label for grades
#' @param LABS_COVARIATES optional label for covariates
#'
#' @returns Given a proper unconstrained parameter vector, it returns a
#' list of the model parameters.
#'
#' @examples
#'
#' seed <- 123
#'
#' ## Setup dimensions
#' external_covariates <- 3L     # number of external covariates
#' years_before <- 6L            # number of possible years in the regime without graduation
#' years_after <- 2L             # number of possible years in the regime with graduation
#' grades <- 4L                  # number of grades
#' exams <- 10L                  # number of exams
#'
#' dim <- 3 * (external_covariates+2) + 2 * (years_before) +
#'         years_after + 3 * exams + exams * grades + 2
#'
#' set.seed(seed)
#' theta <- rnorm(dim)
#'
#' ## Retrieve the list
#' params_list <- paramsVec2list(
#'   THETA = theta,
#'   DIM_EXT = external_covariates,
#'   NYB = years_before,
#'   NYA = years_after,
#'   N_GRADES = grades,
#'   N_EXAMS = exams)
#'
#'
#'
#' @export
paramsVec2list <- function(THETA, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS, LABS_EXAMS=NULL, LABS_GRADES=NULL, LABS_COVARIATES=NULL){

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
  theta_irt <- THETA[(3*(DIM_EXT+2) + 2*(NYB) + NYA+1):(3*(DIM_EXT+2) + 2*(NYB) + NYA + 3*N_EXAMS + N_EXAMS*N_GRADES )]
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
  names(params_list[['CR']][['Dropout']][['Slope_covariates']]) <- LABS_COVARIATES
  params_list[['CR']][['Transfer']] <- list(
    'Slope_ability' = beta_t[1],
    'Slope_speed' = beta_t[2],
    'Slope_covariates' = beta_t[3:(2+DIM_EXT)],
    'Intercepts_year' = extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 5)
  )
  names(params_list[['CR']][['Transfer']][['Slope_covariates']]) <- LABS_COVARIATES

  params_list[['CR']][['Graduation']] <- list(
    'Slope_ability' = beta_g[1],
    'Slope_speed' = beta_g[2],
    'Slope_covariates' = beta_g[3:(2+DIM_EXT)],
    'Intercepts_year' = extract_params_cr(THETA_CR = theta_cr, DIM_EXT = DIM_EXT, NYB = NYB, NYA = NYA, OPTION = 6)
  )
  names(params_list[['CR']][['Graduation']][['Slope_covariates']]) <- LABS_COVARIATES

  # IRT parameters
  params_list[['IRT']][['Exams_slopes']] <- unlist(lapply(1:N_EXAMS,
                                                   function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                  N_GRADES = N_GRADES,
                                                                                  N_EXAMS = N_EXAMS,
                                                                                  OPTION = 1, EXAM = x)))
  names(params_list[['IRT']][['Exams_slopes']]) <- LABS_EXAMS
  #  Reparametrise intercepts
  # examsGradesInt <- matrix(unlist(lapply(
  #   lapply(1:N_EXAMS,
  #          function(x) extract_params_irt(
  #            THETA_IRT = theta_irt, N_GRADES = N_GRADES, N_EXAMS = N_EXAMS, OPTION = 2, EXAM = x
  #          )
  #   ),
  #   function(x) reparInt(x, CON2UN = F)
  # )), N_EXAMS, N_GRADES, byrow = T)
  examsGradesInt <- matrix(unlist(
    lapply(1:N_EXAMS,
           function(x) extract_params_irt(
             THETA_IRT = theta_irt, N_GRADES = N_GRADES, N_EXAMS = N_EXAMS, OPTION = 2, EXAM = x
           )
    )
  ), N_EXAMS, N_GRADES, byrow = T)
  rownames(examsGradesInt) <- LABS_EXAMS
  colnames(examsGradesInt) <- LABS_GRADES

  params_list[['IRT']][['Exams_grades_intercepts']] <- examsGradesInt
  params_list[['IRT']][['Exams_average_time']] <- unlist(lapply(1:N_EXAMS,
                                                         function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                        N_GRADES = N_GRADES,
                                                                                        N_EXAMS = N_EXAMS,
                                                                                        OPTION = 3, EXAM = x)))
  names(params_list[['IRT']][['Exams_average_time']]) <- LABS_EXAMS
  params_list[['IRT']][['Exams_variability_time']] <- unlist(lapply(1:N_EXAMS,
                                                             function(x) extract_params_irt(THETA_IRT = theta_irt,
                                                                                            N_GRADES = N_GRADES,
                                                                                            N_EXAMS = N_EXAMS,
                                                                                            OPTION = 4, EXAM = x)))
  names(params_list[['IRT']][['Exams_variability_time']]) <- LABS_EXAMS

  L <- matrix(c(1,theta_lat[1],0, theta_lat[2]),2,2)
  Slat <- L %*% t(L)
  params_list[['LAT']][['Corr']] <- Slat[2,1]/sqrt(Slat[2,2])
  params_list[['LAT']][['Speed_variability']] <- sqrt(Slat[2,2])
  params_list
}

#' @export
paramsVec2Chol <- function(THETA, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS, LABS_EXAMS=NULL, LABS_GRADES=NULL, LABS_COVARIATES=NULL){

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
  theta_lat <- THETA[(3*(DIM_EXT+2) + 2*(NYB) + NYA + 3*N_EXAMS + N_EXAMS*N_GRADES + 1):(3*(DIM_EXT+2) + 2*(NYB) + NYA +
                                                                                           3*N_EXAMS + N_EXAMS*N_GRADES + 2)]
  L <- matrix(c(1,theta_lat[1], 0,  theta_lat[2]),2,2)
  return(L)
}


#' @export
paramsList2Ivec <- function(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS){

  # handle errors
  checkerList2vec(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS)


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
                 as.vector(t(PARAMS_LIST[['IRT']][['Exams_grades_intercepts']])),
                 PARAMS_LIST[['IRT']][['Exams_average_time']],
                 PARAMS_LIST[['IRT']][['Exams_variability_time']])

  theta_lat <- c(PARAMS_LIST[['LAT']][['Corr']], PARAMS_LIST[['LAT']][['Speed_variability']])
  theta <- c(theta_cr, theta_irt, theta_lat)

  return(theta)
}

#' @export
paramsIvec2list <- function(THETA, HELPER_LIST) {
  PARAMS_LIST <- list(
    CR = list(
      Dropout = list(Slope_ability = NULL, Slope_speed = NULL, Slope_covariates = NULL, Intercepts_year = NULL),
      Transfer = list(Slope_ability = NULL, Slope_speed = NULL, Slope_covariates = NULL, Intercepts_year = NULL),
      Graduation = list(Slope_ability = NULL, Slope_speed = NULL, Slope_covariates = NULL, Intercepts_year = NULL)
    ),
    IRT = list(
      Exams_slopes = NULL,
      Exams_grades_intercepts = NULL,
      Exams_average_time = NULL,
      Exams_variability_time = NULL
    ),
    LAT = list(Corr = NULL, Speed_variability = NULL)
  )

  # Initialize index
  idx <- 1

  # CR parameters
  for (sublist in c('Dropout', 'Transfer', 'Graduation')) {
    for (param in c('Slope_ability', 'Slope_speed', 'Slope_covariates', 'Intercepts_year')) {
      length_param <- length(HELPER_LIST[['CR']][[sublist]][[param]])
      PARAMS_LIST[['CR']][[sublist]][[param]] <- THETA[idx:(idx + length_param - 1)]
      idx <- idx + length_param
    }
  }

  # IRT parameters
  length_param <- length(HELPER_LIST[['IRT']][['Exams_slopes']])
  PARAMS_LIST[['IRT']][['Exams_slopes']] <- THETA[idx:(idx + length_param - 1)]
  idx <- idx + length_param

  matrix_dims <- dim(HELPER_LIST[['IRT']][['Exams_grades_intercepts']])
  length_param <- prod(matrix_dims)
  PARAMS_LIST[['IRT']][['Exams_grades_intercepts']] <- matrix(THETA[idx:(idx + length_param - 1)], nrow = matrix_dims[1], ncol = matrix_dims[2])
  idx <- idx + length_param

  length_param <- length(HELPER_LIST[['IRT']][['Exams_average_time']])
  PARAMS_LIST[['IRT']][['Exams_average_time']] <- THETA[idx:(idx + length_param - 1)]
  idx <- idx + length_param

  length_param <- length(HELPER_LIST[['IRT']][['Exams_variability_time']])
  PARAMS_LIST[['IRT']][['Exams_variability_time']] <- THETA[idx:(idx + length_param - 1)]
  idx <- idx + length_param

  # LAT parameters
  PARAMS_LIST[['LAT']][['Corr']] <- THETA[length(THETA)-1]
  PARAMS_LIST[['LAT']][['Speed_variability']] <- THETA[length(THETA)]

  return(PARAMS_LIST)
}

