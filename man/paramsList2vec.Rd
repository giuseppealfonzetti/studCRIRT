% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractParams.R
\name{paramsList2vec}
\alias{paramsList2vec}
\title{Construct parameter vector from list}
\usage{
paramsList2vec(PARAMS_LIST, DIM_EXT, NYB, NYA, N_GRADES, N_EXAMS)
}
\arguments{
\item{PARAMS_LIST}{Parameter list. Structure to be documented. See the example.}

\item{DIM_EXT}{number of external covariates in the competing risk model.}

\item{NYB}{number of years in the non-graduatable state. Needed for determining how many time-related intercepts in the competing risk model.}

\item{NYA}{number of years in the graduatable state. Needed for determining how many time-related intercepts in the competing risk model.}

\item{N_GRADES}{number of grades modelled.}

\item{N_EXAMS}{number of exams.}
}
\value{
It returns the unconstrained parameter vector.
}
\description{
Given a properly-structered list of parameters,
it returns the implied unconstrained parameter vector.
}
\section{Internal reparameterisations}{

\itemize{
\item Exam-specific intercepts vector. Decreasing from lower grade to higher.
Reparametrised using \link{reparInt}(),
\item Latent correlation. Reparameterised via Fisher's transformation.
\item Standard deviations. Reparameterised via exponential transformation.
}
}

\examples{
seed <- 123

## Setup dimensions
external_covariates <- 3L     # number of external covariates
years_before <- 6L            # number of possible years in the regime without graduation
years_after <- 2L             # number of possible years in the regime with graduation
grades <- 4L                  # number of grades
exams <- 10L                  # number of exams


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

theta <- paramsList2vec(
  PARAMS_LIST = params_list,
  DIM_EXT = external_covariates,
  NYB = years_before,
  NYA = years_after,
  N_GRADES = grades,
  N_EXAMS = exams)


}
