% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{grl_outcomeLik}
\alias{grl_outcomeLik}
\title{Evaluate Outcome Likelihood}
\usage{
grl_outcomeLik(
  OUTCOME,
  YEAR_FIRST,
  YEAR_LAST,
  THETA_CR,
  COVARIATES,
  NYB,
  NYA,
  YEAR_LAST_EXAM = 100L
)
}
\arguments{
\item{OUTCOME}{\code{1} for dropout, \code{2} for transfer, \code{3} for graduation. \code{0} if no outcome is observed.}

\item{YEAR_FIRST}{First year to evaluate.}

\item{YEAR_LAST}{Last year to evaluate.}

\item{THETA_CR}{Portion of the parameter vector related to the competing risk model}

\item{COVARIATES}{The first 2 values refers to ability and speed respectively. Remaining values are external predictors.}

\item{NYB}{Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.}

\item{NYA}{Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.}

\item{YEAR_LAST_EXAM}{Year at which the all exams are completed for the first time.}

\item{LOGFLAG}{Set TRUE to return log value.}
}
\description{
Evaluate Outcome Likelihood
}
