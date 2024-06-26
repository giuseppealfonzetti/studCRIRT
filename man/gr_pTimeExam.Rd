% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gr_pTimeExam}
\alias{gr_pTimeExam}
\title{Evaluate the c.d.f or p.d.f of the last attempt to an exam}
\usage{
gr_pTimeExam(
  EXAM,
  DAY,
  THETA_IRT,
  N_GRADES,
  N_EXAMS,
  SPEED,
  CDFFLAG,
  LOGFLAG = FALSE
)
}
\arguments{
\item{EXAM}{Exam of interest.}

\item{DAY}{Day of interest.}

\item{THETA_IRT}{Portion of the parameter vector related to the IRT model}

\item{N_GRADES}{Number of grades modelled.}

\item{N_EXAMS}{Number of exams.}

\item{SPEED}{speed value.}

\item{CDFFLAG}{\code{TRUE} for c.d.f. of time. \code{FALSE} for p.d.f.}

\item{LOGFLAG}{Set TRUE to return log value.}
}
\description{
Evaluate the c.d.f or p.d.f of the last attempt to an exam
}
