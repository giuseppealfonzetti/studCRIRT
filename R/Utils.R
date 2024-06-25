#' @export
stable_inversion <- function(MAT){
  if(!Matrix::isSymmetric(MAT)){
    warning("Forcing symmetry")
    MAT <- Matrix::forceSymmetric(MAT)
  }

  eig <- eigen(MAT, symmetric = TRUE)
  if(!all(eig$values>0)){
    warning("Negative eigen values detected. Using pracma::_nearest_spd approximation")
    MAT <- pracma::nearest_spd(MAT)
    eig <- eigen(MAT, symmetric = TRUE)
  }

  out <- Matrix::forceSymmetric(
          eig$vectors %*%
          Matrix::Diagonal(x=1/eig$values) %*%
          t(eig$vectors)
    )


  return(as.matrix(out))
}


#' lik_sanity_check(
#'   THETA,
#'   EXTCOVARIATES,
#'   EXAMS_GRADES,
#'   EXAMS_DAYS,
#'   EXAMS_OBSFLAG,
#'   EXAMS_SET,
#'   OUTCOME,
#'   YEAR,
#'   N_GRADES,
#'   N_EXAMS,
#'   NYB,
#'   NYA,
#'   YEAR_LAST_EXAM,
#'   LOGFLAG,
#'   METHOD,
#'   GRID
#' )
