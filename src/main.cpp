#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
#include "extractParams.h"
#include "crMod.h"
#include "irtMod.h"
#include "latMod.h"

//' Evaluate the integrand function
//'
//' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
//' @param EXTCOVARIATES External covariates.
//' @param EXAMS_GRADES Vector of grades.
//' @param EXAMS_DAYS Vector of times.
//' @param EXAMS_OBSFLAG Vector of booleans.
//' `TRUE` elements represent observed exams.
//' `FALSE` elements the unobserved ones.
//' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//' @param YEAR Year of evaluation.
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams modelled
//' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
//' @param ABILITY Ability value.
//' @param SPEED Speed value.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the value of the integrand function,
//' given the parameters and the data of a single observation.
//'
//' @export
// [[Rcpp::export]]
double integrand(
  Eigen::VectorXd& THETA,
  Eigen::VectorXd& EXTCOVARIATES,
  std::vector<unsigned int>& EXAMS_GRADES,
  std::vector<double>& EXAMS_DAYS,
  std::vector<bool>& EXAMS_OBSFLAG,
  const unsigned int OUTCOME,
  const unsigned int YEAR,
  const unsigned int N_GRADES,
  const unsigned int N_EXAMS,
  const unsigned int NYB,
  const unsigned int NYA,
  const double ABILITY,
  const double SPEED,
  const unsigned int YEAR_LAST_EXAM = 100,
  const bool LOGFLAG = false
){
  const unsigned int dim_ext = EXTCOVARIATES.size();
  const unsigned int dim_cr = 3*(dim_ext+2) + 2*(NYB) + NYA;
  const unsigned int dim_irt = 3*N_EXAMS + N_EXAMS*N_GRADES;

  Eigen::VectorXd theta_cr = THETA.segment(0, dim_cr);
  Eigen::VectorXd theta_irt = THETA.segment(dim_cr, dim_irt);

  Eigen::VectorXd covariates(2+EXTCOVARIATES.size());
  covariates << ABILITY, SPEED, EXTCOVARIATES;

  double reparlatcorr = THETA(dim_cr+dim_irt+1);
  double reparspeedva = THETA(dim_cr+dim_irt+2);


  double pLat = latent_distr(ABILITY, SPEED, reparlatcorr, reparspeedva, LOGFLAG);
  double out = pLat;



  std::vector<double> pExams(N_EXAMS);
  for(unsigned int exam = 1; exam <= N_EXAMS; exam++){
    pExams[exam-1]= examLik(
      exam,
      EXAMS_GRADES[exam-1],
      EXAMS_DAYS[exam-1],
      EXAMS_OBSFLAG[exam-1],
      theta_irt,
      N_GRADES,
      N_EXAMS,
      ABILITY,
      SPEED,
      LOGFLAG
    );
    if(LOGFLAG){
      out+=pExams[exam-1];
    }else{
      out*=pExams[exam-1];
    }
  }

  double pCR = outcomeLik(OUTCOME, 1, YEAR, theta_cr, covariates, NYB, NYA, YEAR_LAST_EXAM, LOGFLAG);
  if(LOGFLAG){
    out+=pCR;
  }else{
    out*=pCR;
  }


  // // output list
  // Rcpp::List output =
  //   Rcpp::List::create(
  //     Rcpp::Named("dim_ext") = dim_ext,
  //     Rcpp::Named("dim_cr") = dim_cr,
  //     Rcpp::Named("dim_irt") = dim_irt,
  //     Rcpp::Named("theta_cr") = theta_cr,
  //     Rcpp::Named("theta_irt") = theta_irt,
  //     Rcpp::Named("pLat") = pLat,
  //     Rcpp::Named("pExams") = pExams,
  //     Rcpp::Named("pCR") = pCR,
  //     Rcpp::Named("out") = out
  // );
  return out;
}

