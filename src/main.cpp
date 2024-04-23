#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <RcppNumerical.h>

// [[Rcpp::depends(RcppEigen, RcppNumerical)]]
#include "extractParams.h"
#include "crMod.h"
#include "irtMod.h"
#include "latMod.h"
#include "completeLikelihood.h"



//' Evaluate the complete data likelihood function
//'
//' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
//' @param EXTCOVARIATES External covariates.
//' @param EXAMS_GRADES Vector of grades.
//' @param EXAMS_DAYS Vector of times.
//' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
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
double complete_likelihood(
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

   CompleteL CL(THETA,
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
                YEAR_LAST_EXAM,
                LOGFLAG);

   Eigen::VectorXd latent(2);
   latent << ABILITY, SPEED;
   double out = CL(latent);
   return out;
 }


