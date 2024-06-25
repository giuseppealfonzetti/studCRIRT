#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen, RcppNumerical)]]
#include "extractParams.h"
#include "crMod.h"
#include "irtMod.h"
#include "latMod.h"

// //' Evaluate the complete data likelihood function
// //'
// //' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
// //' @param EXTCOVARIATES External covariates.
// //' @param EXAMS_GRADES Vector of grades.
// //' @param EXAMS_DAYS Vector of times.
// //' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
// //' @param EXAMS_SET Vector filled with booleans.`TRUE` elements represent exams in the study plan. `FALSE` elements non-relevant ones.
// //' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
// //' @param YEAR Year of evaluation.
// //' @param N_GRADES Number of grades modelled.
// //' @param N_EXAMS Number of exams modelled
// //' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
// //' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
// //' @param ABILITY Ability value.
// //' @param SPEED Speed value.
// //' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
// //'
// //' @returns It returns the value of the integrand function,
// //' given the parameters and the data of a single observation.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::List neg_compl_loglik(
//      Eigen::VectorXd& THETA,
//      Eigen::VectorXd& EXTCOVARIATES,
//      std::vector<unsigned int>& EXAMS_GRADES,
//      std::vector<double>& EXAMS_DAYS,
//      std::vector<bool>& EXAMS_OBSFLAG,
//      std::vector<bool>& EXAMS_SET,
//      const unsigned int OUTCOME,
//      const unsigned int YEAR,
//      const unsigned int N_GRADES,
//      const unsigned int N_EXAMS,
//      const unsigned int NYB,
//      const unsigned int NYA,
//      const double ABILITY,
//      const double SPEED,
//      const unsigned int YEAR_LAST_EXAM = 100
// ){
//
//   const unsigned int dim_ext = EXTCOVARIATES.size();
//   const unsigned int dim_cr = 3*(dim_ext+2) + 2*(NYB) + NYA;
//   const unsigned int dim_irt = 3*N_EXAMS + N_EXAMS*N_GRADES;
//
//   Eigen::VectorXd theta_cr = THETA.segment(0, dim_cr);
//   Eigen::VectorXd theta_irt = THETA.segment(dim_cr, dim_irt);
//
//   Eigen::VectorXd covariates(2+EXTCOVARIATES.size());
//   covariates << ABILITY, SPEED, EXTCOVARIATES;
//
//   double reparlatcorr = THETA(dim_cr+dim_irt);
//   double reparspeedva = THETA(dim_cr+dim_irt+1);
//
//
//   double pLat = latent_distr(ABILITY, SPEED, reparlatcorr, reparspeedva, true);
//   double out = pLat;
//
//
//
//   std::vector<double> pExams(N_EXAMS);
//   for(unsigned int exam = 1; exam <= N_EXAMS; exam++){
//     if(EXAMS_SET[exam-1]==true){
//       pExams[exam-1] = examLik(exam,EXAMS_GRADES[exam-1],EXAMS_DAYS[exam-1], double(YEAR*365),
//                                EXAMS_OBSFLAG[exam-1], theta_irt, N_GRADES, N_EXAMS,
//                                ABILITY, SPEED, true);
//     }else{
//       pExams[exam-1]=0;
//     }
//     out+=pExams[exam-1];
//
//   }
//
//   double pCR = outcomeLik(OUTCOME, 1, YEAR, theta_cr, covariates, NYB, NYA, YEAR_LAST_EXAM, true);
//   out+=pCR;
//
//   Rcpp::List output =
//     Rcpp::List::create(
//       Rcpp::Named("nll") = -out,
//       Rcpp::Named("pexams") = pExams,
//       Rcpp::Named("pCR") = pCR,
//       Rcpp::Named("pLat") = pLat
//
//     );
//   return(output);
// }
//
// //' Evaluate the complete data likelihood function
//  //'
//  //' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
//  //' @param EXTCOVARIATES External covariates.
//  //' @param EXAMS_GRADES Vector of grades.
//  //' @param EXAMS_DAYS Vector of times.
//  //' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
//  //' @param EXAMS_SET Vector filled with booleans.`TRUE` elements represent exams in the study plan. `FALSE` elements non-relevant ones.
//  //' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//  //' @param YEAR Year of evaluation.
//  //' @param N_GRADES Number of grades modelled.
//  //' @param N_EXAMS Number of exams modelled
//  //' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//  //' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
//  //' @param ABILITY Ability value.
//  //' @param SPEED Speed value.
//  //' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//  //'
//  //' @returns It returns the value of the integrand function,
//  //' given the parameters and the data of a single observation.
//  //'
//  //' @export
//  // [[Rcpp::export]]
//  Rcpp::List IRT_lat_loglik(
//      Eigen::VectorXd& THETA,
//      Eigen::VectorXd& EXTCOVARIATES,
//      std::vector<unsigned int>& EXAMS_GRADES,
//      std::vector<double>& EXAMS_DAYS,
//      std::vector<bool>& EXAMS_OBSFLAG,
//      std::vector<bool>& EXAMS_SET,
//      const unsigned int OUTCOME,
//      const unsigned int YEAR,
//      const unsigned int N_GRADES,
//      const unsigned int N_EXAMS,
//      const unsigned int NYB,
//      const unsigned int NYA,
//      const double ABILITY,
//      const double SPEED,
//      const unsigned int YEAR_LAST_EXAM = 100
//  ){
//
//    const unsigned int dim_ext = EXTCOVARIATES.size();
//    const unsigned int dim_cr = 3*(dim_ext+2) + 2*(NYB) + NYA;
//    const unsigned int dim_irt = 3*N_EXAMS + N_EXAMS*N_GRADES;
//
//    Eigen::VectorXd theta_cr = THETA.segment(0, dim_cr);
//    Eigen::VectorXd theta_irt = THETA.segment(dim_cr, dim_irt);
//
//
//    double reparlatcorr = THETA(dim_cr+dim_irt);
//    double reparspeedva = THETA(dim_cr+dim_irt+1);
//
//
//    double pLat = latent_distr(ABILITY, SPEED, reparlatcorr, reparspeedva, true);
//    double out = pLat;
//
//
//
//    std::vector<double> pExams(N_EXAMS);
//    for(unsigned int exam = 1; exam <= N_EXAMS; exam++){
//      if(EXAMS_SET[exam-1]==true){
//        pExams[exam-1] = examLik(exam,EXAMS_GRADES[exam-1],EXAMS_DAYS[exam-1], double(YEAR*365),
//                                 EXAMS_OBSFLAG[exam-1], theta_irt, N_GRADES, N_EXAMS,
//                                 ABILITY, SPEED, true);
//      }else{
//        pExams[exam-1]=0;
//      }
//      out+=pExams[exam-1];
//
//    }
//
//
//    Rcpp::List output =
//      Rcpp::List::create(
//        Rcpp::Named("nll") = -out,
//        Rcpp::Named("pexams") = pExams,
//        Rcpp::Named("pLat") = pLat
//      );
//    return(output);
//  }


//' Full model
//'
//' Evaluate the full model log likelihood of a single observation
//'
//' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
//' @param EXTCOVARIATES External covariates.
//' @param EXAMS_GRADES Vector of grades.
//' @param EXAMS_DAYS Vector of times.
//' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
//' @param EXAMS_SET Vector filled with booleans.`TRUE` elements represent exams in the study plan. `FALSE` elements non-relevant ones.
//' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//' @param YEAR Year of evaluation.
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams modelled
//' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
//' @param ABILITY Ability value.
//' @param SPEED Speed value.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//' @param GRFLAG `TRUE` to compute the gradient.
//'
//' @return It returns the value of the integrand function,
//' given the parameters and the data of a single observation.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List full_model(
    Eigen::VectorXd& THETA,
    Eigen::VectorXd& EXTCOVARIATES,
    Eigen::VectorXd& EXAMS_GRADES,
    Eigen::VectorXd& EXAMS_DAYS,
    Eigen::VectorXd& EXAMS_OBSFLAG,
    Eigen::VectorXd& EXAMS_SET,
    const unsigned int OUTCOME,
    const unsigned int YEAR,
    const unsigned int YEAR_LAST_EXAM,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const unsigned int NYB,
    const unsigned int NYA,
    const double ABILITY,
    const double SPEED,
    const bool GRFLAG = true
){

  const unsigned int dim_ext = EXTCOVARIATES.size();
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());

  // Initialize latent distribution
  LAT_DISTR lat_distr(THETA);

  // Initialize conditional IRT model
  IRT_MOD irt_mod(THETA, EXAMS_GRADES, EXAMS_DAYS, EXAMS_OBSFLAG,
                  EXAMS_SET, YEAR, N_GRADES, N_EXAMS,
                  NYB, NYA, dim_ext);

  // Initialize conditional CR model
  CR_MOD cr_mod(THETA, EXTCOVARIATES, OUTCOME, YEAR, N_GRADES,
                N_EXAMS, NYB, NYA, YEAR_LAST_EXAM);

  double ll = lat_distr.ll(ABILITY, SPEED) + irt_mod.ll(ABILITY, SPEED) + cr_mod.ll(ABILITY, SPEED);

  if(GRFLAG){
    gr += irt_mod.grll(ABILITY, SPEED) + cr_mod.grll(ABILITY, SPEED);
    gr.tail(2) += lat_distr.grll(ABILITY, SPEED);
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = gr,
      Rcpp::Named("ll") = ll
    );

  return output;
}


//' Full model
//'
//' Evaluate the full model log likelihood of a single observation
//'
//' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
//' @param EXTCOVARIATES External covariates.
//' @param EXAMS_GRADES Vector of grades.
//' @param EXAMS_DAYS Vector of times.
//' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
//' @param EXAMS_SET Vector filled with booleans.`TRUE` elements represent exams in the study plan. `FALSE` elements non-relevant ones.
//' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//' @param YEAR Year of evaluation.
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams modelled
//' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
//' @param ABILITY Ability value.
//' @param SPEED Speed value.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//' @param GRFLAG `TRUE` to compute the gradient.
//'
//' @return It returns the value of the integrand function,
//' given the parameters and the data of a single observation.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List full_GH_sample(
    Eigen::VectorXd& THETA,
    Eigen::MatrixXd& EXTCOVARIATES,
    Eigen::MatrixXd& EXAMS_GRADES,
    Eigen::MatrixXd& EXAMS_DAYS,
    Eigen::MatrixXd& EXAMS_OBSFLAG,
    Eigen::MatrixXd& EXAMS_SET,
    Eigen::VectorXd& OUTCOME,
    Eigen::VectorXd& YEAR,
    Eigen::VectorXd& YEAR_LAST_EXAM,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd WEIGHTS,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const unsigned int NYB,
    const unsigned int NYA,
    const bool GRFLAG = true
){

  const unsigned int dim_ext = EXTCOVARIATES.cols();

  double ll = 0;
  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

  const unsigned int n = EXAMS_GRADES.rows();
  const unsigned int nq = GRID.rows();

   for(unsigned int i = 0; i < n; i++){
     // Initialize latent distribution
     LAT_DISTR lat_distr(THETA);

     // Initialize conditional IRT model
     IRT_MOD irt_mod(THETA, EXAMS_GRADES.row(i), EXAMS_DAYS.row(i),
                     EXAMS_OBSFLAG.row(i), EXAMS_SET.row(i), YEAR(i),
                     N_GRADES, N_EXAMS, NYB, NYA, dim_ext);

     // Initialize conditional CR model
     CR_MOD cr_mod(THETA, EXTCOVARIATES.row(i), OUTCOME(i), YEAR(i),
                   N_GRADES, N_EXAMS, NYB, NYA, YEAR_LAST_EXAM(i));


     // Initialize grid evaluation for ll and its gradient
     Eigen::VectorXd f(nq);
     Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);


     for(unsigned int point = 0; point < nq; point++){
       f(point) = exp(
         lat_distr.ll(GRID(point, 0), GRID(point, 1)) +
           irt_mod.ll(GRID(point, 0), GRID(point, 1)) +
           cr_mod.ll(GRID(point, 0), GRID(point, 1))
       );

       if(GRFLAG){
         Eigen::VectorXd gr_point = irt_mod.grll(GRID(point, 0), GRID(point, 1)) + cr_mod.grll(GRID(point, 0), GRID(point, 1));
         gr_point.tail(2) = lat_distr.grll(GRID(point, 0), GRID(point, 1));
         gr_point *= f(point);
         gr.col(point) = gr_point;
       }

     }

     double lli = std::max(-10000.0, log(f.dot(WEIGHTS)));
     ll += lli;
     if(GRFLAG){
       grll += gr*WEIGHTS/std::max(1e-16,exp(lli));
     }
   }

   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("gr") = grll,
       Rcpp::Named("ll") = ll
     );

   return output;
 }
