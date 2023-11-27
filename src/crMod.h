#ifndef crMod_H
#define crMod_H
#include "extractParams.h"

//' Evaluate hazard function based on outcome and year
//' @param OUTCOME 1 for dropout, 2 for transfer, 3 for graduation
//' @param YEAR Possible values 1:NYB in case of dropout/transfer; 1:NYA in case of graduation
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES First 2 values for ability and speed. Remaining values are external covariates
//' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts.
//' @returns It returns the hazard probability of the specific outcome and year.
//' @export
// [[Rcpp::export]]
double hazard(
  const unsigned int OUTCOME,
  const unsigned int YEAR,
  Eigen::VectorXd THETA_CR,
  Eigen::VectorXd COVARIATES,
  const unsigned int NYB,
  const unsigned int NYA
){

  double out;
  if((OUTCOME == 1 | OUTCOME == 2) & YEAR > NYB) Rcpp::stop("`YEAR` larger than `NYB`");
  if(OUTCOME == 3 & YEAR > NYA) Rcpp::stop("`YEAR` larger than `NYA`");

  // double expEta;
  // double beta0;
  // Eigen::VectorXd beta;

  if(OUTCOME == 1 | OUTCOME == 2){
    double int1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 4)(YEAR-1);
    Eigen::VectorXd beta1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 1);
    double expEta1 = exp(int1 + beta1.dot(COVARIATES));

    double int2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 5)(YEAR-1);
    Eigen::VectorXd beta2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 2);
    double expEta2 = exp(int2 + beta2.dot(COVARIATES));

    if(OUTCOME == 1){
      // beta0 = int1;
      // beta = beta1;
      // expEta = expEta1;
      out = expEta1 /(1+expEta1+expEta2);
    }else if(OUTCOME == 2){
      // beta0 = int2;
      // beta = beta2;
      // expEta = expEta2;
      out = expEta2 /(1+expEta1+expEta2);
    }
  }

  if(OUTCOME == 3){
    double int3 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 6)(YEAR-1);
    Eigen::VectorXd beta3 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 3);
    double expEta3 = exp(int3 + beta3.dot(COVARIATES));

    // beta0 = int3;
    // beta = beta3;
    // expEta = expEta3;
    out = expEta3 /(1+expEta3);
  }

    // output list
    // Rcpp::List output =
    //   Rcpp::List::create(
    //     Rcpp::Named("beta0") = beta0,
    //     Rcpp::Named("beta") = beta,
    //     Rcpp::Named("expEta") = expEta,
    //     Rcpp::Named("hazard") = out
    //   );
  return(out);
}


#endif
