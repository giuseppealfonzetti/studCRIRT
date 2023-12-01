#ifndef crMod_H
#define crMod_H
#include "extractParams.h"

//' Evaluate hazard function based on outcome and year
//'
//' @param OUTCOME 1 for dropout, 2 for transfer, 3 for graduation
//' @param YEAR Possible values 1:NYB in case of dropout/transfer; 1:NYA in case of graduation
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external covariates
//' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts.
//'
//' @returns It returns the hazard probability of the specific outcome and year.
//'
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

  if(OUTCOME == 1 | OUTCOME == 2){
    const double int1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 4)(YEAR-1);
    const Eigen::VectorXd beta1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 1);
    const double expEta1 = exp(int1 + beta1.dot(COVARIATES));

    const double int2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 5)(YEAR-1);
    const Eigen::VectorXd beta2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 2);
    const double expEta2 = exp(int2 + beta2.dot(COVARIATES));

    if(OUTCOME == 1){
      out = expEta1 /(1+expEta1+expEta2);
    }else if(OUTCOME == 2){
      out = expEta2 /(1+expEta1+expEta2);
    }
  }

  if(OUTCOME == 3){
    const double int3 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 6)(YEAR-1);
    const Eigen::VectorXd beta3 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 3);
    const double expEta3 = exp(int3 + beta3.dot(COVARIATES));

    out = expEta3 /(1+expEta3);
  }

  return(out);
}

//' Evaluate survival function given a the range of years of interest
//'
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param REGIMEFLAG `TRUE` if graduation is possible, `FALSE` otherwise.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//'
//' @returns It returns the probability of survival from `YEAR FIRST` to `YEAR_LAST` included.
//'
//' @export
// [[Rcpp::export]]
double survival(
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::VectorXd THETA_CR,
    Eigen::VectorXd COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM = 100
){

  double out = 1;
  if(YEAR_LAST_EXAM > YEAR_LAST){

    // Regime where graduation is not possible
    if(YEAR_LAST > NYB) Rcpp::stop("Regime 1 mismatch: `YEAR_LAST` > `NYB`");

    for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
      out *= 1 - hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - hazard(2, year, THETA_CR, COVARIATES, NYB, NYA);
    }

  }else if(YEAR_LAST_EXAM <= YEAR_FIRST){

    // Regime where only graduation is possible
    if((YEAR_LAST - YEAR_LAST_EXAM + 1) > NYA) Rcpp::stop("Regime 2 mismatch: `YEAR_LAST` > `YEAR_LAST_EXAM` + `NYA` + 1");

    for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
      out *= 1 - hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA);
    }

  }else if(YEAR_LAST_EXAM <= YEAR_LAST & YEAR_LAST_EXAM > YEAR_FIRST){

    // Regime where graduation is possible from year `YEAR_LAST_EXAM`
    if((YEAR_LAST - YEAR_LAST_EXAM + 1) > NYA) Rcpp::stop("Mixed regime mismatch: `YEAR_LAST` > `YEAR_LAST_EXAM` + `NYA` + 1");

    for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
      out *= 1 - hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - hazard(2, year, THETA_CR, COVARIATES, NYB, NYA);
    }
    for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
      out *= 1 - hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA);
    }

  }

  return(out);
 }

//' Evaluate Outcome Likelihood
//'
//' @param OUTCOME 1 for dropout, 2 for transfer, 3 for graduation
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param OBSFLAG `TRUE` if the outcome is observed, `FALSE` otherwise.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time
//'
//' @export
// [[Rcpp::export]]
double outcomeLik(
    const unsigned int OUTCOME,
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    const bool OBSFLAG,
    Eigen::VectorXd THETA_CR,
    Eigen::VectorXd COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM = 100
){
  double out;

  if(OBSFLAG){
    out = survival(YEAR_FIRST, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM);
  }else{
    out  = survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM);
    out *= hazard(OUTCOME, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA);
  }

  return out;
}
#endif
