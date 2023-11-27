#ifndef irtMod_H
#define irtMod_H
#include "extractParams.h"

//' Evaluate the probability of grades greater or equal than the reference one
//' @param GRADE Grade used as reference
//' @param EXAM Exam of interest
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY Ability value.
//'
//' @returns It returns the probability of obtaining grades higher than `GRADE` on exam `EXAM`.
//' @export
// [[Rcpp::export]]
double pGreaterGrades(
  const unsigned int GRADE,
  const unsigned int EXAM,
  Eigen::VectorXd THETA_IRT,
  const unsigned int N_GRADES,
  const unsigned int N_EXAMS,
  const double ABILITY
){
  if(EXAM > N_EXAMS) Rcpp::stop("`EXAM` larger than `N_EXAMS`");
  if(GRADE > N_GRADES) Rcpp::stop("`GRADE` larger than `N_GRADES`");

  const double intercept = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM)(GRADE-1);
  const double coeff = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM)(0);
  const double expEta = exp(intercept + coeff*ABILITY);
  const double prob = expEta/(1+expEta);

  // // output list
  // Rcpp::List output =
  //   Rcpp::List::create(
  //     Rcpp::Named("intercept") = intercept,
  //     Rcpp::Named("coef") = coeff,
  //     Rcpp::Named("prob") = prob,
  //     // Rcpp::Named("expEta") = expEta,
  //     Rcpp::Named("out") = out
  //   );

  return(prob);
}

//' Evaluate the probability of getting a specific grade
//'
//' @param GRADE Grade used as reference
//' @param EXAM Exam of interest
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY Ability value.
//'
//' @returns It returns the probability of obtaining the grade `GRADE` on exam `EXAM`.
//' @export
// [[Rcpp::export]]
double pGrade(
    const unsigned int GRADE,
    const unsigned int EXAM,
    Eigen::VectorXd THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY
){
  double out;
  if(GRADE==N_GRADES){
    out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
  }else if(GRADE==0){
    out = 1 - pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
  }else if(GRADE<N_GRADES & GRADE >0){
    out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY)-pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
  }

  return(out);
}

//' Evaluate the probability of last attempt to an exam to occur before a certain day
//'
//' @param EXAM Exam of interest.
//' @param DAY Day of interest.
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param SPEED speed value.

//' @export
// [[Rcpp::export]]
double pTimeExam(
    const unsigned int EXAM,
    const double DAY,
    Eigen::VectorXd THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double SPEED
){
  Eigen::VectorXd pars = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 3, EXAM);
  const double mean = pars(0) - SPEED;
  const double sd = 1/pars(1);
  const double out = R::plnorm(DAY, mean, sd, true, false);
  return(out);
}
#endif
