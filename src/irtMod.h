#ifndef irtMod_H
#define irtMod_H
#include "extractParams.h"

//' Evaluate the probability of grades greater or equal than the reference one
//'
//' @param GRADE Grade used as reference
//' @param EXAM Exam of interest
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY Ability value.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the probability of obtaining grades higher than `GRADE` on exam `EXAM`.
//'
//' @export
// [[Rcpp::export]]
double pGreaterGrades(
  const unsigned int GRADE,
  const unsigned int EXAM,
  Eigen::VectorXd& THETA_IRT,
  const unsigned int N_GRADES,
  const unsigned int N_EXAMS,
  const double ABILITY,
  const bool LOGFLAG = false
){
  if(EXAM > N_EXAMS) Rcpp::stop("`EXAM` larger than `N_EXAMS`");
  if(GRADE > N_GRADES) Rcpp::stop("`GRADE` larger than `N_GRADES`");

  const double intercept = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM)(GRADE-1);
  const double coeff = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM)(0);
  // const double expEta = exp(intercept + coeff*ABILITY);
  // const double prob = expEta/(1+expEta);

  double out = -log1pexp(-intercept - coeff*ABILITY);
  if(!LOGFLAG) out = exp(out);

  return(out);
}

//' Evaluate the probability of getting a specific grade
//'
//' @param GRADE Grade used as reference
//' @param EXAM Exam of interest
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY Ability value.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the probability of obtaining the grade `GRADE` on exam `EXAM`.
//' @export
// [[Rcpp::export]]
double pGrade(
    const unsigned int GRADE,
    const unsigned int EXAM,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const bool LOGFLAG = false
){
  double out;

  if(LOGFLAG){
    if(GRADE==N_GRADES){
      out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
    }else if(GRADE==0){
      out = log1mexp(-pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true));
    }else if(GRADE<N_GRADES & GRADE >0){
      out = R::logspace_sub(pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true), pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true));
    }
  }else{
    if(GRADE==N_GRADES){
      out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
    }else if(GRADE==0){
      out = 1 - pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
    }else if(GRADE<N_GRADES & GRADE >0){
      out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY)-pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
    }
  }


  return(out);
}

//' Evaluate the c.d.f or p.d.f of the last attempt to an exam
//'
//' @param EXAM Exam of interest.
//' @param DAY Day of interest.
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param SPEED speed value.
//' @param CDFFLAG `TRUE` for c.d.f. of time. `FALSE` for p.d.f.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @export
// [[Rcpp::export]]
double pTimeExam(
    const unsigned int EXAM,
    const double DAY,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double SPEED,
    const bool CDFFLAG,
    const bool LOGFLAG = false
){
  std::vector<double> pars(2);
  pars[0] = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 3, EXAM)(0);
  pars[1] = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 4, EXAM)(0);
  const double mean = pars[0] - SPEED;
  const double sd = 1/pars[1];
  double out;
  if(CDFFLAG){
    out = R::plnorm(DAY, mean, sd, true, LOGFLAG);
  }else{
    out = R::dlnorm(DAY, mean, sd, LOGFLAG);
  }

  return(out);
}

//' Evaluate exam specific likelihood
//'
//' @param EXAM Exam of interest.
//' @param GRADE Grade of interest.
//' @param DAY Day of interest.
//' @param OBSFLAG TRUE for observed, FALSE for not-observed.
//' @param THETA_IRT Portion of the parameter vector related to the IRT model
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams.
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the probability of observing or not a specific
//' grade on a given exam before a given day conditioned on ability and speed.
//'
//' @export
// [[Rcpp::export]]
double examLik(
    const unsigned int EXAM,
    const unsigned int GRADE,
    const double DAY,
    const bool OBSFLAG,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED,
    const bool LOGFLAG = false
){
  double out, logpExam, logpTime;

  if(LOGFLAG){
    if(OBSFLAG){
      logpExam = pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
      out = logpExam+logpTime;
    }else{
      logpExam = pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
      out = log1mexp(-logpExam-logpTime);
    }
  }else{
    if(OBSFLAG){
      logpExam = pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
      out = exp(logpExam+logpTime);
    }else{
      logpExam = pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
      out = 1 - exp(logpExam+logpTime);
    }
  }


  return(out);

}

#endif
