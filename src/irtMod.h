#ifndef irtMod_H
#define irtMod_H
#include "extractParams.h"
#include "latMod.h"

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

    // avoid log(0)
    out = std::max(-10000.0, out);
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
//' @param MAX_DAY Last day of observation.
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
    const double MAX_DAY,
    const bool OBSFLAG,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED,
    const bool LOGFLAG = false
){
  double out, logpExam, logpTime;

  if(OBSFLAG){
    logpExam = pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
    logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
    out = logpExam+logpTime;
  }else{
    logpExam = pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
    logpTime = pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
    out = log1mexp(-logpExam-logpTime);
    out = std::max(-10000.0, out);  // avoid log(0)
  }

  if(LOGFLAG){
    return(out);
  }else{
    return(exp(out));
  }


  return(out);

}


//






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
Eigen::VectorXd gr_pGreaterGrades(
     const unsigned int GRADE,
     const unsigned int EXAM,
     Eigen::VectorXd& THETA_IRT,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const double ABILITY
){
  if(EXAM > N_EXAMS) Rcpp::stop("`EXAM` larger than `N_EXAMS`");
  if(GRADE > N_GRADES) Rcpp::stop("`GRADE` larger than `N_GRADES`");
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());

  const double intercept = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM)(GRADE-1);
  const double coeff = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM)(0);

  double pr = exp(-log1pexp(-intercept - coeff*ABILITY));
  std::vector<unsigned int> idx_i = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM);
  std::vector<unsigned int> idx_c = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM);

  double val_i = pr - pow(pr,2);
  gr(idx_c[0]) = val_i*ABILITY;

  Eigen::VectorXd tmpi = THETA_IRT.segment(idx_i[0], idx_i[1]);
  tmpi.head(tmpi.size()-1) = tmpi.head(tmpi.size()-1).array().exp();
  tmpi.tail(1) = Eigen::VectorXd::Ones(1);

  if(GRADE>1){
    tmpi.head(GRADE-1).setZero();
  }

  gr.segment(idx_i[0], idx_i[1]) = tmpi*val_i;
  // Eigen::VectorXd tmp2(idx_c[1]); tmp2.setConstant(tmp*ABILITY);
  // gr.segment(idx_c[0], idx_c[1]) = tmp2;


  return(gr);
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
Eigen::VectorXd gr_pTimeExam(
     const unsigned int EXAM,
     const double DAY,
     Eigen::VectorXd& THETA_IRT,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const double SPEED,
     const bool CDFFLAG,
     const bool LOGFLAG = false
){
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());
  std::vector<double> pars(2);
  pars[0] = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 3, EXAM)(0);
  pars[1] = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 4, EXAM)(0);
  const double mean = pars[0] - SPEED;
  const double sd = 1/pars[1];


  std::vector<unsigned int> idx_m = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, 3, EXAM);
  std::vector<unsigned int> idx_v = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, 4, EXAM);


  if(CDFFLAG){
    const double tmp = R::dnorm(pars[1]*(log(DAY)-mean), 0, 1, false);
    gr(idx_m[0]) = -tmp*pars[1];
    gr(idx_v[0]) = tmp*(log(DAY)-mean)*pars[1];
    }else{
      const double tmp = R::dlnorm(DAY, mean, sd, false);
      gr(idx_m[0]) = tmp*pow(pars[1],2)*(log(DAY)-mean);
      gr(idx_v[0]) = (tmp/pars[1] - tmp * pars[1] * pow(log(DAY)-mean, 2))*pars[1];
    }

  return(gr);
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
Eigen::VectorXd gr_pGrade(
    const unsigned int GRADE,
    const unsigned int EXAM,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY
){
  double out;
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());

  if(GRADE==N_GRADES){
       gr = gr_pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
  }else if(GRADE==0){
       gr = - gr_pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
  }else if(GRADE<N_GRADES & GRADE >0){
       gr = gr_pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY) - gr_pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
  }

  return(gr);
}

//' Evaluate exam specific likelihood
//'
//' @param EXAM Exam of interest.
//' @param GRADE Grade of interest.
//' @param DAY Day of interest.
//' @param MAX_DAY Last day of observation.
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
Eigen::VectorXd grl_examLik(
    const unsigned int EXAM,
    const unsigned int GRADE,
    const double DAY,
    const double MAX_DAY,
    const bool OBSFLAG,
    Eigen::VectorXd& THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED,
    const bool LOGFLAG = false
){
  double logp, logpExam, logpTime;
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());


  if(OBSFLAG){
     logpExam = pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, false, true);
     Eigen::VectorXd grTime = gr_pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, false, false);
     Eigen::VectorXd grExam = gr_pGrade(GRADE,EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);

     gr = grTime*exp(logpExam) + exp(logpTime)*grExam;
     logp = logpExam+logpTime;
  }else{
     logpExam = pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES,  N_EXAMS, SPEED, true, true);
     Eigen::VectorXd grExam = gr_pGreaterGrades(1,EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
     Eigen::VectorXd grTime = gr_pTimeExam(EXAM, MAX_DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, true, false);

     gr = -(grTime*exp(logpExam) + exp(logpTime)*grExam);
     logp = log1mexp(-logpExam-logpTime);
     logp = std::max(-10000.0, logp);  // avoid log(0)
   }

   gr/=exp(logp);


   return(gr);

 }


// class IRTloglik
// {
// private:
//   Eigen::VectorXd _theta;
//   Eigen::VectorXd _extcovariates;
//   Eigen::VectorXd _exams_grades;
//   Eigen::VectorXd _exams_days;
//   Eigen::VectorXd _exams_obsflag;
//   Eigen::VectorXd _exams_set;
//   const unsigned int _outcome;
//   const unsigned int _year;
//   const unsigned int _n_grades;
//   const unsigned int _n_exams;
//   const unsigned int _nyb;
//   const unsigned int _nya;
//   const unsigned int _year_last_exam;
// public:
//   IRTloglik(Eigen::VectorXd THETA,
//             Eigen::VectorXd EXTCOVARIATES,
//             Eigen::VectorXd EXAMS_GRADES,
//             Eigen::VectorXd EXAMS_DAYS,
//             Eigen::VectorXd EXAMS_OBSFLAG,
//             Eigen::VectorXd EXAMS_SET,
//             const unsigned int OUTCOME,
//             const unsigned int YEAR,
//             const unsigned int N_GRADES,
//             const unsigned int N_EXAMS,
//             const unsigned int NYB,
//             const unsigned int NYA,
//             const unsigned int YEAR_LAST_EXAM = 100):
//   _theta(THETA), _extcovariates(EXTCOVARIATES), _exams_grades(EXAMS_GRADES), _exams_days(EXAMS_DAYS),
//   _exams_obsflag(EXAMS_OBSFLAG), _exams_set(EXAMS_SET), _outcome(OUTCOME), _year(YEAR), _n_grades(N_GRADES), _n_exams(N_EXAMS),
//   _nyb(NYB), _nya(NYA), _year_last_exam(YEAR_LAST_EXAM){}
//
//   double operator()(const double ABILITY, const double SPEED)
//   {
//
//     const unsigned int dim_ext = _extcovariates.size();
//     const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;
//     const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;
//
//     Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);
//     Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);
//
//     Eigen::VectorXd covariates(2+_extcovariates.size());
//     covariates << ABILITY, SPEED, _extcovariates;
//
//     double reparlatcorr = _theta(dim_cr + dim_irt);
//     double reparspeedva = _theta(dim_cr + dim_irt + 1);
//
//
//     double pLat = latent_distr2(ABILITY, SPEED, reparlatcorr, reparspeedva, 1);
//     double out = pLat;
//
//
//
//     std::vector<double> pExams(_n_exams);
//     for(unsigned int exam = 1; exam <= _n_exams; exam++){
//       if(_exams_set[exam-1]){
//         pExams[exam-1] = examLik(exam,_exams_grades(exam-1),_exams_days(exam-1), double(_year*365),
//                                  _exams_obsflag(exam-1), theta_irt, _n_grades, _n_exams,
//                                  ABILITY, SPEED, 1);
//         out += pExams[exam-1];
//
//       }
//     }
//
//     return out;
//   }
//
//   Eigen::VectorXd gr_ll(const double ABILITY, const double SPEED)
//   {
//
//     const unsigned int dim_ext = _extcovariates.size();
//     const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;
//     const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;
//
//     Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);
//     Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);
//
//     Eigen::VectorXd covariates(2+_extcovariates.size());
//     covariates << ABILITY, SPEED, _extcovariates;
//
//     double lat_par1 = _theta(dim_cr + dim_irt);
//     double lat_par2 = _theta(dim_cr + dim_irt + 1);
//
//
//     Eigen::VectorXd gr_lat = grl_latent_distr2(ABILITY, SPEED, lat_par1, lat_par2);
//     Eigen::VectorXd gr_irt = Eigen::VectorXd::Zero(dim_irt);
//
//     for(unsigned int exam = 1; exam <= _n_exams; exam++){
//       if(_exams_set[exam-1]){
//         gr_irt += grl_examLik(exam,_exams_grades(exam-1),_exams_days(exam-1), double(_year*365),
//                               _exams_obsflag(exam-1), theta_irt, _n_grades, _n_exams,
//                               ABILITY, SPEED);
//       }
//     }
//
//     Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
//     gr.segment(dim_cr, dim_irt) = gr_irt;
//     gr.tail(2) = gr_lat;
//     return gr;
//   }
// };
//
//
// //' Evaluate the complete IRT log likelihood function
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
// //' @param LOGFLAG Set TRUE to return log value.
// //'
// //' @returns It returns the value of the integrand function,
// //' given the parameters and the data of a single observation.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::List IRT_complete_loglik(
//      Eigen::VectorXd& THETA,
//      Eigen::VectorXd& EXTCOVARIATES,
//      Eigen::VectorXd& EXAMS_GRADES,
//      Eigen::VectorXd& EXAMS_DAYS,
//      Eigen::VectorXd& EXAMS_OBSFLAG,
//      Eigen::VectorXd& EXAMS_SET,
//      const unsigned int OUTCOME,
//      const unsigned int YEAR,
//      const unsigned int N_GRADES,
//      const unsigned int N_EXAMS,
//      const unsigned int NYB,
//      const unsigned int NYA,
//      const double ABILITY,
//      const double SPEED,
//      const unsigned int YEAR_LAST_EXAM = 10,
//      const bool LOGFLAG = false
// ){
//
//    IRTloglik LL(THETA,
//                 EXTCOVARIATES,
//                 EXAMS_GRADES,
//                 EXAMS_DAYS,
//                 EXAMS_OBSFLAG,
//                 EXAMS_SET,
//                 OUTCOME,
//                 YEAR,
//                 N_GRADES,
//                 N_EXAMS,
//                 NYB,
//                 NYA,
//                 YEAR_LAST_EXAM);
//
//
//    Rcpp::List output =
//      Rcpp::List::create(
//        Rcpp::Named("gr") = LL.gr_ll(ABILITY, SPEED),
//        Rcpp::Named("ll") = LL(ABILITY, SPEED)
//      );
//    return output;
// }
//
//  //' Evaluate the marginal IRT log likelihood function
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
//  //' @param LOGFLAG Set TRUE to return log value.
//  //'
//  //' @returns It returns the value of the integrand function,
//  //' given the parameters and the data of a single observation.
//  //'
//  //' @export
//  // [[Rcpp::export]]
//  Rcpp::List IRT_GH_loglik(
//      Eigen::VectorXd& THETA,
//      Eigen::VectorXd& EXTCOVARIATES,
//      Eigen::VectorXd& EXAMS_GRADES,
//      Eigen::VectorXd& EXAMS_DAYS,
//      Eigen::VectorXd& EXAMS_OBSFLAG,
//      Eigen::VectorXd& EXAMS_SET,
//      const unsigned int OUTCOME,
//      const unsigned int YEAR,
//      Eigen::MatrixXd GRID,
//      Eigen::VectorXd WEIGHTS,
//      const unsigned int N_GRADES,
//      const unsigned int N_EXAMS,
//      const unsigned int NYB,
//      const unsigned int NYA,
//      const unsigned int YEAR_LAST_EXAM = 10,
//      const bool GRFLAG = true
//  ){
//
//    double ll;
//    Eigen::VectorXd gr_ll = Eigen::VectorXd::Zero(THETA.size());
//
//    IRTloglik LL(THETA,
//                 EXTCOVARIATES,
//                 EXAMS_GRADES,
//                 EXAMS_DAYS,
//                 EXAMS_OBSFLAG,
//                 EXAMS_SET,
//                 OUTCOME,
//                 YEAR,
//                 N_GRADES,
//                 N_EXAMS,
//                 NYB,
//                 NYA,
//                 YEAR_LAST_EXAM);
//
//    Eigen::VectorXd f(GRID.rows());
//    Eigen::MatrixXd gr(THETA.size(), GRID.rows());
//
//    for(unsigned int point = 0; point < GRID.rows(); point++){
//      f(point) = exp(LL(GRID(point, 0), GRID(point, 1)));
//      Eigen::VectorXd gr_point = f(point)*(LL.gr_ll(GRID(point, 0), GRID(point, 1)));
//      gr.col(point) = gr_point;
//    }
//
//    ll = log(f.dot(WEIGHTS));
//    if(GRFLAG){
//     gr_ll = gr*WEIGHTS/exp(ll);
//    }
//
//    Rcpp::List output =
//      Rcpp::List::create(
//        Rcpp::Named("gr") = gr_ll,
//        Rcpp::Named("ll") = ll
//    );
//    return output;
//  }
//
//
//  //' Evaluate the marginal IRT log likelihood function
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
//  //' @param LOGFLAG Set TRUE to return log value.
//  //'
//  //' @returns It returns the value of the integrand function,
//  //' given the parameters and the data of a single observation.
//  //'
//  //' @export
//  // [[Rcpp::export]]
//  double IRT_GH_sloglik(
//      Eigen::VectorXd& THETA,
//      Eigen::MatrixXd& EXTCOVARIATES,
//      Eigen::MatrixXd& EXAMS_GRADES,
//      Eigen::MatrixXd& EXAMS_DAYS,
//      Eigen::MatrixXd& EXAMS_OBSFLAG,
//      Eigen::MatrixXd& EXAMS_SET,
//      Eigen::VectorXd& OUTCOME,
//      Eigen::VectorXd& YEAR,
//      Eigen::VectorXd& YEAR_LAST_EXAM,
//      Eigen::MatrixXd GRID,
//      Eigen::VectorXd WEIGHTS,
//      const unsigned int N_GRADES,
//      const unsigned int N_EXAMS,
//      const unsigned int NYB,
//      const unsigned int NYA
//  ){
//
//
//    const unsigned int n = EXAMS_GRADES.rows();
//    const unsigned int nq = GRID.rows();
//    Eigen::MatrixXd f(n, nq);
//
//    for(unsigned int i = 0; i < n; i++){
//      IRTloglik LL(THETA,
//                   EXTCOVARIATES.row(i),
//                   EXAMS_GRADES.row(i),
//                   EXAMS_DAYS.row(i),
//                   EXAMS_OBSFLAG.row(i),
//                   EXAMS_SET.row(i),
//                   OUTCOME(i),
//                   YEAR(i),
//                   N_GRADES,
//                   N_EXAMS,
//                   NYB,
//                   NYA,
//                   YEAR_LAST_EXAM(i));
//
//      for(unsigned int point = 0; point < GRID.rows(); point++){
//        f(i, point) = exp(LL(GRID(point, 0), GRID(point, 1)));
//      }
//    }
//
//    return (f*WEIGHTS).array().log().sum();
//  }
//
// //' Evaluate the marginal IRT log likelihood function
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
// //' @param LOGFLAG Set TRUE to return log value.
// //'
// //' @returns It returns the value of the integrand function,
// //' given the parameters and the data of a single observation.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::List IRT_GH_sloglik2(
//     Eigen::VectorXd& THETA,
//     Eigen::MatrixXd& EXTCOVARIATES,
//     Eigen::MatrixXd& EXAMS_GRADES,
//     Eigen::MatrixXd& EXAMS_DAYS,
//     Eigen::MatrixXd& EXAMS_OBSFLAG,
//     Eigen::MatrixXd& EXAMS_SET,
//     Eigen::VectorXd& OUTCOME,
//     Eigen::VectorXd& YEAR,
//     Eigen::VectorXd& YEAR_LAST_EXAM,
//     Eigen::MatrixXd GRID,
//     Eigen::VectorXd WEIGHTS,
//     const unsigned int N_GRADES,
//     const unsigned int N_EXAMS,
//     const unsigned int NYB,
//     const unsigned int NYA,
//     const bool GRFLAG = true
// ){
//   double ll = 0;
//   Eigen::VectorXd gr_ll = Eigen::VectorXd::Zero(THETA.size());
//
//   const unsigned int n = EXAMS_GRADES.rows();
//   const unsigned int nq = GRID.rows();
//
//   for(unsigned int i = 0; i < n; i++){
//     IRTloglik LL(THETA,
//                  EXTCOVARIATES.row(i),
//                  EXAMS_GRADES.row(i),
//                  EXAMS_DAYS.row(i),
//                  EXAMS_OBSFLAG.row(i),
//                  EXAMS_SET.row(i),
//                  OUTCOME(i),
//                  YEAR(i),
//                  N_GRADES,
//                  N_EXAMS,
//                  NYB,
//                  NYA,
//                  YEAR_LAST_EXAM(i));
//
//     Eigen::VectorXd f(nq);
//     Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);
//
//     for(unsigned int point = 0; point < GRID.rows(); point++){
//       f(point) = exp(LL(GRID(point, 0), GRID(point, 1)));
//       if(GRFLAG){
//         Eigen::VectorXd gr_point = f(point)*(LL.gr_ll(GRID(point, 0), GRID(point, 1)));
//         gr.col(point) = gr_point;
//       }
//     }
//
//     double lli = log(f.dot(WEIGHTS));
//     ll += lli;
//     if(GRFLAG){
//       gr_ll += gr*WEIGHTS/exp(lli);
//     }
//   }
//
//
//   Rcpp::List output =
//     Rcpp::List::create(
//       Rcpp::Named("gr") = gr_ll,
//       Rcpp::Named("ll") = ll
//     );
//
//   return output;
//
// }



class IRT_MOD
{
private:
  Eigen::VectorXd _theta;
  Eigen::VectorXd _exams_grades;
  Eigen::VectorXd _exams_days;
  Eigen::VectorXd _exams_obsflag;
  Eigen::VectorXd _exams_set;
  const unsigned int _year;
  const unsigned int _n_grades;
  const unsigned int _n_exams;
  const unsigned int _nyb;
  const unsigned int _nya;
  const unsigned int _dim_ext;
public:
  IRT_MOD(Eigen::VectorXd THETA,
          Eigen::VectorXd EXAMS_GRADES,
          Eigen::VectorXd EXAMS_DAYS,
          Eigen::VectorXd EXAMS_OBSFLAG,
          Eigen::VectorXd EXAMS_SET,
          const unsigned int YEAR,
          const unsigned int N_GRADES,
          const unsigned int N_EXAMS,
          const unsigned int NYB,
          const unsigned int NYA,
          const unsigned int DIM_EXT):
  _theta(THETA), _exams_grades(EXAMS_GRADES), _exams_days(EXAMS_DAYS),
  _exams_obsflag(EXAMS_OBSFLAG), _exams_set(EXAMS_SET), _year(YEAR), _n_grades(N_GRADES), _n_exams(N_EXAMS),
  _nyb(NYB), _nya(NYA), _dim_ext(DIM_EXT){}

  double ll(const double ABILITY, const double SPEED);


  Eigen::VectorXd grll(const double ABILITY, const double SPEED);
};

double IRT_MOD::ll(const double ABILITY, const double SPEED) {

  const unsigned int dim_cr = 3*(_dim_ext+2) + 2*(_nyb) + _nya;
  const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;

  Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);

  double out = 0;

  std::vector<double> pExams(_n_exams);
  for(unsigned int exam = 1; exam <= _n_exams; exam++){
    if(_exams_set[exam-1]){
      pExams[exam-1] = examLik(exam,_exams_grades(exam-1),_exams_days(exam-1), double(_year*365),
                               _exams_obsflag(exam-1), theta_irt, _n_grades, _n_exams,
                               ABILITY, SPEED, 1);
      out += pExams[exam-1];

    }
  }

  return out;
}
Eigen::VectorXd IRT_MOD::grll(const double ABILITY, const double SPEED)
{

  const unsigned int dim_cr = 3*(_dim_ext+2) + 2*(_nyb) + _nya;
  const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;

  Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);

  double lat_par1 = _theta(dim_cr + dim_irt);
  double lat_par2 = _theta(dim_cr + dim_irt + 1);

  Eigen::VectorXd gr_irt = Eigen::VectorXd::Zero(dim_irt);

  for(unsigned int exam = 1; exam <= _n_exams; exam++){
    if(_exams_set[exam-1]){
      gr_irt += grl_examLik(exam,_exams_grades(exam-1),_exams_days(exam-1), double(_year*365),
                            _exams_obsflag(exam-1), theta_irt, _n_grades, _n_exams,
                            ABILITY, SPEED);
    }
  }

  Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
  gr.segment(dim_cr, dim_irt) = gr_irt;
  return gr;
}

//' Evaluate the complete IRT log likelihood function
//'
//' @param THETA Suitable parameter vector as provided by \link{paramsList2vec}().
//' @param EXAMS_GRADES Vector of grades.
//' @param EXAMS_DAYS Vector of times.
//' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
//' @param EXAMS_SET Vector filled with booleans.`TRUE` elements represent exams in the study plan. `FALSE` elements non-relevant ones.
//' @param YEAR Year of evaluation.
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams modelled
//' @param NYB Number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA Number of years in the graduatable state. Needed for determining how many time-related intercepts.
//' @param DIM_EXT Number of external covariates in the CR conditional model.
//' @param ABILITY Ability value.
//' @param SPEED Speed value.
//' @param GRFLAG `TRUE` to compute the gradient.
//'
//' @return It returns a list with:
//' \itemize{
//'   \item ll - The log-likelihood of the conditional IRT model.
//'   \item gr - The gradient of the log-likelihood.
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List irt_conditional(
    Eigen::VectorXd& THETA,
    Eigen::VectorXd& EXAMS_GRADES,
    Eigen::VectorXd& EXAMS_DAYS,
    Eigen::VectorXd& EXAMS_OBSFLAG,
    Eigen::VectorXd& EXAMS_SET,
    const unsigned int YEAR,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int DIM_EXT,
    const double ABILITY,
    const double SPEED,
    const bool GRFLAG = true
){
  IRT_MOD irt_mod(THETA, EXAMS_GRADES, EXAMS_DAYS, EXAMS_OBSFLAG,
                  EXAMS_SET, YEAR, N_GRADES, N_EXAMS,
                  NYB, NYA, DIM_EXT);

  double ll = irt_mod.ll(ABILITY, SPEED);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());
  if(GRFLAG){
    gr = irt_mod.grll(ABILITY, SPEED);
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = gr,
      Rcpp::Named("ll") = ll
    );

  return output;
}
#endif

