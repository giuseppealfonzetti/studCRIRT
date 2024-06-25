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
//' @param LOGFLAG Set TRUE to return log value.
//' @returns It returns the hazard probability of the specific outcome and year.
//'
//' @export
// [[Rcpp::export]]
double hazard(
  const unsigned int OUTCOME,
  const unsigned int YEAR,
  Eigen::VectorXd& THETA_CR,
  Eigen::VectorXd& COVARIATES,
  const unsigned int NYB,
  const unsigned int NYA,
  const bool LOGFLAG = false
){

  double out;
  if((OUTCOME == 1 | OUTCOME == 2) & YEAR > NYB) Rcpp::stop("`YEAR` larger than `NYB`");
  if(OUTCOME == 3 & YEAR > (NYB+NYA)) Rcpp::stop("`YEAR` larger than `NYB+NYA`");

  if(OUTCOME == 1 | OUTCOME == 2){
    const double int1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 4)(YEAR-1);
    const Eigen::VectorXd beta1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 1);
    const double eta1 = int1 + beta1.dot(COVARIATES);
    const double expEta1 = exp(int1 + beta1.dot(COVARIATES));

    const double int2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 5)(YEAR-1);
    const Eigen::VectorXd beta2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 2);
    const double eta2 = int2 + beta2.dot(COVARIATES);
    const double expEta2 = exp(int2 + beta2.dot(COVARIATES));

    if(OUTCOME == 1){
      out = eta1-log1pexp(R::logspace_add(eta1, eta2));
      // out = expEta1 /(1+expEta1+expEta2);
    }else if(OUTCOME == 2){
      out = eta2-log1pexp(R::logspace_add(eta1, eta2));
      // out = expEta2 /(1+expEta1+expEta2);
    }
  }

  if(OUTCOME == 3){
    // const double int3 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 6)(YEAR-1);
    // const Eigen::VectorXd beta3 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 3);
    // out = -log1pexp(- int3 - beta3.dot(COVARIATES));
    out = 0; // Probability of graduation when all exams are done is fixed to 1
  }

  if(LOGFLAG){
    return(out);
  }else{
    return(exp(out));
  }
}



//' Evaluate hazard function based on outcome and year
//'
//' @param OUTCOME 1 for dropout, 2 for transfer, 3 for graduation
//' @param YEAR Possible values 1:NYB in case of dropout/transfer; 1:NYA in case of graduation
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external covariates
//' @param NYB number of years in the non-graduatable state. Needed for determining how many time-related intercepts.
//' @param NYA number of years in the graduatable state. Needed for determining how many time-related intercepts.
//' @param LOGFLAG Set TRUE to return log value.
//' @returns It returns the hazard probability of the specific outcome and year.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gr_hazard(
    const unsigned int OUTCOME,
    const unsigned int YEAR,
    Eigen::VectorXd& THETA_CR,
    Eigen::VectorXd& COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA
){

  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_CR.size());
  double lpr;
  if((OUTCOME == 1 | OUTCOME == 2) & (YEAR > NYB)) Rcpp::stop("`YEAR` larger than `NYB`");
  if(OUTCOME == 3 & YEAR > (NYB+NYA)) Rcpp::stop("`YEAR` larger than `NYB+NYA`");

  Eigen::VectorXd gr_tmp = Eigen::VectorXd::Zero(COVARIATES.size()+1);
  if(OUTCOME == 1 | OUTCOME == 2){
    std::vector<unsigned int> idx1_i = extract_params_idx_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 4);
    std::vector<unsigned int> idx1_c = extract_params_idx_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 1);

    const double int1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 4)(YEAR-1);
    const Eigen::VectorXd beta1 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 1);
    const double eta1 = int1 + beta1.dot(COVARIATES);
    const double expEta1 = exp(int1 + beta1.dot(COVARIATES));


    std::vector<unsigned int> idx2_i = extract_params_idx_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 5);
    std::vector<unsigned int> idx2_c = extract_params_idx_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 2);

    const double int2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 5)(YEAR-1);
    const Eigen::VectorXd beta2 = extract_params_cr(THETA_CR, COVARIATES.size()-2, NYB, NYA, 2);
    const double eta2 = int2 + beta2.dot(COVARIATES);
    const double expEta2 = exp(int2 + beta2.dot(COVARIATES));

    double logConstNorm = log1pexp(R::logspace_add(eta1, eta2));
    double lpr1 = eta1-logConstNorm;
    double lpr2 = eta2-logConstNorm;
    double nprod = -exp(lpr1+lpr2);

    if(OUTCOME == 1){
      double tmp = exp(lpr1) - pow(exp(lpr1), 2);
      gr(idx1_i[0]+(YEAR-1)) = tmp;
      gr.segment(idx1_c[0], idx1_c[1]) = tmp*COVARIATES;
      gr.segment(idx2_c[0], idx2_c[1]) = nprod*COVARIATES;
      gr(idx2_i[0]+(YEAR-1)) = nprod;
      // out = expEta1 /(1+expEta1+expEta2);
    }else if(OUTCOME == 2){
      double tmp = exp(lpr2) - pow(exp(lpr2), 2);
      gr(idx2_i[0]+(YEAR-1)) = tmp;
      gr.segment(idx2_c[0], idx2_c[1]) = tmp*COVARIATES;
      gr.segment(idx1_c[0], idx1_c[1]) = nprod*COVARIATES;
      gr(idx1_i[0]+(YEAR-1)) = nprod;

      // out = expEta2 /(1+expEta1+expEta2);
    }
  }

  if(OUTCOME == 3){
    lpr = 0; // Probability of graduation when all exams are done is fixed to 1
  }


   return gr;
 }
//' Evaluate survival function given a the range of years of interest
//'
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the probability of survival from `YEAR FIRST` to `YEAR_LAST` included.
//'
//' @export
// [[Rcpp::export]]
double survival(
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::VectorXd& THETA_CR,
    Eigen::VectorXd& COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM = 100,
    const bool LOGFLAG = false
){

  double logout = 0;
  // double out = 1;
  if(YEAR_LAST_EXAM > YEAR_LAST){

    // Regime where graduation is not possible
    if(YEAR_LAST > NYB) Rcpp::stop("Regime 1 mismatch: `YEAR_LAST` > `NYB`");

    for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
      logout += log1mexp(-R::logspace_add(hazard(1, year, THETA_CR, COVARIATES, NYB, NYA, true), hazard(2, year, THETA_CR, COVARIATES, NYB, NYA, true)));
      // logout += log(1 - hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - hazard(2, year, THETA_CR, COVARIATES, NYB, NYA));
    }

  }else if(YEAR_LAST_EXAM <= YEAR_FIRST){

    // Regime where only graduation is possible
    if((YEAR_LAST - YEAR_LAST_EXAM + 1) > NYA) Rcpp::stop("Regime 2 mismatch: `YEAR_LAST` > `YEAR_LAST_EXAM` + `NYA` - 1");

    for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
      logout += log1mexp(-hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA, true));
      // logout +=  log(1 - hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA));
    }

  }else if(YEAR_LAST_EXAM <= YEAR_LAST & YEAR_LAST_EXAM > YEAR_FIRST){

    // Regime where graduation is possible from year `YEAR_LAST_EXAM`
    if((YEAR_LAST - YEAR_LAST_EXAM + 1) > NYA) Rcpp::stop("Mixed regime mismatch: `YEAR_LAST` > `YEAR_LAST_EXAM` + `NYA` - 1");

    for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
      logout += log1mexp(-R::logspace_add(hazard(1, year, THETA_CR, COVARIATES, NYB, NYA, true), hazard(2, year, THETA_CR, COVARIATES, NYB, NYA, true)));
      // logout += log(1 - hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - hazard(2, year, THETA_CR, COVARIATES, NYB, NYA));
    }
    for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
      logout += log1mexp(-hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA, true));
      // logout += log(1 - hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA));
    }

  }

  logout=std::max(-10000.0, logout);
  if(LOGFLAG){
    return(logout);
  }else{
    return(exp(logout));
  }
 }


//' Evaluate survival function given a the range of years of interest
//'
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @returns It returns the probability of survival from `YEAR FIRST` to `YEAR_LAST` included.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd gr_survival(
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::VectorXd& THETA_CR,
    Eigen::VectorXd& COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM = 100
){

  Eigen::VectorXd grl = Eigen::VectorXd::Zero(THETA_CR.size());
  double logout = 0;
  // double out = 1;
  if(YEAR_LAST_EXAM > YEAR_LAST){
   // Regime where graduation is not possible
   if(YEAR_LAST > NYB) Rcpp::stop("Regime 1 mismatch: `YEAR_LAST` > `NYB`");
     for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
       double logouty = log1mexp(-R::logspace_add(hazard(1, year, THETA_CR, COVARIATES, NYB, NYA, true), hazard(2, year, THETA_CR, COVARIATES, NYB, NYA, true)));
       logout += logouty;
       Eigen::VectorXd gry = - gr_hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - gr_hazard(2, year, THETA_CR, COVARIATES, NYB, NYA);
       grl += gry/exp(logouty);
       // logout += log(1 - hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - hazard(2, year, THETA_CR, COVARIATES, NYB, NYA));
     }
   }else if(YEAR_LAST_EXAM <= YEAR_FIRST){

     // Regime where only graduation is possible
     if((YEAR_LAST - YEAR_LAST_EXAM + 1) > NYA) Rcpp::stop("Regime 2 mismatch: `YEAR_LAST` > `YEAR_LAST_EXAM` + `NYA` - 1");

     for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
       logout += log1mexp(-hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA, true));
       // logout +=  log(1 - hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA));
     }

   }else if(YEAR_LAST_EXAM <= YEAR_LAST & YEAR_LAST_EXAM > YEAR_FIRST){

     // Regime where graduation is possible from year `YEAR_LAST_EXAM`
     if((YEAR_LAST - YEAR_LAST_EXAM + 1) > NYA) Rcpp::stop("Mixed regime mismatch: `YEAR_LAST` > `YEAR_LAST_EXAM` + `NYA` - 1");

     for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
       double logouty = log1mexp(-R::logspace_add(hazard(1, year, THETA_CR, COVARIATES, NYB, NYA, true), hazard(2, year, THETA_CR, COVARIATES, NYB, NYA, true)));
       logout += logouty;
       Eigen::VectorXd gry = - gr_hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - gr_hazard(2, year, THETA_CR, COVARIATES, NYB, NYA);
       grl += gry/exp(logouty);
       // logout += log(1 - hazard(1, year, THETA_CR, COVARIATES, NYB, NYA) - hazard(2, year, THETA_CR, COVARIATES, NYB, NYA));
     }
     for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
       logout += log1mexp(-hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA, true));
       // logout += log(1 - hazard(3, year - YEAR_LAST_EXAM + 1, THETA_CR, COVARIATES, NYB, NYA));
     }

   }

   return(grl*exp(logout));
 }

//' Evaluate Outcome Likelihood
//'
//' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @export
// [[Rcpp::export]]
double outcomeLik(
    const unsigned int OUTCOME,
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::VectorXd& THETA_CR,
    Eigen::VectorXd& COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM = 100,
    const bool LOGFLAG = false
){
  double logout;

  if(OUTCOME==0){
    logout = survival(YEAR_FIRST, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, true);
  }else if(OUTCOME==1|OUTCOME==2){
    if((YEAR_LAST_EXAM<=YEAR_LAST) & !LOGFLAG) return 0;
    if((YEAR_LAST_EXAM<=YEAR_LAST) & LOGFLAG) return -10000;//return R_NegInf;

    logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, true);
    logout += hazard(OUTCOME, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, true);
  }else if(OUTCOME==3){
    if((YEAR_LAST_EXAM>YEAR_LAST) & !LOGFLAG) return 0;
    if((YEAR_LAST_EXAM>YEAR_LAST) & LOGFLAG) return -10000;//return R_NegInf;

    logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, true);
    logout += hazard(3, YEAR_LAST-YEAR_LAST_EXAM+1, THETA_CR, COVARIATES, NYB, NYA, true);
  }

  if(LOGFLAG){
    return(logout);
  }else{
    return(exp(logout));
  }
}

//' Evaluate Outcome Likelihood
//'
//' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//' @param LOGFLAG Set TRUE to return log value.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd grl_outcomeLik(
    const unsigned int OUTCOME,
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    Eigen::VectorXd& THETA_CR,
    Eigen::VectorXd& COVARIATES,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM = 100
){
  Eigen::VectorXd grl = Eigen::VectorXd::Zero(THETA_CR.size());
  double logout;

  if(OUTCOME==0){

    logout = survival(YEAR_FIRST, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, true);
    grl = gr_survival(YEAR_FIRST, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM)/exp(logout);

  }else if(OUTCOME==1 | OUTCOME==2){
    if(YEAR_LAST_EXAM<=YEAR_LAST){
      logout = -10000; //return R_NegInf;
      } else {
        logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, true);
        grl = gr_survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM)/exp(logout);

        double tmp = hazard(OUTCOME, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, true);
        logout += tmp;
        grl += gr_hazard(OUTCOME, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA)/exp(tmp);
      }
  }else if(OUTCOME==3){
    if(YEAR_LAST_EXAM>YEAR_LAST){
      logout = -10000; //return R_NegInf;
    } else{
      logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, true);
      grl = gr_survival(YEAR_FIRST, YEAR_LAST-1, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM)/exp(logout);

      double tmp = hazard(3, YEAR_LAST-YEAR_LAST_EXAM+1, THETA_CR, COVARIATES, NYB, NYA, true);
      logout += tmp;
      grl += gr_hazard(3, YEAR_LAST-YEAR_LAST_EXAM+1, THETA_CR, COVARIATES, NYB, NYA);
    }


  }

  return(grl);
}


class CRloglik
{
private:
  Eigen::VectorXd _theta;
  Eigen::VectorXd _extcovariates;
  Eigen::VectorXd _exams_grades;
  Eigen::VectorXd _exams_days;
  Eigen::VectorXd _exams_obsflag;
  Eigen::VectorXd _exams_set;
  const unsigned int _outcome;
  const unsigned int _year;
  const unsigned int _n_grades;
  const unsigned int _n_exams;
  const unsigned int _nyb;
  const unsigned int _nya;
  const unsigned int _year_last_exam;
public:
  CRloglik(Eigen::VectorXd THETA,
           Eigen::VectorXd EXTCOVARIATES,
           Eigen::VectorXd EXAMS_GRADES,
           Eigen::VectorXd EXAMS_DAYS,
           Eigen::VectorXd EXAMS_OBSFLAG,
           Eigen::VectorXd EXAMS_SET,
           const unsigned int OUTCOME,
           const unsigned int YEAR,
           const unsigned int N_GRADES,
           const unsigned int N_EXAMS,
           const unsigned int NYB,
           const unsigned int NYA,
           const unsigned int YEAR_LAST_EXAM = 100):

  _theta(THETA), _extcovariates(EXTCOVARIATES), _exams_grades(EXAMS_GRADES), _exams_days(EXAMS_DAYS),
  _exams_obsflag(EXAMS_OBSFLAG), _exams_set(EXAMS_SET), _outcome(OUTCOME), _year(YEAR), _n_grades(N_GRADES), _n_exams(N_EXAMS),
  _nyb(NYB), _nya(NYA), _year_last_exam(YEAR_LAST_EXAM){}

  double operator()(const double ABILITY, const double SPEED)
  {

    const unsigned int dim_ext = _extcovariates.size();
    const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;
    const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;

    Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);
    Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);

    Eigen::VectorXd covariates(2+_extcovariates.size());
    covariates << ABILITY, SPEED, _extcovariates;

    double out = outcomeLik(_outcome, 1, _year, theta_cr, covariates, _nyb, _nya, _year_last_exam, true);

    return out;
  }

  Eigen::VectorXd gr_ll(const double ABILITY, const double SPEED)
  {

    const unsigned int dim_ext = _extcovariates.size();
    const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;
    const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;

    Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);
    Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);

    Eigen::VectorXd covariates(2+_extcovariates.size());
    covariates << ABILITY, SPEED, _extcovariates;

    double lat_par1 = _theta(dim_cr + dim_irt);
    double lat_par2 = _theta(dim_cr + dim_irt + 1);


    Eigen::VectorXd gr_cr = grl_outcomeLik(_outcome, 1, _year, theta_cr, covariates, _nyb, _nya, _year_last_exam);
    Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
    gr.segment(0, dim_cr) = gr_cr;
    return gr;
  }
};

class CR_MOD
{
private:
  Eigen::VectorXd _theta;
  Eigen::VectorXd _extcovariates;
  const unsigned int _outcome;
  const unsigned int _year;
  const unsigned int _n_grades;
  const unsigned int _n_exams;
  const unsigned int _nyb;
  const unsigned int _nya;
  const unsigned int _year_last_exam;
public:
  CR_MOD(Eigen::VectorXd THETA,
         Eigen::VectorXd EXTCOVARIATES,
         const unsigned int OUTCOME,
         const unsigned int YEAR,
         const unsigned int N_GRADES,
         const unsigned int N_EXAMS,
         const unsigned int NYB,
         const unsigned int NYA,
         const unsigned int YEAR_LAST_EXAM):

  _theta(THETA), _extcovariates(EXTCOVARIATES), _outcome(OUTCOME), _year(YEAR), _n_grades(N_GRADES), _n_exams(N_EXAMS),
  _nyb(NYB), _nya(NYA), _year_last_exam(YEAR_LAST_EXAM){}

  double ll(const double ABILITY, const double SPEED);

  Eigen::VectorXd grll(const double ABILITY, const double SPEED);
};

double CR_MOD::ll(const double ABILITY, const double SPEED){
  const unsigned int dim_ext = _extcovariates.size();
  const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;
  Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);

  Eigen::VectorXd covariates(2+_extcovariates.size());
  covariates << ABILITY, SPEED, _extcovariates;

  double out = outcomeLik(_outcome, 1, _year, theta_cr, covariates, _nyb, _nya, _year_last_exam, true);

  return out;
}

Eigen::VectorXd CR_MOD::grll(const double ABILITY, const double SPEED)
{

  const unsigned int dim_ext = _extcovariates.size();
  const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;

  Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);

  Eigen::VectorXd covariates(2+_extcovariates.size());
  covariates << ABILITY, SPEED, _extcovariates;


  Eigen::VectorXd gr_cr = grl_outcomeLik(_outcome, 1, _year, theta_cr, covariates, _nyb, _nya, _year_last_exam);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
  gr.segment(0, dim_cr) = gr_cr;
  return gr;
}

//' Evaluate Outcome Likelihood
//'
//' @param OUTCOME  `1` for dropout, `2` for transfer, `3` for graduation. `0` if no outcome is observed.
//' @param YEAR_FIRST First year to evaluate.
//' @param YEAR_LAST Last year to evaluate.
//' @param THETA_CR Portion of the parameter vector related to the competing risk model
//' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
//' @param NYB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
//' @param NYA Total number of years in the graduatable regime. Needed for determining how many time-related intercepts.
//' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
//' @param GRFLAG `TRUE` to compute the gradient.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List cr_conditional(
    Eigen::VectorXd THETA,
    Eigen::VectorXd EXTCOVARIATES,
    Eigen::VectorXd EXAMS_GRADES,
    Eigen::VectorXd EXAMS_DAYS,
    Eigen::VectorXd EXAMS_OBSFLAG,
    Eigen::VectorXd EXAMS_SET,
    const unsigned int OUTCOME,
    const unsigned int YEAR,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int YEAR_LAST_EXAM,
    const double ABILITY,
    const double SPEED,
    const bool GRFLAG = true
){
  CR_MOD cr_mod(THETA,
                EXTCOVARIATES,
                OUTCOME,
                YEAR,
                N_GRADES,
                N_EXAMS,
                NYB,
                NYA,
                YEAR_LAST_EXAM);

  double ll = cr_mod.ll(ABILITY, SPEED);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());
  if(GRFLAG){
    gr = cr_mod.grll(ABILITY, SPEED);
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = gr,
      Rcpp::Named("ll") = ll
    );

  return output;
}
#endif
