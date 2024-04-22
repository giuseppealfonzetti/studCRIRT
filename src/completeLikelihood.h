#ifndef completeLikelihood_H
#define completeLikelihood_H
#include "crMod.h"
#include "irtMod.h"
#include "latMod.h"
// #include <RcppNumerical.h>
//
class CompleteL: public Numer::MFunc
{
private:
  Eigen::VectorXd _theta;
  Eigen::VectorXd _extcovariates;
  std::vector<unsigned int> _exams_grades;
  std::vector<double> _exams_days;
  std::vector<bool> _exams_obsflag;
  const unsigned int _outcome;
  const unsigned int _year;
  const unsigned int _n_grades;
  const unsigned int _n_exams;
  const unsigned int _nyb;
  const unsigned int _nya;
  const unsigned int _year_last_exam;
  const bool _logflag;
public:
  CompleteL(Eigen::VectorXd& THETA,
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
                     const unsigned int YEAR_LAST_EXAM = 100,
                     const bool LOGFLAG = false):
  _theta(THETA), _extcovariates(EXTCOVARIATES), _exams_grades(EXAMS_GRADES), _exams_days(EXAMS_DAYS),
  _exams_obsflag(EXAMS_OBSFLAG), _outcome(OUTCOME), _year(YEAR), _n_grades(N_GRADES), _n_exams(N_EXAMS),
  _nyb(NYB), _nya(NYA), _year_last_exam(YEAR_LAST_EXAM), _logflag(LOGFLAG){}

  double operator()(Numer::Constvec& x)
  {
    // double completeL = integrand(_theta, _extcovariates, _exams_grades, _exams_days,
    //                              _exams_obsflag, _outcome, _year, _n_grades, _n_exams,
    //                              _nyb, _nya, x[0], x[1], _year_last_exam, _logflag
    // );

    const unsigned int dim_ext = _extcovariates.size();
    const unsigned int dim_cr = 3*(dim_ext+2) + 2*(_nyb) + _nya;
    const unsigned int dim_irt = 3*_n_exams + _n_exams*_n_grades;

    Eigen::VectorXd theta_cr = _theta.segment(0, dim_cr);
    Eigen::VectorXd theta_irt = _theta.segment(dim_cr, dim_irt);

    Eigen::VectorXd covariates(2+_extcovariates.size());
    covariates << x[0], x[1], _extcovariates;

    double reparlatcorr = _theta(dim_cr+dim_irt+1);
    double reparspeedva = _theta(dim_cr+dim_irt+2);


    double pLat = latent_distr(x[0], x[1], reparlatcorr, reparspeedva, _logflag);
    double out = pLat;



    std::vector<double> pExams(_n_exams);
    for(unsigned int exam = 1; exam <= _n_exams; exam++){
      pExams[exam-1]= examLik(exam,_exams_grades[exam-1],_exams_days[exam-1],
                   _exams_obsflag[exam-1], theta_irt, _n_grades, _n_exams,
                   x[0], x[1], _logflag);
      if(_logflag){
        out+=pExams[exam-1];
      }else{
        out*=pExams[exam-1];
      }
    }

    double pCR = outcomeLik(_outcome, 1, _year, theta_cr, covariates, _nyb, _nya, _year_last_exam, _logflag);
    if(_logflag){
      out+=pCR;
    }else{
      out*=pCR;
    }
    return out;
  }
};

#endif
