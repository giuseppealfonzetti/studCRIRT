// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// hazard
double hazard(const unsigned int OUTCOME, const unsigned int YEAR, Eigen::VectorXd& THETA_CR, Eigen::VectorXd& COVARIATES, const unsigned int NYB, const unsigned int NYA, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_hazard(SEXP OUTCOMESEXP, SEXP YEARSEXP, SEXP THETA_CRSEXP, SEXP COVARIATESSEXP, SEXP NYBSEXP, SEXP NYASEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type OUTCOME(OUTCOMESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR(YEARSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_CR(THETA_CRSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type COVARIATES(COVARIATESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYB(NYBSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYA(NYASEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(hazard(OUTCOME, YEAR, THETA_CR, COVARIATES, NYB, NYA, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// survival
double survival(const unsigned int YEAR_FIRST, const unsigned int YEAR_LAST, Eigen::VectorXd& THETA_CR, Eigen::VectorXd& COVARIATES, const unsigned int NYB, const unsigned int NYA, const unsigned int YEAR_LAST_EXAM, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_survival(SEXP YEAR_FIRSTSEXP, SEXP YEAR_LASTSEXP, SEXP THETA_CRSEXP, SEXP COVARIATESSEXP, SEXP NYBSEXP, SEXP NYASEXP, SEXP YEAR_LAST_EXAMSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_FIRST(YEAR_FIRSTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_LAST(YEAR_LASTSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_CR(THETA_CRSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type COVARIATES(COVARIATESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYB(NYBSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYA(NYASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_LAST_EXAM(YEAR_LAST_EXAMSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(survival(YEAR_FIRST, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// outcomeLik
double outcomeLik(const unsigned int OUTCOME, const unsigned int YEAR_FIRST, const unsigned int YEAR_LAST, Eigen::VectorXd& THETA_CR, Eigen::VectorXd& COVARIATES, const unsigned int NYB, const unsigned int NYA, const unsigned int YEAR_LAST_EXAM, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_outcomeLik(SEXP OUTCOMESEXP, SEXP YEAR_FIRSTSEXP, SEXP YEAR_LASTSEXP, SEXP THETA_CRSEXP, SEXP COVARIATESSEXP, SEXP NYBSEXP, SEXP NYASEXP, SEXP YEAR_LAST_EXAMSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type OUTCOME(OUTCOMESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_FIRST(YEAR_FIRSTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_LAST(YEAR_LASTSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_CR(THETA_CRSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type COVARIATES(COVARIATESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYB(NYBSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYA(NYASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_LAST_EXAM(YEAR_LAST_EXAMSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(outcomeLik(OUTCOME, YEAR_FIRST, YEAR_LAST, THETA_CR, COVARIATES, NYB, NYA, YEAR_LAST_EXAM, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// cppReparInt
Eigen::VectorXd cppReparInt(const Eigen::VectorXd& X, bool CON2UN);
RcppExport SEXP _studCRIRT_cppReparInt(SEXP XSEXP, SEXP CON2UNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type CON2UN(CON2UNSEXP);
    rcpp_result_gen = Rcpp::wrap(cppReparInt(X, CON2UN));
    return rcpp_result_gen;
END_RCPP
}
// extract_params_idx_cr
std::vector<unsigned int> extract_params_idx_cr(Eigen::VectorXd THETA_CR, const unsigned int DIM_EXT, const unsigned int NYB, const unsigned int NYA, const unsigned int OPTION);
RcppExport SEXP _studCRIRT_extract_params_idx_cr(SEXP THETA_CRSEXP, SEXP DIM_EXTSEXP, SEXP NYBSEXP, SEXP NYASEXP, SEXP OPTIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA_CR(THETA_CRSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type DIM_EXT(DIM_EXTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYB(NYBSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYA(NYASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type OPTION(OPTIONSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_params_idx_cr(THETA_CR, DIM_EXT, NYB, NYA, OPTION));
    return rcpp_result_gen;
END_RCPP
}
// extract_params_idx_irt
std::vector<unsigned int> extract_params_idx_irt(Eigen::VectorXd THETA_IRT, const unsigned int N_GRADES, const unsigned int N_EXAMS, const unsigned int OPTION, const unsigned int EXAM);
RcppExport SEXP _studCRIRT_extract_params_idx_irt(SEXP THETA_IRTSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP OPTIONSEXP, SEXP EXAMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA_IRT(THETA_IRTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type OPTION(OPTIONSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type EXAM(EXAMSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, OPTION, EXAM));
    return rcpp_result_gen;
END_RCPP
}
// extract_params_cr
Eigen::VectorXd extract_params_cr(Eigen::VectorXd THETA_CR, const unsigned int DIM_EXT, const unsigned int NYB, const unsigned int NYA, const unsigned int OPTION);
RcppExport SEXP _studCRIRT_extract_params_cr(SEXP THETA_CRSEXP, SEXP DIM_EXTSEXP, SEXP NYBSEXP, SEXP NYASEXP, SEXP OPTIONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA_CR(THETA_CRSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type DIM_EXT(DIM_EXTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYB(NYBSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYA(NYASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type OPTION(OPTIONSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_params_cr(THETA_CR, DIM_EXT, NYB, NYA, OPTION));
    return rcpp_result_gen;
END_RCPP
}
// extract_params_irt
Eigen::VectorXd extract_params_irt(Eigen::VectorXd THETA_IRT, const unsigned int N_GRADES, const unsigned int N_EXAMS, const unsigned int OPTION, const unsigned int EXAM);
RcppExport SEXP _studCRIRT_extract_params_irt(SEXP THETA_IRTSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP OPTIONSEXP, SEXP EXAMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type THETA_IRT(THETA_IRTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type OPTION(OPTIONSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type EXAM(EXAMSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, OPTION, EXAM));
    return rcpp_result_gen;
END_RCPP
}
// pGreaterGrades
double pGreaterGrades(const unsigned int GRADE, const unsigned int EXAM, Eigen::VectorXd& THETA_IRT, const unsigned int N_GRADES, const unsigned int N_EXAMS, const double ABILITY, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_pGreaterGrades(SEXP GRADESEXP, SEXP EXAMSEXP, SEXP THETA_IRTSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP ABILITYSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type GRADE(GRADESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type EXAM(EXAMSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_IRT(THETA_IRTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const double >::type ABILITY(ABILITYSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// pGrade
double pGrade(const unsigned int GRADE, const unsigned int EXAM, Eigen::VectorXd& THETA_IRT, const unsigned int N_GRADES, const unsigned int N_EXAMS, const double ABILITY, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_pGrade(SEXP GRADESEXP, SEXP EXAMSEXP, SEXP THETA_IRTSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP ABILITYSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type GRADE(GRADESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type EXAM(EXAMSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_IRT(THETA_IRTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const double >::type ABILITY(ABILITYSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// pTimeExam
double pTimeExam(const unsigned int EXAM, const double DAY, Eigen::VectorXd& THETA_IRT, const unsigned int N_GRADES, const unsigned int N_EXAMS, const double SPEED, const bool CDFFLAG, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_pTimeExam(SEXP EXAMSEXP, SEXP DAYSEXP, SEXP THETA_IRTSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP SPEEDSEXP, SEXP CDFFLAGSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type EXAM(EXAMSEXP);
    Rcpp::traits::input_parameter< const double >::type DAY(DAYSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_IRT(THETA_IRTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const double >::type SPEED(SPEEDSEXP);
    Rcpp::traits::input_parameter< const bool >::type CDFFLAG(CDFFLAGSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(pTimeExam(EXAM, DAY, THETA_IRT, N_GRADES, N_EXAMS, SPEED, CDFFLAG, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// examLik
double examLik(const unsigned int EXAM, const unsigned int GRADE, const double DAY, const bool OBSFLAG, Eigen::VectorXd& THETA_IRT, const unsigned int N_GRADES, const unsigned int N_EXAMS, const double ABILITY, const double SPEED, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_examLik(SEXP EXAMSEXP, SEXP GRADESEXP, SEXP DAYSEXP, SEXP OBSFLAGSEXP, SEXP THETA_IRTSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP ABILITYSEXP, SEXP SPEEDSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type EXAM(EXAMSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type GRADE(GRADESEXP);
    Rcpp::traits::input_parameter< const double >::type DAY(DAYSEXP);
    Rcpp::traits::input_parameter< const bool >::type OBSFLAG(OBSFLAGSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA_IRT(THETA_IRTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const double >::type ABILITY(ABILITYSEXP);
    Rcpp::traits::input_parameter< const double >::type SPEED(SPEEDSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(examLik(EXAM, GRADE, DAY, OBSFLAG, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, SPEED, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// latent_distr
double latent_distr(const double ABILITY, const double SPEED, const double REPRHO, const double REPSIGMA, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_latent_distr(SEXP ABILITYSEXP, SEXP SPEEDSEXP, SEXP REPRHOSEXP, SEXP REPSIGMASEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type ABILITY(ABILITYSEXP);
    Rcpp::traits::input_parameter< const double >::type SPEED(SPEEDSEXP);
    Rcpp::traits::input_parameter< const double >::type REPRHO(REPRHOSEXP);
    Rcpp::traits::input_parameter< const double >::type REPSIGMA(REPSIGMASEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(latent_distr(ABILITY, SPEED, REPRHO, REPSIGMA, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// integrand
double integrand(Eigen::VectorXd& THETA, Eigen::VectorXd& EXTCOVARIATES, std::vector<unsigned int>& EXAMS_GRADES, std::vector<double>& EXAMS_DAYS, std::vector<bool>& EXAMS_OBSFLAG, const unsigned int OUTCOME, const unsigned int YEAR, const unsigned int N_GRADES, const unsigned int N_EXAMS, const unsigned int NYB, const unsigned int NYA, const double ABILITY, const double SPEED, const unsigned int YEAR_LAST_EXAM, const bool LOGFLAG);
RcppExport SEXP _studCRIRT_integrand(SEXP THETASEXP, SEXP EXTCOVARIATESSEXP, SEXP EXAMS_GRADESSEXP, SEXP EXAMS_DAYSSEXP, SEXP EXAMS_OBSFLAGSEXP, SEXP OUTCOMESEXP, SEXP YEARSEXP, SEXP N_GRADESSEXP, SEXP N_EXAMSSEXP, SEXP NYBSEXP, SEXP NYASEXP, SEXP ABILITYSEXP, SEXP SPEEDSEXP, SEXP YEAR_LAST_EXAMSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type EXTCOVARIATES(EXTCOVARIATESSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int>& >::type EXAMS_GRADES(EXAMS_GRADESSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type EXAMS_DAYS(EXAMS_DAYSSEXP);
    Rcpp::traits::input_parameter< std::vector<bool>& >::type EXAMS_OBSFLAG(EXAMS_OBSFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type OUTCOME(OUTCOMESEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR(YEARSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_GRADES(N_GRADESSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_EXAMS(N_EXAMSSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYB(NYBSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NYA(NYASEXP);
    Rcpp::traits::input_parameter< const double >::type ABILITY(ABILITYSEXP);
    Rcpp::traits::input_parameter< const double >::type SPEED(SPEEDSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type YEAR_LAST_EXAM(YEAR_LAST_EXAMSEXP);
    Rcpp::traits::input_parameter< const bool >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(integrand(THETA, EXTCOVARIATES, EXAMS_GRADES, EXAMS_DAYS, EXAMS_OBSFLAG, OUTCOME, YEAR, N_GRADES, N_EXAMS, NYB, NYA, ABILITY, SPEED, YEAR_LAST_EXAM, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_studCRIRT_hazard", (DL_FUNC) &_studCRIRT_hazard, 7},
    {"_studCRIRT_survival", (DL_FUNC) &_studCRIRT_survival, 8},
    {"_studCRIRT_outcomeLik", (DL_FUNC) &_studCRIRT_outcomeLik, 9},
    {"_studCRIRT_cppReparInt", (DL_FUNC) &_studCRIRT_cppReparInt, 2},
    {"_studCRIRT_extract_params_idx_cr", (DL_FUNC) &_studCRIRT_extract_params_idx_cr, 5},
    {"_studCRIRT_extract_params_idx_irt", (DL_FUNC) &_studCRIRT_extract_params_idx_irt, 5},
    {"_studCRIRT_extract_params_cr", (DL_FUNC) &_studCRIRT_extract_params_cr, 5},
    {"_studCRIRT_extract_params_irt", (DL_FUNC) &_studCRIRT_extract_params_irt, 5},
    {"_studCRIRT_pGreaterGrades", (DL_FUNC) &_studCRIRT_pGreaterGrades, 7},
    {"_studCRIRT_pGrade", (DL_FUNC) &_studCRIRT_pGrade, 7},
    {"_studCRIRT_pTimeExam", (DL_FUNC) &_studCRIRT_pTimeExam, 8},
    {"_studCRIRT_examLik", (DL_FUNC) &_studCRIRT_examLik, 10},
    {"_studCRIRT_latent_distr", (DL_FUNC) &_studCRIRT_latent_distr, 5},
    {"_studCRIRT_integrand", (DL_FUNC) &_studCRIRT_integrand, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_studCRIRT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
