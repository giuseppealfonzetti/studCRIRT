#ifndef extractParams_H
#define extractParams_H

//' Extract parameters related to the competing risk model
//' @param THETA_CR Parameter vector related to the competing risk model
//' @param DIM_EXT Number of external covariates in the competing risk model.
//' @param NYB Number of years in the non-graduatable state. Needed for
//' determining how many time-related intercepts in the competing risk model.
//' @param NYA Number of years in the graduatable state.
//'  Needed for determining how many time-related intercepts in the competing risk model.
//' @param OPTION It selects the parameters of interest.
//' `1` for beta_d, `2` for beta_t, `3` for beta_g, `4` for beta0_d,
//' `5` for beta0_t, `6` for beta0_g.
//'
// [[Rcpp::export]]
Eigen::VectorXd extract_params_cr(
    Eigen::VectorXd THETA_CR,
    const unsigned int DIM_EXT,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int OPTION
){
  Eigen::VectorXd out;

  switch(OPTION){
    case 1: // beta_d
      out = THETA_CR.segment(0, DIM_EXT + 2);
      break;
    case 2: // beta_t
      out = THETA_CR.segment(DIM_EXT + 2, DIM_EXT + 2);
      break;
    case 3: // beta_g
      out = THETA_CR.segment(2*DIM_EXT + 4, DIM_EXT + 2);
      break;
    case 4: // beta0_d
      out = THETA_CR.segment(3*DIM_EXT + 6, NYB);
      break;
    case 5: // beta0_t
      out = THETA_CR.segment(3*DIM_EXT + NYB + 6, NYB);
      break;
    case 6: // beta0_g
      out = THETA_CR.segment(3*DIM_EXT + 2*NYB + 6, NYA);
      break;
  }

  return out;
}

//' Extract parameters related to the IRT model
//' @param THETA_IRT parameter vector related to irt model
//' @param N_GRADES number of grades modelled.
//' @param N_EXAMS number of exams.
//' @param OPTION Select parameters of interest. `1` for exam-grades slopes,
//' `2` for exam-grade intercepts, `3` for exam-speed parameters.
//' @param EXAM exam of interest. Posible values in `1:N_EXAMS`.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd extract_params_irt(
     Eigen::VectorXd THETA_IRT,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const unsigned int OPTION,
     const unsigned int EXAM
 ){
   Eigen::VectorXd out;

   switch(OPTION){
   case 1: // coeff
     out = THETA_IRT.segment(EXAM-1, 1);
     break;
   case 2: // intercept
     out = THETA_IRT.segment(N_EXAMS + (EXAM-1)*N_GRADES, N_GRADES);
     break;
   case 3: // speed-related
    out = Eigen::VectorXd::Zero(2);
    out(0) = THETA_IRT(N_EXAMS + N_EXAMS*N_GRADES + EXAM - 1);
    out(1) = THETA_IRT(2*N_EXAMS + N_EXAMS*N_GRADES + EXAM - 1);
    break;
    }


   return(out);
 }

#endif
