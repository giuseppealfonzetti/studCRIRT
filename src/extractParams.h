#ifndef extractParams_H
#define extractParams_H

//' Intercepts reparameterisation
//'
//' @param X Intercepts or unconstrained parameter values for grades low-to-high.
//' @param CON2UN TRUE if going from the constrained space to the unconstrained one.
//'
//' @returns It allows to go back and forth from constrained intercepts
//' to unconstrained parameters.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd cppReparInt(const Eigen::VectorXd& X, bool CON2UN = true) {
  Eigen::VectorXd out(X.size());

  if (!CON2UN) {
    // From unconstrained to constrained
    Eigen::VectorXd reverse = X.reverse();
    Eigen::VectorXd param(X.size());
    param[0] = reverse[0];
    Eigen::VectorXd tmp = reverse.tail(reverse.size() - 1).array().exp();
    param.segment(1, X.size()-1) = tmp;
    out = param;
    for (unsigned int i = 1; i < X.size(); i++) {
      out[i] += out[i - 1];
    }
    return out.reverse();
  } else {

    // From constrained to unconstrained
    Eigen::VectorXd reverse = X.reverse();
    Eigen::VectorXd tmp(X.size()-1);
    for (unsigned int i = 0; i < (X.size()-1); i++) {
      tmp[i] = reverse[i + 1] - reverse[i];
    }
    tmp = tmp.array().log();
    out.segment(0, X.size()-1) = tmp.reverse();
    out(X.size()-1) = X(X.size()-1);

    return out;

  }

}

//' Extract parameters indexes related to the competing risk model
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
//' @returns It return a vector with two values representing the starting
//' index and the length of the segment of THETA_CR corresponding to
//' the parameter of interest
// [[Rcpp::export]]
std::vector<unsigned int> extract_params_idx_cr(
    Eigen::VectorXd THETA_CR,
    const unsigned int DIM_EXT,
    const unsigned int NYB,
    const unsigned int NYA,
    const unsigned int OPTION
){
  std::vector<unsigned int> out(2);

  switch(OPTION){
  case 1: // beta_d
    out[0] = 0;
    out[1] = DIM_EXT + 2;
    break;
  case 2: // beta_t
    // out = THETA_CR.segment(DIM_EXT + 2, DIM_EXT + 2);
    out[0] = DIM_EXT + 2;
    out[1] = DIM_EXT + 2;
    break;
  case 3: // beta_g
    // out = THETA_CR.segment(2*DIM_EXT + 4, DIM_EXT + 2);
    out[0] = 2*DIM_EXT + 4;
    out[1] = DIM_EXT + 2;
    break;
  case 4: // beta0_d
    // out = THETA_CR.segment(3*DIM_EXT + 6, NYB);
    out[0] = 3*DIM_EXT + 6;
    out[1] = NYB;
    break;
  case 5: // beta0_t
    // out = THETA_CR.segment(3*DIM_EXT + NYB + 6, NYB);
    out[0] = 3*DIM_EXT + NYB + 6;
    out[1] = NYB;
    break;
  case 6: // beta0_g
    // out = THETA_CR.segment(3*DIM_EXT + 2*NYB + 6, NYA);
    out[0] = 3*DIM_EXT + 2*NYB + 6;
    out[1] = NYA;
    break;
  }

  return out;
}


//' Extract parameters related to the IRT model
//' @param THETA_IRT parameter vector related to irt model
//' @param N_GRADES number of grades modelled.
//' @param N_EXAMS number of exams.
//' @param OPTION Select parameters of interest. `1` for exam-grades slopes,
//' `2` for exam-grade intercepts, `3` for time-mean parameter,
//' `4` for time-var parameter.
//' @param EXAM exam of interest. Possible values in `1:N_EXAMS`.
//'
//' @returns It return a vector with two values representing the starting
//' index and the length of the segment of THETA_IRT corresponding to
//' the parameter of interest.
//'
//' @export
// [[Rcpp::export]]
std::vector<unsigned int> extract_params_idx_irt(
     Eigen::VectorXd THETA_IRT,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const unsigned int OPTION,
     const unsigned int EXAM
 ){
  std::vector<unsigned int> out(2);

   switch(OPTION){
   case 1: // slope
     // out = THETA_IRT.segment(EXAM-1, 1);
     out[0] = EXAM-1;
     out[1] = 1;
     break;
   case 2: // intercept
     // out = THETA_IRT.segment(N_EXAMS + (EXAM-1)*N_GRADES, N_GRADES);
     out[0] = N_EXAMS + (EXAM-1)*N_GRADES;
     out[1] = N_GRADES;
     break;
   case 3: // time mean //speed-related
     // out = Eigen::VectorXd::Zero(2);
     // out(0) = THETA_IRT(N_EXAMS + N_EXAMS*N_GRADES + EXAM - 1);
     // out(1) = THETA_IRT(2*N_EXAMS + N_EXAMS*N_GRADES + EXAM - 1);
     out[0] = N_EXAMS + N_EXAMS*N_GRADES + EXAM - 1;
     out[1] = 1;
     break;
   case 4: // time var
     out[0] = 2*N_EXAMS + N_EXAMS*N_GRADES + EXAM - 1;
     out[1] = 1;
     break;
   }


   return(out);
 }

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
  std::vector<unsigned int> idx = extract_params_idx_cr(THETA_CR, DIM_EXT, NYB, NYA, OPTION);
  Eigen::VectorXd out = THETA_CR.segment(idx[0], idx[1]);

  return out;
}

//' Extract parameters related to the IRT model
//' @param THETA_IRT parameter vector related to irt model
//' @param N_GRADES number of grades modelled.
//' @param N_EXAMS number of exams.
//' @param OPTION Select parameters of interest. `1` for exam-grades slopes,
//' `2` for exam-grade intercepts, `3` for time-mean parameter,
//' `4` for time-var parameter.
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
  std::vector<unsigned int> idx = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, OPTION, EXAM);
  Eigen::VectorXd out = THETA_IRT.segment(idx[0], idx[1]);

  return(out);
}
#endif
