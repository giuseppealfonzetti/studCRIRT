#ifndef extractParams_H
#define extractParams_H

//' Extract parameters related to the competing risk model
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
