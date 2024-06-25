#ifndef latMod_H
#define latMod_H
#include "extractParams.h"

const double pi = 3.1415926535897;
const double neglog2pi = -1.837877;

//' Joint density of speed and ability
//'
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param REPRHO atanh transformation of the latent correlation.
//' @param REPSIGMA log standard deviation of speed.
//' @param LOGFLAG `TRUE` for log density.
//'
//' @export
// [[Rcpp::export]]
double latent_distr(
  const double ABILITY,
  const double SPEED,
  const double REPRHO,
  const double REPSIGMA,
  const bool LOGFLAG = false
){
  double logout, out;
  double rho = tanh(REPRHO);
  double sig = exp(REPSIGMA);
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd cov(2,2); cov << 1, rho*sig, rho*sig, pow(sig, 2);
  double det = cov.determinant();
  Eigen::MatrixXd inv_cov = cov.inverse();
  logout = neglog2pi-0.5*log(det) - 0.5 * double(lat.transpose()*inv_cov*lat);

  if(LOGFLAG){
    out = logout;
  }else{
    out = exp(logout);
  }

  return out;
}

//' Joint density of speed and ability
//'
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param PAR1 L[2,1].
//' @param PAR2 L[2,2].
//' @param LOGFLAG `TRUE` for log density.
//'
//' @export
// [[Rcpp::export]]
double latent_distr2(
     const double ABILITY,
     const double SPEED,
     const double PAR1,
     const double PAR2,
     const bool LOGFLAG = false
 ){
   double logout, out;
   Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
   Eigen::MatrixXd L{{1,0},{PAR1, PAR2}};//(2,2); L << 1, 0, PAR1, PAR2;;
   Eigen::MatrixXd cov = L*L.transpose();
   double det = cov.determinant();
   Eigen::MatrixXd inv_cov = cov.inverse();
   logout = neglog2pi-0.5*log(det) - 0.5 * double(lat.transpose()*inv_cov*lat);

   if(LOGFLAG){
     out = logout;
   }else{
     out = exp(logout);
   }

   return out;
 }


//' Gradient of the log joint density of speed and ability
//'
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param PAR1 L[2,1].
//' @param PAR2 L[2,2].
//' @param LOGFLAG `TRUE` for log density.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd grl_latent_distr2(
    const double ABILITY,
    const double SPEED,
    const double PAR1,
    const double PAR2
){
  double logout, out;
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{PAR1, PAR2}};//(2,2); L << 1, 0, PAR1, PAR2;
  Eigen::MatrixXd iL = L.inverse();
  Eigen::MatrixXd Z1 = Eigen::MatrixXd::Zero(2,2); Z1(1,0) = 1;
  Eigen::MatrixXd Z2 = Eigen::MatrixXd::Zero(2,2); Z2(1,1) = 1;

  Eigen::MatrixXd M1 = -iL.transpose()*Z1.transpose()*iL.transpose()*iL;
  Eigen::MatrixXd M2 = -iL.transpose()*Z2.transpose()*iL.transpose()*iL;

  gr(0) = -(iL*Z1).trace() - .5*lat.transpose()*(M1 + M1.transpose())*lat;
  gr(1) = -(iL*Z2).trace() - .5*lat.transpose()*(M2 + M2.transpose())*lat;

  // Rcpp::List output =
  //   Rcpp::List::create(
  //     Rcpp::Named("gr") = gr,
  //     Rcpp::Named("M1") = M1,
  //     Rcpp::Named("M2") = M2,
  //     Rcpp::Named("L") = L,
  //     Rcpp::Named("Z1") = Z1,
  //     Rcpp::Named("Z2") = Z2
  //   );
  return(gr);
 }

class LAT_DISTR{
  private:
    Eigen::VectorXd _theta;

  public:
    LAT_DISTR(Eigen::VectorXd THETA):
    _theta(THETA){}

  double ll(const double ABILITY, const double SPEED);

  Eigen::VectorXd grll(const double ABILITY, const double SPEED);
};

double LAT_DISTR::ll(const double ABILITY, const double SPEED){
  double ll;
  Eigen::VectorXd pars = _theta.tail(2);
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{pars(0), pars(1)}};
  Eigen::MatrixXd cov = L*L.transpose();
  double det = cov.determinant();
  Eigen::MatrixXd inv_cov = cov.inverse();
  ll = neglog2pi-0.5*log(det) - 0.5 * double(lat.transpose()*inv_cov*lat);
  return(ll);
}

Eigen::VectorXd LAT_DISTR::grll(const double ABILITY, const double SPEED){
  Eigen::VectorXd pars = _theta.tail(2);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{pars(0), pars(1)}};
  Eigen::MatrixXd iL = L.inverse();
  Eigen::MatrixXd Z1 = Eigen::MatrixXd::Zero(2,2); Z1(1,0) = 1;
  Eigen::MatrixXd Z2 = Eigen::MatrixXd::Zero(2,2); Z2(1,1) = 1;

  Eigen::MatrixXd M1 = -iL.transpose()*Z1.transpose()*iL.transpose()*iL;
  Eigen::MatrixXd M2 = -iL.transpose()*Z2.transpose()*iL.transpose()*iL;

  gr(0) = -(iL*Z1).trace() - .5*lat.transpose()*(M1 + M1.transpose())*lat;
  gr(1) = -(iL*Z2).trace() - .5*lat.transpose()*(M2 + M2.transpose())*lat;
  return(gr);
}

//' Latent distribution
//'
//' Provides a wrapper for the cpp class `LAT_DISTR`, which computes
//' log density and gradient of the joint distribution of ability and speed.
//'
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param THETA parameters vector. The last two elements are used as (L[2,1], L[2,2]), where L is the lower triangular Cholesky of the latent covariance matrix.
//' @param GRFLAG `TRUE` to compute the gradient.
//'
//' @return It returns a list with:
//' \itemize{
//'   \item ll - The log-likelihood of the latent distribution.
//'   \item gr - The gradient of the log-likelihood.
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lat_distr(
    const double ABILITY,
    const double SPEED,
    Eigen::VectorXd THETA,
    const bool GRFLAG = true
){

  LAT_DISTR lat(THETA);

  double ll = lat.ll(ABILITY, SPEED);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(2);
  if(GRFLAG){
    gr = lat.grll(ABILITY, SPEED);
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = gr,
      Rcpp::Named("ll") = ll
    );

  return output;
 }

#endif
