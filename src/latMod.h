#ifndef latMod_H
#define latMod_H
#include "extractParams.h"

const double pi = 3.1415926535897;
const double neglog2pi = -1.837877;

//' Joint density of speed and ability
//'
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param RHO Latent correlation.
//' @param SIGMA standard deviation of speed.
//' @param LOGFLAG `TRUE` for log density.
//'
//' @export
// [[Rcpp::export]]
double latent_distr(
  const double ABILITY,
  const double SPEED,
  const double RHO,
  const double SIGMA,
  const bool LOGFLAG = false
){
  double logout, out;
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd cov(2,2); cov << 1, RHO*SIGMA, RHO*SIGMA, pow(SIGMA, 2);
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


#endif
