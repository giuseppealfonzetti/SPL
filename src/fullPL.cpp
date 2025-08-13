#include "fullPL.h"


// Compute log likelihood for all the pairs
double cpp_llikFullPool2D(
    Eigen::Map<Eigen::VectorXd> Y,
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::MatrixXi> DICT1,
    Eigen::Map<Eigen::MatrixXi> DICT2,
    Eigen::Map<Eigen::VectorXd> THETA,
    const std::string LINK,
    const int NCAT
){

  int nthr = 0;
  if(NCAT>2) nthr=NCAT-1;
  const int d = THETA.size();
  const int p = d-2-nthr;
  std::vector<int> pool1(DICT1.rows());
  std::vector<int> pool2(DICT2.rows());
  std::iota (std::begin(pool1), std::end(pool1), 0);
  std::iota (std::begin(pool2), std::end(pool2), 0);

  double ll1, ll2;
  if(LINK=="probit"){
    ll1 = pl::probit::llikPool1D(Y, X, DICT1, pool1, THETA.segment(0, p), THETA(p));
    ll2 = pl::probit::llikPool1D(Y, X, DICT2, pool2, THETA.segment(0, p), THETA(p+1));
  }else if(LINK=="logit"){
    ll1 = pl::logit::llikPool1D(Y, X, DICT1, pool1, THETA.segment(0, p), THETA(p), THETA(p+1));
    ll2 = pl::logit::llikPool1D(Y, X, DICT2, pool2, THETA.segment(0, p), THETA(p+1), THETA(p));
  }else if(LINK=="ordprobit"){

    Eigen::VectorXd tau(nthr+2);
    tau << -100, THETA.segment(0,nthr), 100;

    ll1 = pl::ordprobit::llikPool1D(Y, X, DICT1, pool1, tau, THETA.segment(nthr, p), THETA(nthr+p));
    ll2 = pl::ordprobit::llikPool1D(Y, X, DICT2, pool2, tau, THETA.segment(nthr, p), THETA(nthr+p+1));
  }else{
    Rcpp::stop("Link not supported");
  }

  return ll1+ll2;
}

// Compute gradient of log likelihood for all the pairs
Eigen::VectorXd cpp_grllFullPool2D(
    Eigen::Map<Eigen::VectorXd> Y,
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::MatrixXi> DICT1,
    Eigen::Map<Eigen::MatrixXi> DICT2,
    Eigen::Map<Eigen::VectorXd> THETA,
    const std::string LINK,
    const int NCAT
){
  int nthr = 0;
  if(NCAT>2) nthr=NCAT-1;
  const int d = THETA.size();
  const int p = d-2-nthr;
  std::vector<int> pool1(DICT1.rows());
  std::vector<int> pool2(DICT2.rows());
  std::iota (std::begin(pool1), std::end(pool1), 0);
  std::iota (std::begin(pool2), std::end(pool2), 0);

  Eigen::VectorXd g1, g2;
  Eigen::VectorXd out = Eigen::VectorXd::Zero(d);

  if(LINK=="probit"){
    g1 = pl::probit::grllPool1D(Y, X, DICT1, pool1, THETA.segment(0, p), THETA(p));
    g2 = pl::probit::grllPool1D(Y, X, DICT2, pool2, THETA.segment(0, p), THETA(p+1));
    out << (g1.segment(0, p)+g2.segment(0, p)), g1(p) , g2(p);

  }else if(LINK=="logit"){
    g1 = pl::logit::grllPool1D(Y, X, DICT1, pool1, THETA.segment(0, p), THETA(p), THETA(p+1));
    g2 = pl::logit::grllPool1D(Y, X, DICT2, pool2, THETA.segment(0, p), THETA(p+1), THETA(p));
    g2.tail(2)=g2.tail(2).reverse();
    out = g1+g2;
  }else if(LINK=="ordprobit"){
    Eigen::VectorXd tau(nthr+2);
    tau << -100, THETA.segment(0,nthr), 100;
    g1 = pl::ordprobit::grllPool1D(Y, X, DICT1, pool1, tau, THETA.segment(nthr, p), THETA(nthr+p));
    g2 = pl::ordprobit::grllPool1D(Y, X, DICT2, pool2, tau, THETA.segment(nthr, p), THETA(nthr+p+1));
    out << (g1.segment(0, nthr+p)+g2.segment(0, nthr+p)), g1.tail(1) , g2.tail(1);

  }else{
    Rcpp::stop("Link not supported");
  }

  return out;

}