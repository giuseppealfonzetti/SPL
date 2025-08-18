#ifndef utils_H
#define utils_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>

namespace utils{

  // bivariate pnorm
  //   double pmvnorm_cpp(Rcpp::NumericVector& upper, double rho, double abseps = 1e-3);

  //pbv_rcpp_pbvnorm0 (Drezner & Wesolowksy, 1990, JCSC)
  double pbvnorm(Rcpp::NumericVector& UPPER, double R);

  //
  double compute_double_mean(double su, double so, double eta1, double eta2);

  //
  double compute_marginal_mean(double su, double so, double eta);

  //
  Eigen::VectorXd compute_double_mean_gradient(
      double su,
      double so,
      double eta1,
      double eta2,
      Eigen::VectorXd v1,
      Eigen::VectorXd v2,
      int dimBeta);

  //
  Eigen::VectorXd compute_marginal_mean_gradient(
      double su,
      double so,
      double eta,
      Eigen::VectorXd v,
      int dimBeta);


   Eigen::MatrixXi get_dict(
      const std::vector<std::vector<int>> LIST,
      const int NPAIRS);


   std::vector<int> get_pair_idx_nested(
    const int IDX,
    const std::vector<int> VALS);


   std::vector<int> get_pair_idx(
    const int IDX,
    const std::vector<std::vector<int>> LIST,
    const std::vector<int> CUMPAIRS); 
    
    void in_place_sample(
        std::vector<int>& VEC,
        const int K,
        const int SEED
    );

    std::vector<int> pool_with_replacement(
        const int N,
        const int K,
        const int SEED
    );

    std::vector<int> pool_without_replacement(
        const int N,
        const int K,
        const int SEED
    );

}

// Construct dictionary of pair indices where row number is their idx
// [[Rcpp::export]]
Eigen::MatrixXi cpp_get_dict(
    const std::vector<std::vector<int>> LIST,
    const int NPAIRS);


#endif