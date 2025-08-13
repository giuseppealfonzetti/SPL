#ifndef fullPL_H
#define fullPL_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "pl.h"

// Compute log likelihood for all the pairs
// [[Rcpp::export]]
double cpp_llikFullPool2D(
    Eigen::Map<Eigen::VectorXd> Y,
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::MatrixXi> DICT1,
    Eigen::Map<Eigen::MatrixXi> DICT2,
    Eigen::Map<Eigen::VectorXd> THETA,
    const std::string LINK,
    const int NCAT=2
);

// Compute gradient of log likelihood for all the pairs
// [[Rcpp::export]]
Eigen::VectorXd cpp_grllFullPool2D(
    Eigen::Map<Eigen::VectorXd> Y,
    Eigen::Map<Eigen::MatrixXd> X,
    Eigen::Map<Eigen::MatrixXi> DICT1,
    Eigen::Map<Eigen::MatrixXi> DICT2,
    Eigen::Map<Eigen::VectorXd> THETA,
    const std::string LINK,
    const int NCAT=2
);

#endif