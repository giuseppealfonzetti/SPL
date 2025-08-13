#ifndef pl_H
#define pl_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "pairs.h"

namespace pl{

  namespace probit{
        // Pairwise log-likelihood for all the pairs passed by DICT
        double llikPool1D(
            const Eigen::Ref<const Eigen::VectorXd> Y,
            const Eigen::Ref<const Eigen::MatrixXd> X,
            const Eigen::Ref<const Eigen::MatrixXi> DICT,
            const std::vector<int> SEL_IDX,
            const Eigen::Ref<const Eigen::VectorXd> B,
            const double RHO
        );

        // Gradient of Pairwise log-likelihood for all the pairs passed by DICT
        Eigen::VectorXd grllPool1D(
            const Eigen::Ref<const Eigen::VectorXd> Y,
            const Eigen::Ref<const Eigen::MatrixXd> X,
            const Eigen::Ref<const Eigen::MatrixXi> DICT,
            const std::vector<int> SEL_IDX,
            const Eigen::Ref<const Eigen::VectorXd> B,
            const double RHO
        );
  }  

  namespace logit{
    // Pairwise log-likelihood for all the pairs passed by DICT
    double llikPool1D(
        const Eigen::Ref<const Eigen::VectorXd> Y,
        const Eigen::Ref<const Eigen::MatrixXd> X,
        const Eigen::Ref<const Eigen::MatrixXi> DICT,
        const std::vector<int> SEL_IDX,
        const Eigen::Ref<const Eigen::VectorXd> B,
        const double SU,
        const double SO
    );

    // Gradient of Pairwise log-likelihood for all the pairs passed by DICT
    Eigen::VectorXd grllPool1D(
        const Eigen::Ref<const Eigen::VectorXd> Y,
        const Eigen::Ref<const Eigen::MatrixXd> X,
        const Eigen::Ref<const Eigen::MatrixXi> DICT,
        const std::vector<int> SEL_IDX,
        const Eigen::Ref<const Eigen::VectorXd> B,
        const double SU,
        const double SO
    );
  }

  namespace ordprobit{
    // Pairwise log-likelihood for all the pairs passed by DICT
    double llikPool1D(
        const Eigen::Ref<const Eigen::VectorXd> Y,
        const Eigen::Ref<const Eigen::MatrixXd> X,
        const Eigen::Ref<const Eigen::MatrixXi> DICT,
        const std::vector<int> SEL_IDX,
        const Eigen::Ref<const Eigen::VectorXd> TAU,
        const Eigen::Ref<const Eigen::VectorXd> B,
        const double RHO
    );

    // Gradient of Pairwise log-likelihood for all the pairs passed by DICT
    Eigen::VectorXd grllPool1D(
        const Eigen::Ref<const Eigen::VectorXd> Y,
        const Eigen::Ref<const Eigen::MatrixXd> X,
        const Eigen::Ref<const Eigen::MatrixXi> DICT,
        const std::vector<int> SEL_IDX,
        const Eigen::Ref<const Eigen::VectorXd> TAU,
        const Eigen::Ref<const Eigen::VectorXd> B,
        const double RHO
    );
  }
}

# endif