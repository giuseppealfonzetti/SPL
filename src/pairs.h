#ifndef pairs_H
#define pairs_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.h"

namespace pairs{
    namespace probit{
        // log-likelihood of single pair
        double llikSinglePair(
            const int Y1,
            const int Y2,
            const double ETA1,
            const double ETA2,
            const double RHO
        );

        // gradient of log-likelihood single pair
        Eigen::VectorXd grllSinglePair(
            const int Y1,
            const int Y2,
            const double ETA1,
            const double ETA2,
            const Eigen::Ref<const Eigen::VectorXd> V1,
            const Eigen::Ref<const Eigen::VectorXd> V2,
            const double RHO,
            const int DIMBETA
        );
    }

    namespace logit{
        //su pertains to the shared object, so is the other
        double llikSinglePair(
            const int Y1,
            const int Y2,
            const double ETA1,
            const double ETA2,
            const double SU,
            const double SO
        );

        Eigen::VectorXd grllSinglePair(
            const int Y1,
            const int Y2,
            const double ETA1,
            const double ETA2,
            const Eigen::Ref<const Eigen::VectorXd> V1,
            const Eigen::Ref<const Eigen::VectorXd> V2,
            const double SU,
            const double SO,
            const int DIMBETA
        );


    }

    namespace ordprobit{
        double llikSinglePair(
        const int Y1,
        const int Y2,
        const double TAU1U,
        const double TAU1L,
        const double TAU2U,
        const double TAU2L,
        const double RHO
        );

        Eigen::VectorXd grllSinglePair(
            const int Y1,
            const int Y2,
            const Eigen::Ref<const Eigen::VectorXd> V1,
            const Eigen::Ref<const Eigen::VectorXd> V2,
            const double TAU1U,
            const double TAU1L,
            const double TAU2U,
            const double TAU2L,
            const double RHO,
            const int DIMBETA,
            const int DIMTAU
        );
    }
}

# endif