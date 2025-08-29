#include "pl.h"

double pl::probit::llikPool1D(
    const Eigen::Ref<const Eigen::VectorXd> Y,
    const Eigen::Ref<const Eigen::MatrixXd> X,
    const Eigen::Ref<const Eigen::MatrixXi> DICT,
    const std::vector<int> SEL_IDX,
    const Eigen::Ref<const Eigen::VectorXd> B,
    const double RHO
){
    double out = 0;

    const int n_pairs = SEL_IDX.size();

    for(int sel_idx=0; sel_idx<n_pairs; sel_idx++){
        const int pair_idx  = SEL_IDX.at(sel_idx);
        const int indfirst  = DICT(pair_idx,0) - 1;
        const int indsecond = DICT(pair_idx,1) - 1;
        const double eta1 = X.row(indfirst)*B;
        const double eta2 = X.row(indsecond)*B;
        const int y1 = Y(indfirst);
        const int y2 = Y(indsecond);
        out += pairs::probit::llikSinglePair(y1, y2, eta1, eta2, RHO);
    }

    return(out);
}


Eigen::VectorXd pl::probit::grllPool1D(
    const Eigen::Ref<const Eigen::VectorXd> Y,
    const Eigen::Ref<const Eigen::MatrixXd> X,
    const Eigen::Ref<const Eigen::MatrixXi> DICT,
    const std::vector<int> SEL_IDX,
    const Eigen::Ref<const Eigen::VectorXd> B,
    const double RHO
){
    const int p = B.size();
    Eigen::VectorXd  out(p + 1);
    const int n_pairs = SEL_IDX.size();

    for(int sel_idx=0; sel_idx<n_pairs; sel_idx++){
    const int pair_idx  = SEL_IDX.at(sel_idx);
    const int indfirst  = DICT(pair_idx,0) - 1;
    const int indsecond = DICT(pair_idx,1) - 1;
    const double eta1 = X.row(indfirst)*B;
    const double eta2 = X.row(indsecond)*B;
    const int y1 = Y(indfirst);
    const int y2 = Y(indsecond);
    out += pairs::probit::grllSinglePair(
        y1, y2, eta1, eta2,
        X.row(indfirst),
        X.row(indsecond),
        RHO,
        p);

    }

    return(out);
}


double pl::logit::llikPool1D(
    const Eigen::Ref<const Eigen::VectorXd> Y,
    const Eigen::Ref<const Eigen::MatrixXd> X,
    const Eigen::Ref<const Eigen::MatrixXi> DICT,
    const std::vector<int> SEL_IDX,
    const Eigen::Ref<const Eigen::VectorXd> B,
    const double SU,
    const double SO
){
    double out = 0;

    const int n_pairs = SEL_IDX.size();

    for(int sel_idx=0; sel_idx<n_pairs; sel_idx++){
    const int pair_idx  = SEL_IDX.at(sel_idx);
    const int indfirst  = DICT(pair_idx,0) - 1;
    const int indsecond = DICT(pair_idx,1) - 1;
    const double eta1 = X.row(indfirst)*B;
    const double eta2 = X.row(indsecond)*B;
    const int y1 = Y(indfirst);
    const int y2 = Y(indsecond);
    out += pairs::logit::llikSinglePair(y1, y2, eta1, eta2, SU, SO);
    }

    return(out);
}

Eigen::VectorXd pl::logit::grllPool1D(
    const Eigen::Ref<const Eigen::VectorXd> Y,
    const Eigen::Ref<const Eigen::MatrixXd> X,
    const Eigen::Ref<const Eigen::MatrixXi> DICT,
    const std::vector<int> SEL_IDX,
    const Eigen::Ref<const Eigen::VectorXd> B,
    const double SU,
    const double SO
){
    const int p = B.size();
    Eigen::VectorXd  out(p + 2);
    const int n_pairs = SEL_IDX.size();

    for(int sel_idx=0; sel_idx<n_pairs; sel_idx++){
    const int pair_idx  = SEL_IDX.at(sel_idx);
    const int indfirst  = DICT(pair_idx,0) - 1;
    const int indsecond = DICT(pair_idx,1) - 1;
    const double eta1 = X.row(indfirst)*B;
    const double eta2 = X.row(indsecond)*B;
    const int y1 = Y(indfirst);
    const int y2 = Y(indsecond);
    out += pairs::logit::grllSinglePair(
        y1, y2, eta1, eta2,
        X.row(indfirst),
        X.row(indsecond),
        SU, SO,
        p);

    }

    return(out);
}

double pl::ordprobit::llikPool1D(
    const Eigen::Ref<const Eigen::VectorXd> Y,
    const Eigen::Ref<const Eigen::MatrixXd> X,
    const Eigen::Ref<const Eigen::MatrixXi> DICT,
    const std::vector<int> SEL_IDX,
    const Eigen::Ref<const Eigen::VectorXd> TAU,
    const Eigen::Ref<const Eigen::VectorXd> B,
    const double RHO
){
    double out = 0;

    const int n_pairs = SEL_IDX.size();

    for(int sel_idx=0; sel_idx<n_pairs; sel_idx++){
    const int pair_idx  = SEL_IDX.at(sel_idx);
    const int indfirst  = DICT(pair_idx,0) - 1;
    const int indsecond = DICT(pair_idx,1) - 1;
    const double eta1 = X.row(indfirst)*B;
    const double eta2 = X.row(indsecond)*B;
    const int y1 = Y(indfirst);
    const int y2 = Y(indsecond);
    double tau1u = TAU(y1) - eta1;
    double tau1l = TAU(y1 - 1) - eta1;
    double tau2u = TAU(y2) - eta2;
    double tau2l = TAU(y2-1) - eta2;
    out += pairs::ordprobit::llikSinglePair(y1, y2, tau1u, tau1l, tau2u, tau2l, RHO);
    }

    return(out);
}

// Gradient of Pairwise log-likelihood for all the pairs passed by DICT
Eigen::VectorXd pl::ordprobit::grllPool1D(
    const Eigen::Ref<const Eigen::VectorXd> Y,
    const Eigen::Ref<const Eigen::MatrixXd> X,
    const Eigen::Ref<const Eigen::MatrixXi> DICT,
    const std::vector<int> SEL_IDX,
    const Eigen::Ref<const Eigen::VectorXd> TAU,
    const Eigen::Ref<const Eigen::VectorXd> B,
    const double RHO
){
    const int p = B.size();
    const int dimtau = TAU.size()-2;
    Eigen::VectorXd  out(dimtau + p + 1);
    const int n_pairs = SEL_IDX.size();

    for(int sel_idx=0; sel_idx<n_pairs; sel_idx++){
    const int pair_idx  = SEL_IDX.at(sel_idx);
    const int indfirst  = DICT(pair_idx,0) - 1;
    const int indsecond = DICT(pair_idx,1) - 1;
    const double eta1 = X.row(indfirst)*B;
    const double eta2 = X.row(indsecond)*B;
    const int y1 = Y(indfirst);
    const int y2 = Y(indsecond);
    double tau1u = TAU(y1) - eta1;
    double tau1l = TAU(y1 - 1) - eta1;
    double tau2u = TAU(y2) - eta2;
    double tau2l = TAU(y2-1) - eta2;
    out += pairs::ordprobit::grllSinglePair(
        y1, y2,
        X.row(indfirst),
        X.row(indsecond),
        tau1u, tau1l, tau2u, tau2l,
        RHO,
        p,
        dimtau);

    }

    return(out);
}