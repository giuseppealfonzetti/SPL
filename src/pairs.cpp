#include "pairs.h"

double pairs::probit::llikSinglePair(
    const int Y1,
        const int Y2,
        const double ETA1,
        const double ETA2,
        const double RHO
    ){

      double out;
      Rcpp::NumericVector upper = {ETA1, ETA2};
      double I12 = utils::pbvnorm(upper, RHO);
      if(Y1 == 1 & Y2 == 1) out = log(I12);
      if(Y1 == 1 & Y2 == 0) out = log(R::pnorm(ETA1, 0.0, 1.0, 1, 0) - I12);
      if(Y1 == 0 & Y2 == 1) out = log(R::pnorm(ETA2, 0.0, 1.0, 1, 0) - I12);
      if(Y1 == 0 & Y2 == 0) out = log(1.0 - R::pnorm(ETA1, 0.0, 1.0, 1, 0) -  R::pnorm(ETA2, 0.0, 1.0, 1, 0) + I12);

      return(out);
}
    
Eigen::VectorXd pairs::probit::grllSinglePair(
    const int Y1,
    const int Y2,
    const double ETA1,
    const double ETA2,
    const Eigen::Ref<const Eigen::VectorXd> V1,
    const Eigen::Ref<const Eigen::VectorXd> V2,
    const double RHO,
    const int DIMBETA
){
    Eigen::VectorXd  ggamma = Eigen::VectorXd::Zero(DIMBETA);
    Eigen::VectorXd  out = Eigen::VectorXd::Zero(DIMBETA + 1);
    double grho = 0.0;

    Rcpp::NumericVector upper = {ETA1, ETA2};
    double rho1 = 1.0 - pow(RHO, 2.0);
    double I12 = utils::pbvnorm(upper, RHO);
    double I1 = R::pnorm(ETA1, 0.0, 1.0, 1, 0);
    double I2 = R::pnorm(ETA2, 0.0, 1.0, 1, 0);
    double d1 = R::dnorm(ETA1, 0.0, 1.0, 0);
    double d2 = R::dnorm(ETA2, 0.0, 1.0, 0);
    double arg12 = (ETA1 - RHO * ETA2) / sqrt(rho1);
    double arg21 = (ETA2 - RHO * ETA1) / sqrt(rho1);
    Eigen::VectorXd I12_gamma = (d1 * R::pnorm(arg21, 0.0, 1.0, 1, 0)) * V1 + (d2 * R::pnorm(arg12, 0.0, 1.0, 1, 0)) * V2;
    double arg = pow(ETA1, 2.0) +  pow(ETA2, 2.0) - 2.0 * RHO * ETA1 * ETA2;
    arg *= -0.5;
    double I12_rho = exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1)) ;
    if(Y1 == 1 & Y2 == 1) {
      ggamma =  I12_gamma / I12;
      grho = I12_rho / I12;
    }
    if(Y1 == 1 & Y2 == 0) {
      ggamma = (d1 * V1 - I12_gamma) / (I1 - I12);
      grho =- I12_rho /(I1 - I12);
    }
    if(Y1 == 0 & Y2 == 1)
    {
      ggamma = (d2 * V2 - I12_gamma) / (I2 - I12);
      grho =- I12_rho /(I2 - I12);
    }
    if(Y1 == 0 & Y2 == 0){
      ggamma = (-d1 * V1 - d2 * V2 + I12_gamma) / (1.0 - I1 - I2 + I12);
      grho = I12_rho / (1.0 - I1 - I2 + I12);
    }

    out << ggamma, grho;
    return(out);
}        

double pairs::logit::llikSinglePair(
    const int Y1,
    const int Y2,
    const double ETA1,
    const double ETA2,
    const double SU,
    const double SO
){

    double out;
    Rcpp::NumericVector upper = {ETA1, ETA2};
    double I12 = utils::compute_double_mean(SU, SO, ETA1, ETA2);

    if(Y1 == 1 & Y2 == 1) out = log(I12);
    if(Y1 == 1 & Y2 == 0) out = log(utils::compute_marginal_mean(SU, SO, ETA1) - I12);
    if(Y1 == 0 & Y2 == 1) out = log(utils::compute_marginal_mean(SU, SO, ETA2) - I12);
    if(Y1 == 0 & Y2 == 0) out = log(1.0 - utils::compute_marginal_mean(SU, SO, ETA1) - utils::compute_marginal_mean(SU, SO, ETA2) + I12);

    return(out);
}

Eigen::VectorXd pairs::logit::grllSinglePair(
    const int Y1,
    const int Y2,
    const double ETA1,
    const double ETA2,
    const Eigen::Ref<const Eigen::VectorXd> V1,
    const Eigen::Ref<const Eigen::VectorXd> V2,
    const double SU,
    const double SO,
    const int DIMBETA
){
    Eigen::VectorXd  out = Eigen::VectorXd::Zero(DIMBETA + 2);
    double I12 = utils::compute_double_mean(SU, SO, ETA1, ETA2);
    Eigen::VectorXd  I12_d = utils::compute_double_mean_gradient(SU, SO, ETA1,  ETA2, V1, V2, DIMBETA);

    if(Y1 == 1 & Y2 == 1) {
        out  +=  I12_d / I12;
    }

    if(Y1 == 1 & Y2 == 0) {
        double I1 = utils::compute_marginal_mean(SU, SO, ETA1);
        Eigen::VectorXd  I1_d = utils::compute_marginal_mean_gradient(SU, SO, ETA1, V1, DIMBETA);
        out += (I1_d - I12_d) / (I1 - I12);
    }

    if(Y1 == 0 & Y2 == 1)
    {
        double I2 = utils::compute_marginal_mean(SU, SO, ETA2);
        Eigen::VectorXd  I2_d = utils::compute_marginal_mean_gradient(SU, SO, ETA2, V2, DIMBETA);
        out += (I2_d - I12_d) / (I2 - I12);
    }

    if(Y1 == 0 & Y2 == 0){
        double I1 = utils::compute_marginal_mean(SU, SO, ETA1);
        double I2 = utils::compute_marginal_mean(SU, SO, ETA2);
        Eigen::VectorXd  I1_d = utils::compute_marginal_mean_gradient(SU, SO, ETA1, V1, DIMBETA);
        Eigen::VectorXd  I2_d = utils::compute_marginal_mean_gradient(SU, SO, ETA2, V2, DIMBETA);
        out += (- I1_d - I2_d + I12_d) / (1.0 - I1 - I2 + I12);
    }

    return(out);
}

double pairs::ordprobit::llikSinglePair(
    const int Y1,
    const int Y2,
    const double TAU1U,
    const double TAU1L,
    const double TAU2U,
    const double TAU2L,
    const double RHO
){
    
    Rcpp::NumericVector upper11 = {TAU1L, TAU2L};
    double I11 = utils::pbvnorm(upper11, RHO);
    Rcpp::NumericVector upper12 = {TAU1L, TAU2U};
    double I12 = utils::pbvnorm(upper12, RHO);
    Rcpp::NumericVector upper21 = {TAU1U, TAU2L};
    double I21 = utils::pbvnorm(upper21, RHO);
    Rcpp::NumericVector upper22 = {TAU1U, TAU2U};
    double I22 = utils::pbvnorm(upper22, RHO);
    double out = log(I22 + I11 - I21 - I12); 

    return out;
}

Eigen::VectorXd pairs::ordprobit::grllSinglePair(
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
){

    int K = DIMTAU + 1;
    Eigen::VectorXd out = Eigen::VectorXd::Zero(DIMBETA+DIMTAU+1);
    double rho1 = 1.0 - pow(RHO, 2.0);

    // I11
    Rcpp::NumericVector upper11 = {TAU1L, TAU2L};
    double I11 = utils::pbvnorm(upper11, RHO);
    Eigen::VectorXd dt11_1(DIMTAU);
    if(Y1 > 1) dt11_1(Y1 - 2) += 1.0;
    Eigen::VectorXd I11_tau = R::pnorm((TAU2L - RHO * TAU1L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1L, 0.0, 1.0, 0) * dt11_1; 
    Eigen::VectorXd dt11_2(DIMTAU);
    if(Y2 > 1) dt11_2(Y2 - 2) += 1.0;
    I11_tau += R::pnorm((TAU1L - RHO * TAU2L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2L, 0.0, 1.0, 0) * dt11_2;
    Eigen::VectorXd I11_beta =  -R::pnorm((TAU2L - RHO * TAU1L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1L, 0.0, 1.0, 0) * V1;
    I11_beta += -R::pnorm((TAU1L - RHO * TAU2L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2L, 0.0, 1.0, 0)  * V2;
    double arg = pow(TAU1L, 2.0) +  pow(TAU2L, 2.0) - 2.0 * RHO * TAU1L * TAU2L;
    arg *= -0.5;
    double I11_rho = exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1)) ;

    // I12
    Rcpp::NumericVector upper12 = {TAU1L, TAU2U};
    double I12 = utils::pbvnorm(upper12, RHO);
    Eigen::VectorXd dt12_1(DIMTAU);
    if(Y1 > 1) dt12_1(Y1 - 2) += 1.0;
    Eigen::VectorXd I12_tau = R::pnorm((TAU2U - RHO * TAU1L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1L, 0.0, 1.0, 0) * dt12_1; 
    Eigen::VectorXd dt12_2(DIMTAU);
    if(Y2 < K) dt12_2(Y2 - 1) += 1.0;
    I12_tau += R::pnorm((TAU1L - RHO * TAU2U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2U, 0.0, 1.0, 0) * dt12_2;
    Eigen::VectorXd I12_beta =  -R::pnorm((TAU2U - RHO * TAU1L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1L, 0.0, 1.0, 0) * V1;
    I12_beta += -R::pnorm((TAU1L - RHO * TAU2U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2U, 0.0, 1.0, 0)  * V2;
    arg = pow(TAU1L, 2.0) +  pow(TAU2U, 2.0) - 2.0 * RHO * TAU1L * TAU2U;
    arg *= -0.5;
    double I12_rho = exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1));

    // I21
    Rcpp::NumericVector upper21 = {TAU1U, TAU2L};
    double I21 = utils::pbvnorm(upper21, RHO);
    Eigen::VectorXd dt21_1(DIMTAU);
    if(Y1 < K) dt21_1(Y1 - 1) += 1.0;
    Eigen::VectorXd I21_tau = R::pnorm((TAU2L - RHO * TAU1U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1U, 0.0, 1.0, 0) * dt21_1; 
    Eigen::VectorXd dt21_2(DIMTAU);
    if(Y2 > 1) dt21_2(Y2 - 2) += 1.0;
    I21_tau += R::pnorm((TAU1U - RHO * TAU2L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2L, 0.0, 1.0, 0) * dt21_2;
    Eigen::VectorXd I21_beta =  -R::pnorm((TAU2L - RHO * TAU1U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1U, 0.0, 1.0, 0) * V1;
    I21_beta += -R::pnorm((TAU1U - RHO * TAU2L) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2L, 0.0, 1.0, 0)  * V2;
    arg = pow(TAU1U, 2.0) +  pow(TAU2L, 2.0) - 2.0 * RHO * TAU1U * TAU2L;
    arg *= -0.5;
    double I21_rho = exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1)) ;

    // I22
    Rcpp::NumericVector upper22 = {TAU1U, TAU2U};
    double I22 = utils::pbvnorm(upper22, RHO);
    Eigen::VectorXd dt22_1(DIMTAU);
    if(Y1 < K) dt22_1(Y1 - 1) += 1.0;
    Eigen::VectorXd I22_tau = R::pnorm((TAU2U - RHO * TAU1U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1U, 0.0, 1.0, 0) * dt22_1; 
    Eigen::VectorXd dt22_2(DIMTAU);
    if(Y2 < K) dt22_2(Y2 - 1) += 1.0;
    I22_tau += R::pnorm((TAU1U - RHO * TAU2U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2U, 0.0, 1.0, 0) * dt22_2;
    Eigen::VectorXd I22_beta =  -R::pnorm((TAU2U - RHO * TAU1U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU1U, 0.0, 1.0, 0) * V1;
    I22_beta += -R::pnorm((TAU1U - RHO * TAU2U) / sqrt(rho1), 0.0, 1.0, 1, 0) * R::dnorm(TAU2U, 0.0, 1.0, 0)  * V2;
    arg = pow(TAU1U, 2.0) +  pow(TAU2U, 2.0) - 2.0 * RHO * TAU1U * TAU2U;
    arg *= -0.5;
    double I22_rho = exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1)) ;

    // put everything together 
    double den = I11 + I22 - I12 - I21;

    out.segment(0, DIMTAU) = (I11_tau + I22_tau - I12_tau - I21_tau) / den;
    out.segment(DIMTAU, DIMBETA) = (I11_beta + I22_beta - I12_beta - I21_beta) / den;
    out(DIMTAU+DIMBETA) = (I11_rho + I22_rho - I12_rho - I21_rho) / den;
    return(out);

}