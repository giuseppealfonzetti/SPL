#include "utils.h"

#define N_NODES 5
double weights[] = {0.044333151939163, 0.294973376977114, 0.429812481900555, 0.207589505757111, 0.023291483426056};
double nodes[] = {1.122716044601626, 0.803901112364718, 0.576197647416307, 0.410913788475565, 0.291479895013730};



// // bivariate pnorm
// double utils::pmvnorm_cpp(Rcpp::NumericVector& upper, double rho, double abseps){
//     int n = 2;
//     int nu = 0;
//     int maxpts = 25000;      // default in mvtnorm: 25000
//     double releps = 0;      // default in mvtnorm: 0
//     int rnd = 1;            // Get/PutRNGstate
//     double* upper_ = upper.begin();
//     int infin[2] = {0,0};
//     double lower[2] = {0,1};
//     double delta[2] = {0,0};
//     double corr[1];
//     corr[0] = rho;
//     double error = 0;
//     double value = 0;
//     int inform = 0;

//     mvtnorm_C_mvtdst(&n, &nu, lower, upper_, infin, corr, delta, &maxpts,
//                      &abseps, &releps, &error, &value, &inform, &rnd);
//     return(value);
// }  

double utils::pbvnorm(Rcpp::NumericVector& UPPER, double R){
  double H1 = UPPER[0];
  double HK = UPPER[1];
  int NX=5;
  std::vector<double> X(NX);
  std::vector<double> W(NX);
  // data
  X[0]=.04691008;
  X[1]=.23076534;
  X[2]=.5;
  X[3]=.76923466;
  X[4]=.95308992;
  W[0]=.018854042;
  W[1]=.038088059;
  W[2]=.0452707394;
  W[3]=.038088059;
  W[4]=.018854042;
  // declarations
  double bv = 0;
  double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
  double cor_max = 0.7;
  double bv_fac1 = 0.13298076;
  double bv_fac2 = 0.053051647;
  // computation
  double h2 = HK;
  double h12 = (H1*H1+h2*h2)/2;
  double r_abs = std::abs(R);
  if (r_abs > cor_max){
    r2 = 1.0 - R*R;
    r3 = std::sqrt(r2);
    if (R<0){
      h2 = -h2;
    }
    h3 = H1*h2;
    h7 = std::exp( -h3 / 2.0);
    if ( r_abs < 1){
      h6 = std::abs(H1-h2);
      h5 = h6*h6 / 2.0;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8.0;
      ab = 3.0 - 2.0 * aa * h5;
      bv = bv_fac1*h6*ab*(1-R::pnorm(h6, 0, 1, 1, 0))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
      for (int ii=0; ii<NX; ii++){
        r1 = r3*X[ii];
        rr = r1*r1;
        r2 = std::sqrt( 1.0 - rr);
        bv += - W[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
      }
    }
    h11 = std::min(H1,h2);
    bv = bv*r3*h7 + R::pnorm(h11, 0, 1, 1, 0);
    if (R < 0){
      bv = R::pnorm(H1, 0, 1, 1, 0) - bv;
    }

  } else {
    h3=H1*h2;
    for (int ii=0; ii<NX; ii++){
      r1 = R*X[ii];
      rr2 = 1.0 - r1*r1;
      bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
    }
    bv = R::pnorm(H1, 0, 1, 1, 0)*R::pnorm(h2, 0, 1, 1, 0) + R*bv;
  }
  return bv;
}

//
double utils::compute_double_mean(double su, double so, double eta1, double eta2){
    double out = 0.0;
    double sAB2 = pow(su, 2.0) + pow(so, 2.0);
    for(int i=0; i<N_NODES; i++)
      for(int j=0; j<N_NODES; j++){
        Rcpp::NumericVector den = {sqrt(1.0 + pow(nodes[i], 2.0) * sAB2),  sqrt(1.0 + pow(nodes[j], 2.0) * sAB2)};
        Rcpp::NumericVector upper = {eta1 * nodes[i]/den[0], eta2 * nodes[j]/den[1]};
        double rho = nodes[i] * nodes[j] * pow(su, 2.0) / (den[0] * den[1]);
        out += weights[i]  * weights[j] * utils::pbvnorm(upper, rho);
      }
      return(out);
  }

//
double utils::compute_marginal_mean(double su, double so, double eta)
{
    double out = 0.0;
    double sAB2 = pow(su, 2.0) + pow(so, 2.0);
    for(int i=0; i<N_NODES; i++){
      double den = sqrt(1 + pow(nodes[i], 2.0) * sAB2);
      double upper = eta * nodes[i] / den;
      out += weights[i] * R::pnorm(upper, 0.0, 1.0, 1, 0);
    }
    return(out);
}

//
Eigen::VectorXd utils::compute_double_mean_gradient(
      double su,
      double so,
      double eta1,
      double eta2,
      Eigen::VectorXd v1,
      Eigen::VectorXd v2,
      int dimBeta)
{
    Eigen::VectorXd out =Eigen::VectorXd::Zero(dimBeta + 2);
    double sAB2 = pow(su, 2.0) + pow(so, 2.0);
    for(int i=0; i<N_NODES; i++)
      for(int j=0; j<N_NODES; j++){
        double a = nodes[i] * eta1 / sqrt(1.0 + pow(nodes[i], 2.0) * sAB2);
        double b = nodes[j] * eta2 / sqrt(1.0 + pow(nodes[j], 2.0) * sAB2);
        Rcpp::NumericVector den = {sqrt(1.0 + pow(nodes[i], 2.0) * sAB2),  sqrt(1.0 + pow(nodes[j], 2.0) * sAB2)};
        double rho = nodes[i] * nodes[j] * pow(su, 2.0) / (den[0] * den[1]);
        double rho1 = 1.0 - pow(rho, 2.0);
        out.segment(0,dimBeta) +=  weights[i]  * weights[j] * R::dnorm(a, 0.0, 1.0, 0) * R::pnorm( (b  - a * rho) / sqrt(rho1) , 0.0, 1.0, 1, 0) *  nodes[i] * v1 /den[0];
        out.segment(0,dimBeta) +=  weights[i]  * weights[j] * R::dnorm(b, 0.0, 1.0, 0) * R::pnorm( (a  - b * rho) / sqrt(rho1) , 0.0, 1.0, 1, 0) *  nodes[j] * v2 /den[1];
        out(dimBeta) -= su *  weights[i]  * weights[j] * R::dnorm(a, 0.0, 1.0, 0) *  R::pnorm( (b  - a * rho) / sqrt(rho1), 0.0, 1.0, 1, 0) * pow(nodes[i], 3.0) * eta1 / pow(den[0], 3.0);
        out(dimBeta) -= su * weights[i]  * weights[j] * R::dnorm(b, 0.0, 1.0, 0) *  R::pnorm( (a  - b * rho) / sqrt(rho1), 0.0, 1.0, 1, 0) * pow(nodes[j], 3.0) * eta2 / pow(den[1], 3.0);
        double arg = pow(a, 2.0) +  pow(b, 2.0) - 2.0 * rho * a * b;
        arg *= -0.5;
        double Den = pow(den[0] * den[1], 2.0);
        double Densu = 2.0 * su * (pow(nodes[i], 2.0) +  pow(nodes[j], 2.0) + pow(nodes[i], 2.0) * pow(nodes[j], 2.0) * 2 * sAB2);
        double drhosu = nodes[i] * nodes[j] * (2.0 * su * sqrt(Den) - pow(su,2.0) * 0.5 * (1 / sqrt(Den)) * Densu) / Den;
        out(dimBeta) += weights[i]  * weights[j] * exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1)) *  drhosu;
        out(dimBeta + 1) -=  so * weights[i]  * weights[j] * R::dnorm(a, 0.0, 1.0, 0) * R::pnorm( (b  - a * rho) / sqrt(rho1), 0.0, 1.0, 1, 0) * pow(nodes[i], 3.0) * eta1 / pow(den[0], 3.0);
        out(dimBeta + 1) -=  so * weights[i]  * weights[j] * R::dnorm(b, 0.0, 1.0, 0) * R::pnorm( (a  - b * rho) / sqrt(rho1), 0.0, 1.0, 1, 0) * pow(nodes[j], 3.0) * eta2 / pow(den[1], 3.0);
        double Denso = 2.0 * so * (pow(nodes[i], 2.0) +  pow(nodes[j], 2.0) + pow(nodes[i], 2.0) * pow(nodes[j], 2.0) * 2 * sAB2);
        double drhoso = nodes[i] * nodes[j] * pow(su,2.0) * (-0.5  * (1 / sqrt(Den)) * Denso) / Den;
        out(dimBeta + 1) += weights[i]  * weights[j] * exp(arg / rho1) / (2.0 * M_PI * sqrt(rho1)) *  drhoso;
      }
      return(out);
    
}

//
Eigen::VectorXd utils::compute_marginal_mean_gradient(
      double su,
      double so,
      double eta,
      Eigen::VectorXd v,
      int dimBeta)
{
    Eigen::VectorXd out =Eigen::VectorXd::Zero(dimBeta + 2);
    double sAB2 = pow(su, 2.0) + pow(so, 2.0);
    for(int i=0; i<N_NODES; i++){
      double a = nodes[i] * eta / sqrt(1.0 + pow(nodes[i], 2.0) * sAB2);
      out.segment(0,dimBeta) += weights[i] * R::dnorm(a, 0.0, 1.0, 0) * nodes[i] / sqrt(1.0 + pow(nodes[i], 2.0) * sAB2) * v;
      out(dimBeta) -=   weights[i] *  R::dnorm(a, 0.0, 1.0, 0) * eta * pow(nodes[i], 3.0) * su / pow(1.0 + pow(nodes[i], 2.0) * sAB2, 1.5);
      out(dimBeta + 1) -= weights[i] * R::dnorm(a, 0.0, 1.0, 0) * eta * pow(nodes[i], 3.0) * so / pow(1.0 + pow(nodes[i], 2.0) * sAB2, 1.5);
    }
    return(out);
}

Eigen::MatrixXi utils::get_dict(
      const std::vector<std::vector<int>> LIST,
      const int NPAIRS){

  Eigen::MatrixXi out(NPAIRS, 2);
  int idx=0;
  for(int i=0; i<LIST.size(); i++){
    const std::vector<int> veci = LIST.at(i);
    const int ni=veci.size();
    for(int j=0; j<ni-1; j++){
      for(int jp=j+1; jp< ni; jp++){
        out(idx, 0) = veci.at(j);
        out(idx, 1) = veci.at(jp);
        idx++;
      }
    }
  }

  return out;
}


std::vector<int> utils::get_pair_idx_nested(
  const int IDX,
  const std::vector<int> VALS){
  // Note on pairs identification:
  // For each level of the factor variable representing the
  // random effect of reference, Pairs are ordered by row
  // in a upper triangular matrix
  // where rows number represent the first index
  // and columns number the second, with row<column

  // Number of rows involved
  const int n = VALS.size();

  // Compute row index i as the minimum row such that IDX < cumulative numbe of pairs
  // That is, ceiling of the minimum solution of i*n-i(i+1)/2=IDX
  // Note that it is computed with 1-based counting
  const int i = ceil((2*n-1-sqrt(pow(2*n-1,2)-8*IDX))/2);

  //Evaluate the cumulative number of pairs at row i-1
  const int s = (i-1)*n-i*(i-1)/2;

  // Identify column
  const int j = i + IDX - s;

  // Rcpp::Rcout<<"i:"<<i<<", j:"<<j<<"\n";
  std::vector<int> out = {VALS.at(i-1), VALS.at(j-1)};

  return out;
}


std::vector<int> utils::get_pair_idx(
    const int IDX,
    const std::vector<std::vector<int>> LIST,
    const std::vector<int> CUMPAIRS){

  auto it = std::find_if(CUMPAIRS.begin(), CUMPAIRS.end(), [IDX](double val) {
    return val >= IDX;});

  const int l_idx = std::distance(CUMPAIRS.begin(), it);

  // Rcpp::Rcout<<"l_idx:"<<l_idx<<"\n";
  int prev_npairs = 0;
  if(l_idx>0) prev_npairs=CUMPAIRS.at(l_idx-1);
  // Rcpp::Rcout<<"prev_npairs:"<<prev_npairs<<"\n";

  const int IIDX = IDX-prev_npairs;
  // Rcpp::Rcout<<"IIDX:"<<IIDX<<"\n";

  std::vector<int> out = utils::get_pair_idx_nested(IIDX, LIST.at(l_idx));

  return out;
}

void utils::in_place_sample(
    std::vector<int>& VEC,
    const int K,
    const int SEED
){
  int n = VEC.size();
  std::mt19937 randomizer(SEED);
  std::uniform_int_distribution<int> sampler; 

  for (int i = 0; i < K; i++) {
      sampler.param(std::uniform_int_distribution<int>::param_type(i, n - 1));
      int j = sampler(randomizer);
      std::swap(VEC[i], VEC[j]);
  }

  // return VEC;

}

std::vector<int> utils::pool_with_replacement(
    const int N,
    const int K,
    const int SEED
){
  std::mt19937 randomizer(SEED);
  std::uniform_int_distribution<int> sampler(0, N-1); 
  std::vector<int> pool(K);
  for (int i = 0; i < K; i++) {
        int j = sampler(randomizer);
        pool[i] = j;
    }
  return pool;
}

std::vector<int> utils::pool_without_replacement(
    const int N,
    const int K,
    const int SEED
){
  std::mt19937 randomizer(SEED);
  std::uniform_int_distribution<int> sampler; 
  std::vector<int> reservoir(N) ;
  std::iota(std::begin(reservoir), std::end(reservoir), 0);
  std::vector<int> pool(K);

  for (int i = 0; i < K; i++) {
      sampler.param(std::uniform_int_distribution<int>::param_type(i, N - 1));
      int j = sampler(randomizer);
      std::swap(reservoir[i], reservoir[j]);
      pool[i] = j;
  }
  return pool;
}

// Construct dictionary of pair indices where row number is their idx
Eigen::MatrixXi cpp_get_dict(
    const std::vector<std::vector<int>> LIST,
    const int NPAIRS
){

  return utils::get_dict(LIST, NPAIRS);
}