#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <random>

#include "utils.h"
#include "pl.h"

//' Stochastic approximator
//' @export
// [[Rcpp::export]]
Rcpp::List cpp_SA2(
    Eigen::Map<Eigen::VectorXd> Y,
    Eigen::Map<Eigen::MatrixXd> X,
    const std::string LINK,
    Eigen::Map<Eigen::MatrixXi> DICT1,
    Eigen::Map<Eigen::MatrixXi> DICT2,
    Eigen::Map<Eigen::VectorXd> START,
    const double STEP0,
    const double STEP1=1,
    const double STEP2=1e-20,
    const double STEP3=.75,
    const int SCHEDULE=2,
    const int UPDATE=0,
    const int SWITCH=0,
    const double AD1=.9,
    const double AD2=.999,
    const int PAIRS_PER_ITERATION=8,
    const int BURNE=2,
    const int MAXE=4,
    const bool ISH=true,
    const int UPE=1e5,
    const int SEED=123,
    const int VERBOSE=1,
    const int NCAT = 2
){

  ///////////////
  // READ DIMS //
  ///////////////
  int nthr = 0;
  if(NCAT>2) nthr=NCAT-1;
  const int n = Y.size();
  const int d = START.size();
  const int p = d-2-nthr;
  const int pairs1 = DICT1.rows();
  const int pairs2 = DICT2.rows();
  const double prob = static_cast<double>(2*PAIRS_PER_ITERATION)/static_cast<double>(pairs1+pairs2);
  const double scale = 1/(static_cast<double>(n));
  int update_flag = UPDATE;


  /////////////////////
  // INITIALISATIONS //
  ////////////////////

  // Initialize vectors of indexes for all the pairs
  std::vector<int> pool1(pairs1) ;
  std::vector<int> pool2(pairs2) ;
  std::iota(std::begin(pool1), std::end(pool1), 0);
  std::iota(std::begin(pool2), std::end(pool2), 0);

  // Initialise vector of trajectories to output to R
  std::vector<Eigen::VectorXd> path_theta;
  std::vector<Eigen::VectorXd> path_avtheta;
  std::vector<Eigen::VectorXd> path_avtheta2;
  std::vector<int>             path_iters2;

  std::vector<int>             path_iters;
  std::vector<double>          path_nll;
  std::vector<int>             path_iters_nll;





  // double nll;
  int burnt = 1;
  Eigen::VectorXd theta   = START;
  Eigen::VectorXd avtheta = START;
  Eigen::VectorXd mtheta  = Eigen::VectorXd::Zero(theta.size());
  Eigen::VectorXd vtheta  = Eigen::VectorXd::Zero(theta.size());

  //////////////////////////////
  // STORE INITIAL QUANTITIES //
  /////////////////////////////
  Eigen::VectorXd prev_theta = theta;
  path_iters.push_back(0);
  path_theta.push_back(theta);
  path_avtheta.push_back(avtheta);


  path_iters2.push_back(0);
  path_avtheta2.push_back(avtheta);

  //////////////////
  // OPTIMISATION //
  //////////////////
  int t=1;
  int last_iter=1;
  int upe = std::min(static_cast<double>(pairs1),
                     static_cast<double>(pairs2))/static_cast<double>(PAIRS_PER_ITERATION);

  upe = std::min(static_cast<double>(upe), static_cast<double>(UPE));
  if(VERBOSE>0) Rcpp::Rcout << "Starting...\nUpdates per cycle: "<<upe <<", pairs per dimension: "<<PAIRS_PER_ITERATION<<"\n";

  double ll=1;
  for(int epoch=0; epoch < MAXE; epoch++){
    Rcpp::checkUserInterrupt();

    Eigen::VectorXd prev_epoch_theta = theta;
    Eigen::VectorXd prev_epoch_avtheta = avtheta;

    /////////////////////
    // SAMPLING SCHEME //
    /////////////////////
    // Set-up the randomizer
    if(ISH || epoch>0){
      std::mt19937 randomizer(SEED + epoch);
      std::shuffle(pool1.begin(), pool1.end(), randomizer);
      std::shuffle(pool2.begin(), pool2.end(), randomizer);
    }


    int idx_start=0;

    /////////////////////
    // EPOCH UPDATES.  //
    /////////////////////
    for(int te = 0; te < upe; te++){
        Rcpp::checkUserInterrupt();

        // Initialize empty contributions for iteration
        Eigen::VectorXd iter_grad(d);

        ///////////////////////////
        // GRADIENT COMPUTATION  //
        ///////////////////////////
        std::vector<int> iter_chosen_pairs1;
        std::vector<int> iter_chosen_pairs2;

        // Select the pairs from shuffled indices
        for(int draw = 0; draw < PAIRS_PER_ITERATION; draw++){
            int pair_index1 = pool1.at(idx_start+draw);
            // if(VERBOSE) Rcpp::Rcout << pair_index1<<",";
            int pair_index2 = pool2.at(idx_start+draw);
            iter_chosen_pairs1.push_back(pair_index1);
            iter_chosen_pairs2.push_back(pair_index2);
        }

        idx_start += PAIRS_PER_ITERATION;

        Eigen::VectorXd g1, g2;
        if(LINK=="probit"){
            g1 = pl::probit::grllPool1D(Y, X, DICT1, iter_chosen_pairs1, theta.segment(0, p), theta(p));
            g2 = pl::probit::grllPool1D(Y, X, DICT2, iter_chosen_pairs2, theta.segment(0, p), theta(p+1));
            iter_grad << (g1.segment(0, p)+g2.segment(0, p)), g1(p) , g2(p);
        }else if(LINK=="logit"){
            g1 = pl::logit::grllPool1D(Y, X, DICT1, iter_chosen_pairs1, theta.segment(0, p), theta(p), theta(p+1));
            g2 = pl::logit::grllPool1D(Y, X, DICT2, iter_chosen_pairs2, theta.segment(0, p), theta(p+1),  theta(p));
            g2.tail(2)=g2.tail(2).reverse();
            iter_grad = g1+g2;
        }else if(LINK=="ordprobit"){
            Eigen::VectorXd tau(nthr+2);
            tau << -100, theta.segment(0,nthr), 100;
            g1 = pl::ordprobit::grllPool1D(Y, X, DICT1, iter_chosen_pairs1, tau, theta.segment(nthr, p), theta(nthr+p));
            g2 = pl::ordprobit::grllPool1D(Y, X, DICT2, iter_chosen_pairs2, tau, theta.segment(nthr, p), theta(nthr+p+1));
            iter_grad << (g1.segment(0, nthr+p)+g2.segment(0, nthr+p)), g1.tail(1) , g2.tail(1);
        }else{
            Rcpp::stop("Link not supported");
        }

        iter_grad *= scale;


        /////////////
        // UPDATE  //
        /////////////
        // Stepsize scheduling
        double step_length;
        switch (SCHEDULE){
            case 0:
                step_length=STEP0;
                break;
            case 1:
                step_length=STEP0 * STEP1 * pow(t, -STEP3);
                break;
            case 2:
                step_length=STEP0*STEP1 * pow(1 + STEP2*STEP0*(t), -STEP3);
                break;
        }

        switch(update_flag){
            case 0:{
                theta += step_length * iter_grad;
                break;
            }
            case 1:{
                mtheta = AD1*mtheta+(1-AD1)*iter_grad;
                vtheta = AD2*vtheta+(1-AD2)*Eigen::VectorXd(iter_grad.array().square());
                double step_adam = step_length*(1-AD1)*pow((1-pow(AD2, t))/(1-AD2), .5);
                theta += step_adam * Eigen::VectorXd(mtheta.array()/(vtheta.array().sqrt()+1e-8));
                break;
            }
            case 2:{
                theta += step_length * Eigen::VectorXd(iter_grad.array()/(vtheta.array().sqrt()+1e-8));
                break;
            }
        }



        



      if(epoch < BURNE){
        avtheta = theta;
      }else{
        avtheta = ( (t - burnt) * avtheta + theta ) / (t - burnt + 1);
      }

      if(t%10000==0){
        path_iters2.push_back(t);
        path_avtheta2.push_back(avtheta);
        if(VERBOSE>1) Rcpp::Rcout<<"Iter: "<< t<<"\n";
        }


      t++;
    }

    // if(VERBOSE>1) Rcpp::Rcout<<"\n";



    ///////////////////
    // STORE EPOCH   //
    ///////////////////
    path_iters.push_back(t);
    path_theta.push_back(theta);
    path_avtheta.push_back(avtheta);
    // path_nll.push_back(nll);

    ////////////
    // REPORT //
    ////////////
    Eigen::VectorXd prev_epoch_thetapdiff = (theta.array()-prev_epoch_theta.array())/prev_epoch_theta.array();
    if(VERBOSE>0) Rcpp::Rcout << "End of cycle: "<<epoch<< "| mean abs theta pdiff from prev cycle: "<< prev_epoch_thetapdiff.array().abs().mean()<<"\n";


    /////////////////
    // AUTO BURNIN //
    /////////////////
    if(epoch==(BURNE-1)){
      if(VERBOSE>0) Rcpp::Rcout << "Burn-in ended\n";
      burnt=t;
      if(SWITCH) update_flag=2;
    }



    /////////////
    // EXITING //
    /////////////
    if(epoch==(MAXE-1)){
      if(VERBOSE>0) Rcpp::Rcout << "Ended\n";
      last_iter = t;
    }

  }




  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("burnt") = burnt,
      Rcpp::Named("scale") = scale,
      Rcpp::Named("pool1") = pool1,
      Rcpp::Named("pool2") = pool2,
      Rcpp::Named("path_theta") = path_theta,
      Rcpp::Named("path_avtheta") = path_avtheta,
      Rcpp::Named("path_avtheta2") = path_avtheta2,
      Rcpp::Named("path_iters2") = path_iters2,
      Rcpp::Named("path_iters") = path_iters,
      Rcpp::Named("path_nll") = path_nll,
      Rcpp::Named("path_iters_nll") = path_iters_nll,
      Rcpp::Named("last_iter") = last_iter,
      Rcpp::Named("theta") = theta,
      Rcpp::Named("avtheta") = avtheta
    );
  return(output);
}