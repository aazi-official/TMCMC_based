#ifndef TMCMCHEADERSEEN
#define TMCMCHEADERSEEN
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <functional>
#include <string>
#include <random>
#include <mpi.h>
#include <assert.h>

extern double* boundary;
const double SMALL_LPRIOR = -1.0e30;
const double Pi = 3.14159265358979323846;
const double BETA_MAX = 1.0;
const int MAX_TMCMC_ITER = 1e3;
extern int nc;
extern double accThreshold;
extern int maxL, minL;
typedef std::vector<double>               RealVector;
typedef std::vector<std::vector<double> > RealMatrix;
typedef std::vector<int>                  IntVector;
typedef std::vector<std::vector<int> >    IntMatrix;

RealVector flatten(const RealMatrix& v);

void readPriorParams(double *alphas, double *betas, char *types, int ndim, std::string priorParamsFile, std::string priorTypesFile);
void genPriorSampsVarL(std::default_random_engine &generator,
  std::uniform_real_distribution<double> &u_distribution,
  RealMatrix &spls, IntVector &splDim, 
  int nsamps, double *generalPrior);
void genPriorSampsFixL(std::default_random_engine &generator,
  std::uniform_real_distribution<double> &u_distribution,
  RealMatrix &spls, IntVector &splDim, int dim,
  int nsamps, double *generalPrior);
void full_prior(int nLayers, double* generalPrior, std::vector<char> &Ptypes, 
  RealVector &PriorMin, RealVector &PriorMax);
void full_prior_noSigma(int nLayers, double* generalPrior, std::vector<char> &Ptypes, 
  RealVector &PriorMin, RealVector &PriorMax);
void outProcToFile(const RealVector spls, const int ndim, const int nspl, std::string fname);
void outProcToFile(const RealMatrix& spls, std::string fname);
void shuffle_spls(RealVector &spls, RealVector &llik, RealVector &lprior, std::default_random_engine &generator);
double birthDeath(RealVector& splCand, double* spls, const char BorD, const int Zmin, const int Zmax, 
                  double jumpStd, std::default_random_engine &generator, double lambda,
                  std::normal_distribution<double> &n_distribution);
void cholesky(RealVector &A, int n);
void sort_pair(RealMatrix& spls, IntVector& dims);
int factorial(int n);
void smp_jump(RealMatrix splCand, RealVector jumpProb, RealVector candDimLocal,
  RealVector spls_flat,  RealVector splsDims, char* BorD, double* Zrange, 
  double jumpStd, std::default_random_engine &generator, 
  std::normal_distribution<double> &n_distribution,
  IntVector UniqueDims, RealMatrix CVMATS);
/*************************************************
Will be redefined for likelihood PDF
*************************************************/
class LogDensity{
public:
  LogDensity(){};
  ~LogDensity(){};
  virtual double eval(RealVector&){return 0.0;};
  virtual double eval_fixSigma(RealVector&){return 0.0;};
};

/*************************************************
Prior PDF (Uniform or Gaussian)
*************************************************/
class LogPrior{
public:
  LogPrior(){};
  ~LogPrior(){};

  static double eval(RealVector& x, double* alphas, double* betas, char* types){

    int ndim = x.size();
    double log_prior = 0.0;

    for (int j = 0; j < ndim; j++) {
      switch (types[j]){
        case 'U':
          // std::cout << betas[j] << std::endl;
          if (x[j] < alphas[j] || x[j] > betas[j])
            return SMALL_LPRIOR;
          else
            log_prior += log(1.0/(betas[j]-alphas[j]));
          break;

        case 'G':
          log_prior += -0.5*log(2.0*Pi) - log(betas[j]) -0.5*pow((x[j]-alphas[j])/betas[j],2.0);
          break;
        case 'J':
          log_prior += log(sqrt(2/pow(x[j],2)));
          break;
      }
    }
    return log_prior;
  }

  static double eval_dim(RealVector& x, double* alphas, double* betas, char* types, double lambda){

    int ndim = x.size();
    double log_prior = 0.0;
    double log_prior_dim;
    if (lambda == 0)
      log_prior_dim = 0;
    else
      log_prior_dim = log(lambda*exp(-ndim*lambda));
    for (int j = 0; j < ndim; j++) {
      switch (types[j]){
        case 'U':
          // std::cout << betas[j] << std::endl;
          if (x[j] < alphas[j] || x[j] > betas[j])
            return SMALL_LPRIOR;
          else
            log_prior += log(1.0/(betas[j]-alphas[j]));
          break;

        case 'G':
          log_prior += -0.5*log(2.0*Pi) - log(betas[j]) -0.5*pow((x[j]-alphas[j])/betas[j],2.0);
          break;
        case 'J':
          log_prior += log(sqrt(2/pow(x[j],2)));
          break;
      }
    }
    return log_prior+log_prior_dim;
  }
};

double tmcmc(int nsamps, int fixInit, double lambda,
          int iseed, int ndim, double cv,
          int write_flag, int pid, int np,
          double* generalPrior,
          double* Zrange, double jumpStd, double* BDprobs,
          LogDensity* logLikelihood);

#endif