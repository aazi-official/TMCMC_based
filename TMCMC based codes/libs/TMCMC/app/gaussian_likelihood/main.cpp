// Description: Driver application for TMCMC applied to user-defined likelihood
// and independent priors (uniform or Gaussian as specified by priorParamsFile
// and priorTypesFile)

#include <math.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "tmcmc.h"
#include <mpi.h>

typedef std::vector<double> RealVector;

/*************************************************
Define Likelihood function - 3d bimodal likelihood
To be redefined by user
Must inherit from LogDenisty class (in tmcmc.h)
and redefine eval
*************************************************/
class LogLikelihood : public LogDensity {
public:
  LogLikelihood() : LogDensity() {};

  double eval(std::vector<double>& x){

    int ndim = x.size();
    double loglik = 0.0;
    double gamma = 1.0; //standard deviation

    for (int j = 0; j < ndim; j++) {
      loglik += -0.5*log(2.0*PI) - log(gamma) - 0.5*pow(x[j]/gamma,2.0);
    }
    return loglik;
  }
};

int main (int argc, char *argv[]) {

  int pid; // MPI process ID
  int ierr;
  int np; // number of MPI processes

  ierr = MPI_Init ( &argc, &argv );

  if ( ierr != 0 )
  {
    std::cout << "\n";
    std::cout << "Fatal error!\n";
    std::cout << "MPI_Init returned nonzero ierr.\n";
    exit ( 1 );
  }

  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &np );

  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &pid );

  double chainDim = 3; // number of parameters to be inferred
  double tmcmc_cv = 0.2; // Maximum TMCMC CoV for weights of intermediate samples

  // Prior types, in single column format, one line per parameter,
  // U for uniform, G for Gaussian 
  std::string priorTypesFile = "prior_types.dat";
  // Prior hyperparameters, 2-column format, one line per parameter,
  // lower and upper bounds for uniform prior PDF,
  // mean and standard deviation for Gaussian prior PDF
  std::string priorParamsFile = "prior_params.dat";

  int write_flag = 1; // provide intermediate TMCMC artifacts
  int nsamps = 1e4; // number of TMCMC samples per process

  // random number generator seed (different for each process)
  int ran_seed = 100*pid+4;

  // model evidence, to be returned by TMCMC
  double logevid;

  // Likelihood instance
  LogLikelihood logLikelihood;

  // Run TMCMC, get evidence
  logevid = tmcmc(nsamps, ran_seed, chainDim,
  tmcmc_cv, write_flag, pid, np,
  priorParamsFile, priorTypesFile, &logLikelihood);

  // write the log-evidence to file
  if (pid == 0){
    FILE *myfile  = fopen("log_evidence.dat","w") ;
    fprintf(myfile,"%18.12e \n",logevid);
    fclose(myfile);
  }

  MPI_Finalize();

  return 0;
}