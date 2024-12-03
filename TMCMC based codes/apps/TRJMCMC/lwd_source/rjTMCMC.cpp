// Description: Driver application for TMCMC applied to the LWD inverse problems.

#include <math.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "tmcmc.h"
#include <mpi.h>
#include "Curves.h"
#include <time.h>


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
  double TVD;
  double Dip;
  double* meas;
  int numLayers;
  int ncurves;

  void read_data(char* dataFile){
    TVD = 0;
    ncurves = 72;
    meas = (double*)malloc(ncurves*sizeof(double));
    FILE* fp = fopen(dataFile,"r");
    fscanf(fp,"%d", &numLayers);
    fscanf(fp,"%lf", &Dip);
    for(int i = 0; i < ncurves; i++){
      fscanf(fp,"%lf",&meas[i]);
    }
  }
  double eval(std::vector<double>& x){

    double* curves = (double*)malloc(ncurves*sizeof(double));
    double* input = (double*)malloc(x.size()*sizeof(double));
    int nLayers = (x.size()-1)/2;
    for(int i = 0; i < x.size(); i++){
      if(i < nLayers+2)
        input[i] = pow(10,x[i]);
      else{
        input[i] = x[i];
        if(i > nLayers+2){
          // Check if negative thickness exist
          // Give very small likelihood if thickness <=0 and quit.
          if(input[i]<=0)
          	return -1e30;

          else
          	input[i] = input[i] + input[i-1];
        }
      }
    }
    forward(curves, nLayers, &input[2], &input[2], &input[nLayers+2], TVD, Dip, 2);
    forward(&curves[24], nLayers, &input[2], &input[2], &input[nLayers+2], TVD, Dip, 6);
    forward(&curves[48], nLayers, &input[2], &input[2], &input[nLayers+2], TVD, Dip, 24);

    double* diffscurve = (double*)malloc(ncurves*sizeof(double));
    double loglik = 0;
    double* tmp = (double*)malloc(ncurves*sizeof(double));
    double sigma2;
    for(int i = 0;i < ncurves; i++){
      diffscurve[i] = meas[i] - curves[i];
      if(i%2 == 0) // odd
        sigma2 = pow(input[0],2);
      else 
        sigma2 = pow(input[1],2);
      tmp[i] = pow(diffscurve[i],2)/sigma2;   
      loglik+=tmp[i];
      loglik+=log(2*3.1415926*sigma2);
    }
    if(std::isnan(loglik))
      return -1e30;
    else
      return -.5*loglik;
  }

  double eval_fixSigma(std::vector<double>& x){

    double* curves = (double*)malloc(ncurves*sizeof(double));
    double* input = (double*)malloc(x.size()*sizeof(double));
    int nLayers = (x.size()+1)/2;
    for(int i = 0; i < x.size(); i++){
      if(i < nLayers)
        input[i] = pow(10,x[i]);
      else{
        input[i] = x[i];
        if(i > nLayers){
          // Check if negative thickness exist
          // Give very small likelihood if thickness <=0 and quit.
          if(input[i]<=0)
          	return -1e30;

          else
          	input[i] = input[i] + input[i-1];
        }
      }
    }
    forward(curves, nLayers, &input[0], &input[0], &input[nLayers], TVD, Dip, 2);
    forward(&curves[24], nLayers, &input[0], &input[0], &input[nLayers], TVD, Dip, 6);
    forward(&curves[48], nLayers, &input[0], &input[0], &input[nLayers], TVD, Dip, 24);

    double* diffscurve = (double*)malloc(ncurves*sizeof(double));
    double loglik = 0;
    double* tmp = (double*)malloc(ncurves*sizeof(double));
    double sigma2;
    for(int i = 0;i < ncurves; i++){
      diffscurve[i] = meas[i] - curves[i];
      if(i%2 == 0) // odd
        sigma2 = pow(0.75,2);
      else 
        sigma2 = pow(0.125,2);
      tmp[i] = pow(diffscurve[i],2)/sigma2;   
      loglik+=tmp[i];
      loglik+=log(2*3.1415926*sigma2);
    }
    if(std::isnan(loglik))
      return -1e30;
    else
      return -.5*loglik;
  }

};




double accThreshold;
int maxL;
int minL;
int nc;
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

  int write_flag = 1; // provide intermediate TMCMC artifacts
  //int nsamps = 1e3; // number of TMCMC samples per process
  int chainDim;
  double tmcmc_cv;
  int nsamps;
  // random number generator seed (different for each process)
  srand(time(0)); 
  int ran_seed = pid+rand();
  int fixInit;
  double lambda;
  FILE* fp = fopen("algo_config.dat","r");
  fscanf(fp,"%d %lf %d %lf %d %lf %d", 
    &chainDim, &tmcmc_cv,&nsamps, &accThreshold, &fixInit, &lambda, &nc);
  fclose(fp);
  
  // Read parameter boundaries to a extern array
 
  // double rwSigma[3];
  // fp = fopen("rwSigma.dat","r");
  // for(int i = 0; i < 3; i++)
  //   fscanf(fp, "%lf", &rwSigma[i]);
  // fclose(fp);

  double jumpStd;
  double BDprob[2];
  fp = fopen("rjSettings.dat","r");
  fscanf(fp, "%lf %lf %lf %d %d", &jumpStd, &BDprob[0], &BDprob[1], &minL,&maxL);
  fclose(fp);

  double generalPrior[10];
  // for(int i = 0; i < 8; i++)
  //   fscanf(fp, "%lf", &generalPrior[i]);
  // fclose(fp);
  fp = fopen("generalPrior.dat","r");
  int i = 0;
  while(!feof(fp))
  {
    fscanf(fp, "%lf", &generalPrior[i]);
    i += 1;
  }

  double Zrange[2];
  fp = fopen("Zrange.dat","r");
  fscanf(fp, "%lf %lf", &Zrange[0], &Zrange[1]);
  fclose(fp);
  // model evidence, to be returned by TMCMC
  double logevid;

  // Likelihood instance
  LogLikelihood logLikelihood;
  logLikelihood.read_data("measurements.dat");

  // Run TMCMC, get evidence
  logevid = tmcmc(nsamps, fixInit,lambda, ran_seed, chainDim,
    tmcmc_cv, write_flag, pid, np,
    generalPrior, Zrange, jumpStd, BDprob,
    &logLikelihood);


  // write the log-evidence to file
  if (pid == 0){
    FILE *myfile  = fopen("log_evidence.dat","w") ;
    fprintf(myfile,"%18.12e \n",logevid);
    fclose(myfile);
  }

  MPI_Finalize();

  return 0;
}