This folder contains a paralle C++ implementation of Transitional Markov chain
Monte Carlo (TMCMC) sampling scheme for Bayesian parameter inference and model
evidence computation.

Author(s):  Moe Khalil, Sandia National Laboratories, mkhalil@sandia.gov
            Cosmin Safta, Sandia National Laboratories

Contents:
---------

C++ source code:
  /lib

C++ example:
  /app/bimodal_likelihood: A three-dimensional inference problem with a bimodal
    likelihood function
  /app/gaussian_likelihood: A three-dimensional inference problem with a Gaussian
    likelihood function
    
Include and Library folders: (do not delete!)
  /include
  /lib

Compiling the library and example:
----------------------

Set Environment for Unix / Mac OS X:

1. Edit your Bash startup file in your favorite text editor. For Linux, this is
~/.bashrc. OS X terminal runs a login shell, and so the start up file may be
~/.bashrc, ~/.bash_profile, ~/.bash_login, or ~/.profile. See the manpage for
Bash for more information about the differences between login and non-login
shells.

2. Create an environment variable to point to installed TMCMC library:

export TMCMC_PATH=/Users/mkhalil/TMCMC_Sandia (as an example)

3. Compile library (must have mpic++ compiler installed)
For the library:
- Go to $(TMCMC_PATH)/lib
- run 'make all'

4. As a toy problem:
- Go to $(TMCMC_PATH)/app/bimodal_likelihood
- run 'make'
- run driver application with 'mpirun -np 2 tmcmc_bimodal'

5. Examine results:
- samples.dat.i contains the intermediate TMCMC samples at the ith iteration
- logprior.dat.i contains the natural logarithm of prior PDF values for the intermediate TMCMC samples at the ith iteration
- loglik.dat.i contains the natural logarithm of likelihood function values for the intermediate TMCMC samples at the ith iteration
- log_evidence.dat contains the natural logarithm of model evidence