echo "======== Compiling the forward library ======== "
cd libs
cd Forward/source
make
cd ../../../

echo "========== Compiling the TMCMC library ========="
cd libs
cd TMCMC/lib
module load openmpi/1.10.7
make
make -f with_bound.Makefile
cd ../../../

echo "========= Compiling the GTMCMC library ========="
cd libs
cd GTMCMC/lib
module load openmpi/1.10.7
make -f with_bound.Makefile
cd ../../../

echo "========= Compiling the T-RJMCMC library ========"
cd libs
cd T-RJMCMC/lib
module load openmpi/1.10.7
make
cd ../../../



