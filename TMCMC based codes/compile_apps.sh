module load openmpi/1.10.7

echo "======== Compling the TMCMC LWD app ========"
cd apps/TMCMC/lwd_source
make
cd ../../../

echo "======== Compling the GTMCMC LWD app ========"
cd apps/GTMCMC/lwd_source
make
cd ../../../

echo "======== Compling the TRJMCMC LWD app ========"
cd apps/TRJMCMC/lwd_source
make -f rjTMCMC.Makefile
cd ../../../
