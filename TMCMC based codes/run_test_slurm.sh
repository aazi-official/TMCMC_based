cd apps/TMCMC/lwd_test
# if you have a server with SLURM, run the command below
sbatch run.sbatch
# otherwise, run a bash script. Same for all other tests.
# bash run.sh

cd ../../GTMCMC/lwd_test
sbatch run.sbatch
# bash run.sh

cd ../../GTMCMC/lwd_test
sbatch run.sbatch
# bash run.sh