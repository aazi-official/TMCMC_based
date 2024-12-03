#!/bin/bash

#SBATCH -J GTMCMC     # Job name
#SBATCH -o GTMCMC%j.log # Slurm out file name

#SBATCH -N 15                  # number of nodes
#SBATCH -p xfu
#SBATCH --ntasks-per-node=16    # processors per node
#SBATCH -t 120:00:00              # time in hours
#SBATCH --mem-per-cpu=2048      # number of MB per processor
#SBATCH --mail-user=your_email_address
#SBATCH --mail-type=begin        # email me when the job starts
#SBATCH --mail-type=end          # email me when the job finishes

module load openmpi/1.10.7
module load python/3.7

outDir=outDir
mkdir $outDir
mkdir $outDir/samples
mkdir $outDir/loglike
mkdir $outDir/logproposal
for i in {1..1}
do
	mkdir point$i
	cd point$i
	sed $i,$i'!d' ../measurements/layer_n_dip.dat >> measurements.dat
	sed $i,$i'!d' ../measurements/measurements.dat >> measurements.dat

	cp ../parameters/boundary.dat ./
	cp ../parameters/prior_types.dat prior_types.dat
	cp ../parameters/prior_params.dat prior_params.dat
	cp ../parameters/algo_config.dat algo_config.dat
	cp ../parameters/proposal_params.dat ./proposal_params.dat
	cp ../parameters/proposal_types.dat ./proposal_types.dat
	mpirun -np 240 ../../lwd_obj/gtmcmc_lwd
	
	
	mkdir samples
	mkdir logproposal
	mkdir loglike
	mkdir logprior
	mv samples.dat* ./samples/
	cd samples
	ls -ltr ./ | awk 'END{print $NF}'|xargs -i -t cp {} ../
	cd ../
	filename=$(find samples.dat.*)
    var=$(echo $filename | tr -cd "[0-9]")  
    outName="iter_"$i".txt"      
    echo $var >> ../iterations.txt
	mv samples.dat.* samples.dat
	python3 ../get_matrices.py
	mv samples.dat point$i.dat
	mv logproposal.dat* ./logproposal/
	mv loglik.dat* ./loglike/
	mv logprior.dat* ./logprior
	cp point$i.dat ../$outDir/samples
	cd ../
done
mv point* $outDir


