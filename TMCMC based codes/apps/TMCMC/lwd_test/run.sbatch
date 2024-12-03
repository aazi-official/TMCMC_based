#!/bin/bash

#SBATCH -J tmcmc     # Job name
#SBATCH -o tmcmc%j.log # Slurm out file name

#SBATCH -N 15                  # number of nodes
#SBATCH -p xfu
#SBATCH --ntasks-per-node=16    # processors per node
#SBATCH -t 168:00:00              # time in hours
#SBATCH --mem-per-cpu=2048      # number of MB per processor
#SBATCH --mail-user=your_email_address
#SBATCH --mail-type=begin        # email me when the job starts
#SBATCH --mail-type=end          # email me when the job finishes

# How to conver DOS line breakers to Linux line breakers:
# vi run_tmcmc.sbatch -c "set ff=unix" -c ":wq"
module load openmpi/1.10.7

outDir=outDir
mkdir $outDir
mkdir $outDir/samples
mkdir $outDir/loglike

for i in {1..1}
do
	mkdir point_$i
	cd point_$i
    echo "0.2" > cv.dat
    cp ../parameters/prior_params.dat ./prior_params.dat
    cp ../parameters/prior_types.dat ./prior_types.dat
	sed $i,$i'!d' ../measurements/layer_n_dip_5L_1.dat >> measurements.dat
	sed $i,$i'!d' ../measurements/measurements_5L_1.dat >> measurements.dat
	cp ../parameters/boundary.dat ./
	cp ../parameters/algo_config.dat ./
	mpirun -np 240 ../../lwd_obj/tmcmc_lwd

	mkdir samples
	mkdir loglike
	mkdir logprior
	mv samples.dat* ./samples/
	cd samples
	ls -ltr ./ | awk 'END{print $NF}'|xargs -i -t cp {} ../
	cd ../
	mv loglik.dat* ./loglike/
	cd loglike
	ls -ltr ./ | awk 'END{print $NF}'|xargs -i -t cp {} ../
	cd ../
	filename=$(find samples.dat.*)
    var=$(echo $filename | tr -cd "[0-9]")       
    echo $var >> ../iterations.txt
	mv samples.dat.* samples.dat
	mv samples.dat point$i.dat
	mv loglik.dat* ../$outDir/loglike/point$i.dat
	cp point$i.dat ../$outDir/samples
	cd ../
done
mv point_* $outDir


