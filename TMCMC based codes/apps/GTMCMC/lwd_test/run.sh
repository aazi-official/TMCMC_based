
module load openmpi/1.10.7
module load python/3.7

outDir=outDir
mkdir $outDir
mkdir $outDir/samples
mkdir $outDir/loglike
mkdir $outDir/logproposal
for i in {1..2}
do
	mkdir point$i
	cd point$i
	sed $i,$i'!d' ../measurements/layer_n_dip.dat >> measurements.dat
	sed $i,$i'!d' ../measurements/measurements.dat >> measurements.dat

	cp ../parameters/boundary.dat ./
	cp ../parameters/prior_types.dat prior_types.dat
	cp ../parameters/prior_params.dat prior_params.dat
	cp ../parameters/algo_config.dat ./
	if [ $i != 1 ]
	then
		prevID=`expr $i - 1`
		cp ../point$prevID/spls_cov_sqrt.dat ./proposal_cov_sqrt.dat
		cp ../point$prevID/spls_cov_inv.dat ./proposal_cov_inv.dat
		cp ../point$prevID/spls_mean.dat ./proposal_mean.dat
		# Need additional importance distribution parametesr
		mpirun -np 240 ../../lwd_obj/gtmcmc_cov_lwd

	else
		cp ../parameters/prior_params.dat ./proposal_params.dat
		cp ../parameters/prior_types.dat ./proposal_types.dat
		mpirun -np 240 ../../lwd_obj/gtmcmc_lwd
	fi
	
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


