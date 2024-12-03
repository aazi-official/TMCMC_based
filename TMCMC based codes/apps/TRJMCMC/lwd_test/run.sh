module load openmpi/1.10.7
module load python/3.7

outDir=outputDir
mkdir $outDir
mkdir $outDir/samples
mkdir $outDir/loglike
mkdir $outDir/dimensions

for i in {1..1}
do
	mkdir temppoint_$i
	cd temppoint_$i
	# number of layers and dip angle
	echo "5 85" >> measurements.dat
	# the measurements_file.dat contains measurements
	# each row for one problem
	# fetch the i-th row
	sed $i,$i'!d' ../measurements/measurements_5L_10m.dat >> measurements.dat
	# hyper-parameter
	cp ../parameters/rjSettings_11.dat ./rjSettings.dat
	# Zrange used to generate new boundaries
	cp ../parameters/Zrange.dat ./Zrange.dat
	# Prior for each parameter
	cp ../parameters/generalPrior.dat ./generalPrior.dat
	cp ../parameters/algo_config.dat ./
	# parameters: dim cv nsample accThreshold fixInitFlag lambda             
	#echo "11 0.2 50 0.1 0 0 1" > algo_config.dat
	
	mpirun -np 240 ../../lwd_obj/rjTMCMC
	
	# post-processing
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
	mv dimensions.dat ../$outDir/dimensions/dimensions$i.dat
	mv loglik.dat* ../$outDir/loglike/point$i.dat
	cp point$i.dat ../$outDir/samples
	# cp dimensions$i.dat ../$outDir/dimensions
	cd ../
done
mv temppoint_* $outDir
