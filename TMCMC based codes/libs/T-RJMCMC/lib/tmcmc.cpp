#include "tmcmc.h"

/*
	Currently the best version with using covariance matrix.
*/

double tmcmc(int nsamps, int fixInit, double lambda,
          int iseed, int ndim, double cv,
          int write_flag, int pid, int np,
          double* generalPrior,
          double* Zrange, double jumpStd, double* BDprobs,
          LogDensity* logLikelihood) {
  /* reverse-jump TMCMC Algorithm
      Input: nsamps - number of Samples
             iseed - random seed
             fixInit - fixed initial dimensionality = ndim if 1
             lambda - parameter of exponential distribution for dimensionality
             ndim - dimensionality
             cv - Coefficient of Variance threshold for adapting Beta
             write_flag - provide intermediate TMCMC artifacts
             pid - MPI process ID
             np - number of MPI processes
             generalPrior - range of noise standard deviation, logRh, Zbed1, and thickness (vector of length 8)
             rwSigma - initial random walk step size (vector of length 8)
             Zrange - max and min of Zbed generated
             BDprobs - probability of layer birth, death, of unchange
             logLikelihood - likelihood instance
      Output: evid - asymptotically unbiased model evidence estimator
  */

  int MFactor = 1;
  
  int CATSteps = 1;
  int nspl = nsamps*np;
  int nSteps;
  int nsplSt;
  double covTime = 0;
  RealMatrix spls;
  RealVector spls_flat;
  RealVector jumpP;
  RealVector lprior;
  RealVector llik;
  RealVector phi;
  llik.resize(nsamps);
  lprior.resize(nsamps);
  jumpP.resize(nsamps);
  spls_flat.resize(nsamps*ndim);
  phi.resize(ndim);

  RealVector all_spls;
  RealVector all_lprior;
  RealVector all_llik;
  RealVector all_jumpP(nspl, 0);
  // all_spls.resize(nsamps*np*ndim);
  all_lprior.resize(nspl);
  all_llik.resize(nspl);


  RealVector splSt, llikSt, lpriorSt;
  
  RealVector cvmat;
  RealMatrix CVMATS;
  IntVector UniqueDims;
  IntMatrix smpIDs;
  RealVector splsComp;
  RealVector splSave, llikSave, lpriorSave;
  IntVector splCard;

  // For vectors of varying lenght
  // Initially, all samples have the same dimensionality
  IntVector splDims(nspl, ndim);
  IntVector splDimsNew;
  IntVector splDimsLocal;

  RealMatrix splNew(nspl);
  RealMatrix splNext(nspl);


  /* initialize random number generator */
  std::default_random_engine generator(iseed);
  std::uniform_real_distribution<double> u_distribution(0.0,1.0);
  std::normal_distribution<double> n_distribution(0.0,1.0);

  RealMatrix AccRatio(maxL-minL+1);
  // Roberts and Rosenthal 2011 - Initial gamma
  RealVector gm(maxL-minL+1);
  RealVector gm2(maxL-minL+1);
  for(int i = 0; i < maxL-minL+1; i++){
    gm[i] = 2.38/sqrt((i + minL)*2+1);
    gm2[i] = gm[i]*gm[i];
  }
  // double gm = 2.38 / sqrt(ndim);
  // double gm2 = gm * gm;

  //===========================================================================
  // Read prior parameters
  double* alphas = new double[ndim];
  double* betas = new double[ndim];
  char* types = new char[ndim];

  //===========================================================================
  // Generate prior parameters, local to each process
  // Initial parameters can have different dimensionality
  if (fixInit == 0)
  	genPriorSampsVarL(generator, u_distribution, spls, splDimsLocal, nsamps, generalPrior);
  else
  	genPriorSampsFixL(generator, u_distribution, spls, splDimsLocal, ndim, nsamps, generalPrior);
  int initTotal = std::accumulate(splDimsLocal.begin(), splDimsLocal.end(), 0);
  
  RealVector PriorMin_;
  RealVector PriorMax_;
  std::vector<char> Ptypes_;
  int nLayer = 0;
  sort_pair(spls, splDimsLocal);
  spls_flat = flatten(spls);
  for (int i = 0; i < nsamps; i++) {
    nLayer = (spls[i].size()-1)/2;
    full_prior(nLayer, generalPrior, Ptypes_, PriorMin_, PriorMax_);

    // -------------------------------------------------
    lprior[i] = LogPrior::eval_dim(spls[i],&PriorMin_[0],&PriorMax_[0],&Ptypes_[0], lambda);
    lprior[i] += log(factorial(nLayer-1));
    llik[i] = logLikelihood->eval(spls[i]);
  }
  // std::cout<<"Initial Prior and Likelihood computation Done\n";
  IntVector recvCount;
  IntVector Displs;
  int totalCount;
  if(pid == 0)
    recvCount.resize(np);
  // Gather the data length of each process
  MPI_Gather(&initTotal, 1, MPI_INT, &recvCount[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(pid == 0){
    splDims.resize(nspl);
    totalCount = std::accumulate(recvCount.begin(), recvCount.end(),0);
    all_spls.resize(totalCount);
    Displs.resize(np);
    Displs[0] = 0;
    for(int i = 1; i < np; i++)
      Displs[i] = Displs[i-1] + recvCount[i-1];
  }
  // Gather spls from all processes
  MPI_Gatherv(&spls_flat[0], initTotal, MPI_DOUBLE,
    &all_spls[0], &recvCount[0], &Displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&splDimsLocal[0],nsamps, MPI_INT, 
    &splDims[0], nsamps, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&lprior[0], nsamps, MPI_DOUBLE,
    &all_lprior[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&llik[0], nsamps, MPI_DOUBLE,
    &all_llik[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Output first set of samples to file*/
  if (pid == 0){
    // outProcToFile(all_spls,ndim,nspl,std::string("samples.dat.0"));
    int offset = 0;
    splNew.resize(nspl);
    for(int i = 0; i < nspl; i++){
      splNew[i].clear();
      for(int j = 0; j < splDims[i]; j++)
        splNew[i].push_back(all_spls[offset+j]);
      offset += splDims[i];
    }
    outProcToFile(splNew, std::string("samples.dat.0"));
    outProcToFile(all_llik,1,nspl,std::string("loglik.dat.0") );
    outProcToFile(all_lprior,1,nspl,std::string("logprior.dat.0"));
    std::cout<<"Process 0 gathers initial information Done\n";
  }

  RealVector Sm;
  std::vector<double> accRatio(maxL-minL+1,1.0);
  //double accRatio = 1.0;
  int iter = 0;
  bool check_final_iter;
  double beta = 0.0, dBeta = 0.0, evid = 0.0, dBetaFinal;
  double max_llik;
  double dBeta_prev = dBeta;
  double dBeta_slope;
  do { // Start algorithm
    iter++;
    
    if (pid == 0){
      /* shuffle samples */
      //shuffle_spls(all_spls,all_llik,all_lprior,generator);

      /* compute weights */
      RealVector w(nspl,0.0);
      RealVector w_norm(nspl,0.0);
      double wsum, wmean, w2mean, wstd;

      dBeta = std::min(BETA_MAX,1.0-beta);
      // dBeta_slope = -0.9;

      /* used to normalize log_likelihood values for weight computations */
      max_llik = *std::max_element(all_llik.begin(), all_llik.end());

      /* Adapt delta beta as needed */
      do {
        dBetaFinal = dBeta;
        for (int j=0; j < nspl; j++)
          w[j] = exp(dBeta*(all_llik[j]-max_llik));

        wsum   = std::accumulate(w.begin(), w.end(), 0.0);
        wmean  = wsum / w.size();
        w2mean = std::inner_product(w.begin(), w.end(), w.begin(), 0.0)/ w.size();
        wstd   = sqrt(w2mean- pow(wmean, 2));

        if (wstd/wmean > (cv + 1.0) || wstd == 0) dBeta *= 0.9;
        else if (wstd/wmean > (cv + 0.5) || wstd == 0) dBeta *= 0.95;
        else if (wstd/wmean > (cv + 0.05) || wstd == 0) dBeta *= 0.99;
        else if (wstd/wmean > (cv + 0.005) || wstd == 0) dBeta *= 0.999;
        else if (wstd/wmean > (cv + 0.0005) || wstd == 0) dBeta *= 0.9999;
        else if (wstd/wmean > (cv + 0.00005) || wstd == 0) dBeta *= 0.99999;

        if (dBeta < 1.0e-10)
          break;

      } while (wstd/wmean > (cv + 0.00005) || wstd == 0);
      dBeta = dBetaFinal;
      if (write_flag == 1){
        std::cout<<"DBeta: " << dBeta<<" Wmean: "<<wmean;
        std::cout<<" Wstd: "<<wstd<<" Cv: "<<wstd/wmean<<std::endl;
      }
     
      beta += dBeta;
      evid += log(wsum) + dBeta*max_llik - log(w.size());
      if (write_flag == 1){
        std::cout<<"Iteration "<<iter<<" Beta= "<<beta;
        std::cout<<" wMean= "<< wmean << " Evid=   " << evid<<std::endl<<std::flush;
      }

      /* Save mean ll for Bayes factor */
      Sm.push_back(wmean);

      /* rescale w and do cumulative sum */
      RealVector wb(nspl);
      std::transform(w.begin(),w.end(),w.begin(),
        std::bind(std::multiplies<double>(), 1.0/wsum, std::placeholders::_1));


      std::partial_sum(&w[0],&w[0]+nspl,&wb[0]);
      wb.insert ( wb.begin() , 0.0 );
      double t1 = MPI_Wtime();
      int nTypes; 
      std::cout<< "Computing cholesky decomposition..." << std::endl;
      UniqueDims = splDims;         
      std::sort(UniqueDims.begin(), UniqueDims.end());
      auto last = std::unique(UniqueDims.begin(), UniqueDims.end());
      UniqueDims.erase(last, UniqueDims.end()); // UniqueDims now contains all unique dimensions in splDims.
      nTypes = UniqueDims.size();
      CVMATS.resize(nTypes); // Covariance matrix for each unique dimensionality
      smpIDs.resize(nTypes);
      for(int i = 0; i < nTypes;i++)
        CVMATS[i].resize(UniqueDims[i]*UniqueDims[i]);
      for(int i = 0; i < nspl; i++)
        for(int j = 0; j < nTypes; j++){
          if(splDims[i] == UniqueDims[j]){
            smpIDs[j].push_back(i);
            break;
          }
        }

      RealMatrix Theta0(nTypes);
      for(int k = 0; k < nTypes; k++){
        Theta0[k].resize(UniqueDims[k]);
        for (int i=0; i<UniqueDims[k]; i++ ) {
          Theta0[k][i] = 0.0;
          for (int j=0; j< smpIDs[k].size(); j++ ) {
            Theta0[k][i] += splNew[smpIDs[k][j]][i]*w[smpIDs[k][j]];
          }
        }

        std::fill(CVMATS[k].begin(), CVMATS[k].end(), 0.0);
        for (int j=0; j < smpIDs[k].size(); j++) {
          for (int i1=0; i1<UniqueDims[k]; i1++) {
            for (int i2=0; i2<UniqueDims[k]; i2++) {
              double outp=w[smpIDs[k][j]] * (splNew[smpIDs[k][j]][i1]-Theta0[k][i1]) *
                                (splNew[smpIDs[k][j]][i2]-Theta0[k][i2]);
              CVMATS[k][i1*UniqueDims[k]+i2] += outp;
              CVMATS[k][i2*UniqueDims[k]+i1] += outp;
              }
          }
        }
        /* Control Parameter, Covariance rescaling */
        double gm2_c;
        for(int i = 0; i < maxL-minL+1; i++){
          if(UniqueDims[k] == (i+minL)*2+1){
            gm2_c = gm2[i];
            break;
          }
        }
        std::transform(CVMATS[k].begin(), CVMATS[k].end(), CVMATS[k].begin(),
        std::bind(std::multiplies<double>(), gm2_c, std::placeholders::_1));
        /* Cholesky factorization of the proposal covariance, in-place */
        int chol_info=0;
        char lu='L';
        cholesky(CVMATS[k], UniqueDims[k]);
      }
      double t2 = MPI_Wtime();
      covTime += t2-t1;

      /* generate random samples into [0,1] */
      RealVector spl01(nspl);
      for (int j=0; j<nspl; j++)
        spl01[j] = u_distribution(generator);

      /* get bin IDs and count */
      IntVector pos(nspl,0);
      for (int j=0; j < nspl; j++) {
        pos[j] = std::lower_bound(wb.begin(),wb.end(),spl01[j])-wb.begin()-1;
      }

      /* Count number of times a sample is "picked" by PRNG */
      IntVector splPos, splCount;
      for (int j=0; j < nspl; j++) {
        int icount=0;
        for (int ispl=0; ispl<nspl; ispl++) {
          if (pos[ispl]==j) icount += MFactor;
        }
        if (icount>0) {
          splPos.push_back(j);
          splCount.push_back(icount);
        }
      }

      /* Initialize samples that were retained, cardinality, and
      likelihood values */
      splCard.clear();
      nsplSt = splPos.size();
      splDimsNew.clear();

      splSt.clear();
      llikSt.clear();
      lpriorSt.clear();
      /* Resampling Step */
      int total_count = 0;

      for (int ispl=0; ispl<nsplSt; ispl++) {

        int isplCount = splCount[ispl];
        for (size_t i = 0; i < isplCount; ++i) {
        	splNext[total_count] = splNew[splPos[ispl]];
        	splDimsNew.push_back(splNext[total_count].size());         	          
          splCard.push_back(CATSteps);
          llikSt.push_back(all_llik[splPos[ispl]]);
          lpriorSt.push_back(all_lprior[splPos[ispl]]);
          total_count += 1;
        }
      }
      splDims = splDimsNew;
      splNew = splNext;
      if(iter>1)
      	splSt = flatten(splNext);
      nSteps = *std::max_element(splCard.begin(), splCard.end());
    }

    MPI_Bcast(&nSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&evid, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Run single steps of the Markov chains at a time, for the chains
        that need several jumps */
    IntVector candDim(nspl,0);
    IntVector candShift(nspl, 0);
    for (int isbSteps=0; isbSteps < nSteps; isbSteps++ ){
      if (pid == 0){
        nsplSt = nspl; // Post resampling, nspl # of chains
        // Sample from U(0,1), decide whether to change the dimensionality or not
        double jumpFactor;
        // birth or death or unchange
        char* BorD = new char[nsplSt];
        for(int i = 0; i < nsplSt; i++){
        	jumpFactor = u_distribution(generator);
        	BorD[i] = 'u';        	
			    if(BDprobs[0] < jumpFactor && jumpFactor < BDprobs[0]+BDprobs[1] && splDims[i] != maxL*2+1 ){ 	
    			  BorD[i] = 'b';	
        	}
        	else if(jumpFactor > BDprobs[0]+BDprobs[1] && splDims[i] != minL*2+1 ){
    			  BorD[i] = 'd';
        	}
        }
        
        RealMatrix splCand(nsplSt);

        int dataShift = 0;
        std::cout<<"Generating candidate samples for iteration "<< iter << std::endl;
        for (int ispl = 0; ispl < nsplSt; ispl++) {
          splCand[ispl].resize(splDims[ispl]);
          /* generate candidate */
          if(BorD[ispl] != 'u'){
          	all_jumpP[ispl] = birthDeath(splCand[ispl], &splNew[ispl][0], BorD[ispl], Zrange[0], Zrange[1], 
          						jumpStd, generator, lambda, n_distribution);
          	candDim[ispl] = splCand[ispl].size();
          }
          else{
            for(int i = 0; i < UniqueDims.size(); i++){
                if(splDims[ispl] == UniqueDims[i]){
                  cvmat = CVMATS[i];
                  break;
                }
            }
            all_jumpP[ispl] = 0;
          	candDim[ispl] = splDims[ispl];
          	RealVector xi(candDim[ispl]);
          	double Lnrv = 0.0;
          	int nLayersCand = 0;
          	for (size_t i=0; i < candDim[ispl]; i++)  
            {
              xi[i] = n_distribution(generator);
              splCand[ispl][i] = splNew[ispl][i];
              Lnrv=0.0;
              for (size_t j=0; j < (i+1); ++j) {
                    Lnrv += cvmat[j*candDim[ispl]+i] * xi[j];
                }
                splCand[ispl][i] += Lnrv;
            }
          }
          dataShift += splDims[ispl];
          if(ispl > 0)
          	candShift[ispl] = candShift[ispl-1] + candDim[ispl-1];
		    /* done generating candidate */
        }
        int candLenFlat = std::accumulate(candDim.begin(), candDim.end(), 0);
        // RealVector splCandFlat;
        // splCandFlat = flatten(splCand);
        // outProcToFile(splCand, std::string("candidate.dat.")+std::to_string(iter));
        /* Compute new likelihoods */
        splsComp.clear();
        // int compCount=0;
        for (int ispl=0; ispl<nsplSt; ispl++) {
          for (int i=0; i<splCand[ispl].size(); i++ ) {
            splsComp.push_back(splCand[ispl][i]);
          }
          // compCount++;
        }
      }
      // MPI_Bcast(&candShift[0], nspl, MPI_INT, 0, MPI_COMM_WORLD);
      // MPI_Scatter(&all_jumpP[0], nsamps, MPI_DOUBLE, &jumpP[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      //MPI_Scatter(&candDim[0], nsamps, MPI_INT, &candDimLocal[0], nsamps, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&candDim[0], nspl, MPI_INT, 0, MPI_COMM_WORLD);
      IntVector candDimLocal(&candDim[pid*nsamps], &candDim[(pid+1)*nsamps]);

  	  int* sendCounts = new int[np];
  	  int* displs = new int[np];
  	  int offset1 = 0;
  	  for(int i = 0; i < np; i++){
  	  	displs[i] = offset1;
  	  	sendCounts[i] = std::accumulate(&candDim[i*nsamps], &candDim[(i+1)*nsamps], 0);
  	  	offset1+=sendCounts[i];
  	  		
  	  }
  	  

  	  spls_flat.resize(sendCounts[pid]);

  	  MPI_Barrier(MPI_COMM_WORLD);
      MPI_Scatterv(&splsComp[0], sendCounts, displs, MPI_DOUBLE,
          &spls_flat[0], sendCounts[pid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      

      int offset = 0;
      RealVector PriorMin;
      RealVector PriorMax;
      int nLayer;
      std::vector<char> Ptypes;

      for (int i = 0; i < nsamps; i++) {
      	phi.resize(candDimLocal[i]);
        for (int j = 0; j < candDimLocal[i]; j++ ){
          phi[j] = spls_flat[offset+j];
        }
        nLayer = (phi.size()-1)/2;
        offset += candDimLocal[i];
        llik[i] = logLikelihood->eval(phi);
        full_prior(nLayer, generalPrior, Ptypes, PriorMin, PriorMax);
        // -------------------------------------------------
        lprior[i] = LogPrior::eval_dim(phi,&PriorMin[0],&PriorMax[0],&Ptypes[0],lambda);
        lprior[i] += log(factorial(nLayer-1));
        // -------------------------------------------------
        
      }

      MPI_Gather(&lprior[0], nsamps, MPI_DOUBLE,
       &all_lprior[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&llik[0], nsamps, MPI_DOUBLE,
        &all_llik[0], nsamps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (pid == 0){
      	
        /* decide who jumps */
        int icomp=0;
        
        int nTypes; 
        nTypes = maxL - minL + 1; // all possible number of layers
        IntVector acceptCount(nTypes, 0);
        IntVector nSplsEachDim(nTypes, 0);
        RealVector llikNew(nsplSt), lpriorNew(nsplSt); // number of samples for each category
        for (int ispl=0; ispl<nsplSt; ispl++) {
          // Count number of samples in each category
          for(int i = minL; i < maxL+1; i++){
            if(candDim[ispl] == i*2+1)
              nSplsEachDim[i-minL]++;
          }
          double alpha = u_distribution(generator);
          double AcceptRatio = -1.0;
          AcceptRatio = beta * (all_llik[icomp] - llikSt[ispl] + all_lprior[icomp] - lpriorSt[ispl] + all_jumpP[ispl]);                       

          if (log(alpha) < AcceptRatio) { // Accept proposal
          	splNew[ispl].resize(candDim[ispl]);
          	// std::cout<<"sample dim: "<< splNew[ispl].size() << std::endl;
          	splDims[ispl] = candDim[ispl];
            for (int i=0; i < candDim[ispl]; i++) {
              splNew[ispl][i] = splsComp[candShift[icomp]+i];
            }
            lpriorNew[ispl] = all_lprior[icomp];
            llikNew[ispl] = all_llik[icomp];
            // Number of samples accepted in each category
            for(int i = minL; i < maxL+1; i++){
              if(splDims[ispl] == i*2+1)
                acceptCount[i - minL]++;
            }
            

          } else { // Reject Proposal
            lpriorNew[ispl] = lpriorSt[ispl];
            llikNew[ispl] = llikSt[ispl];
          }
          icomp++;
          if(splNew[ispl].size()!= splDims[ispl])
          	std::cout<<"Error: dimensionality doesn't match" << std::endl;
        }

        outProcToFile(splNew,
                      std::string("samples.dat.")+std::to_string(iter));
        // /* Clear Proposals for next iteration */
        splSt = flatten(splNew);
        llikSt = llikNew;
        lpriorSt = lpriorNew;
        

        /* Reduce length of chains remaining in samples */
        for (int ispl=0; ispl<splCard.size();) {
          if (splCard[ispl]==0)
            splCard.erase(splCard.begin()+ispl) ;
          else
            ispl++;
        }

        double G = 2.1;
        for(int i = 0; i < nTypes; i++){
          if(nSplsEachDim[i] != 0 ){
            accRatio[i] = (double) acceptCount[i] / (double) nSplsEachDim[i];
            gm[i] = gm[i] * exp(G * (accRatio[i] - 0.234));
            gm2[i] = pow(gm[i], 2.0);
          }
          else
            accRatio[i] = 2;
        }

        nsplSt = llikSt.size();
        assert(splCard.size()==llikSt.size());
        //assert(splSt.size()==llikSt.size()*ndim);
      }
    }

    if (pid == 0){
      // Rescaling based on Catanach (Thesis 2017)
      // Assumptions made on optimal acceptance rate
      //std::cout<<" Accept ratio: "<< accRatio << std::endl;
      // double G = 2.1;
      // if(accRatio < accThreshold){
      // 	for(int i = 0; i < 3; i++)
    		// rwSigma[i] = rwSigma[i]*exp(G * (accRatio - 0.234));
      // }

      
     //  for(int i = 0; i < 3; i++)
    	// rwSigma[i] = rwSigma[i]*exp(G * (accRatio - 0.234));

      /* Set samples for next temperature iteration */
      for(int i = 0; i < maxL-minL+1; i++)
        AccRatio[i].push_back(accRatio[i]);
      all_spls = splSt;
      all_llik = llikSt;
      all_lprior = lpriorSt;
      assert(all_llik.size()==nspl);
      assert(all_lprior.size()==nspl);
      //assert(all_spls.size()==nspl*ndim);    
      outProcToFile(llikSt,1,   nspl,
                      std::string("loglik.dat." )+std::to_string(iter));
      std::string dimNanme = std::string("dimensions.dat.")+std::to_string(iter);
      FILE* fp=fopen(const_cast<char*>(dimNanme.c_str()),"w");
      for(int i = 0; i < nspl; i++)
        fprintf(fp,"%d\n", splDims[i]);
      fclose(fp);
    }

    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  } while ( ( iter < MAX_TMCMC_ITER ) && (beta<1-1.e-10) );

  if (pid == 0){
    if (write_flag == 1){
      std::cout << "TMCMC Algorithm Done" << std::endl;
    }
    std::cout<<"Time for covariance matrix computing"<< covTime<<std::endl;
  for(int i = 0; i < maxL-minL+1; i++)
    outProcToFile(AccRatio[i],1, iter,std::string("AccRatio.dat" ) + std::to_string(i+minL));
  }
  
  return (evid);
} 
/* done tmcmc */

void cholesky(RealVector &A, int n) {
  RealVector L;
  L.resize(n*n);

  for (int i = 0; i < n; i++)
      for (int j = 0; j < (i+1); j++) {
          double s = 0;
          for (int k = 0; k < j; k++)
              s += L[i * n + k] * L[j * n + k];
          L[i * n + j] = (i == j) ?
                         sqrt(A[i * n + i] - s) :
                         (1.0 / L[j * n + j] * (A[i * n + j] - s));
      }

  for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
          A[i * n + j] = L[i * n + j];
      }
}

RealVector flatten(const RealMatrix& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size();
    RealVector result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}




void readPriorParams(double *alphas, double *betas, char *types, int ndim, std::string priorParamsFile, std::string priorTypesFile) {
  std::string line;
  std::ifstream DAT;
  double temp;
  std::string::size_type sz;
  int i = 0;

  // Prior types
  DAT.open(priorTypesFile);
  while (std::getline(DAT, line)) {
      types[i] = line.at(0);
      i = i+1;
  }
  DAT.close();

  // Prior parameters
  i = 0;
  DAT.open(priorParamsFile);
  while (std::getline(DAT, line)) {
      try {
        temp = std::stod (line,&sz);
      } catch (const std::exception& ex) {
        break;
      }
      line = line.substr(sz);
      alphas[i] = temp;
      try {
        temp = std::stod (line,&sz);
      } catch (const std::exception& ex) {
        break;
      }
      line = line.substr(sz);
      betas[i] = temp;
      i = i+1;
  }
  DAT.close();
}


void genPriorSampsFixL(std::default_random_engine &generator,
  std::uniform_real_distribution<double> &u_distribution,
  RealMatrix &spls, IntVector &splDim, int dim,
  int nsamps, double *generalPrior) {

  spls.resize(nsamps);
  splDim.resize(nsamps);
  int nL;
  for (int i = 0; i < nsamps; i++){
    splDim[i] = dim;
    nL = (dim-1)/2;
    spls[i].resize(splDim[i]);
    for (int j = 0; j < 2; j++) {
      spls[i][j] = generalPrior[0] + (generalPrior[1]-generalPrior[0])*u_distribution(generator);
    }
    for (int j = 2; j < 2+nL; j++) {
      spls[i][j] = generalPrior[2] + (generalPrior[3]-generalPrior[2])*u_distribution(generator);
    }
    spls[i][2+nL] = generalPrior[4] + (generalPrior[5]-generalPrior[4])*u_distribution(generator);
    for (int j = 3+nL; j < splDim[i]; j++) {
      spls[i][j] = generalPrior[6] + (generalPrior[7]-generalPrior[6])*u_distribution(generator);
    }
  }
}

void genPriorSampsVarL(std::default_random_engine &generator,
  std::uniform_real_distribution<double> &u_distribution,
  RealMatrix &spls, IntVector &splDim, 
  int nsamps, double *generalPrior) {

  spls.resize(nsamps);
  splDim.resize(nsamps);
  int nL;
  for (int i = 0; i < nsamps; i++){
    if(maxL == minL)
      nL = minL;
    else
    {
      nL = rand()%(maxL - minL);
      nL += minL;
    }
    splDim[i] = nL*2 + 1;
    spls[i].resize(splDim[i]);
    for (int j = 0; j < 2; j++) {
      spls[i][j] = generalPrior[0] + (generalPrior[1]-generalPrior[0])*u_distribution(generator);
    }
    for (int j = 2; j < 2+nL; j++) {
      spls[i][j] = generalPrior[2] + (generalPrior[3]-generalPrior[2])*u_distribution(generator);
    }
    spls[i][2+nL] = generalPrior[4] + (generalPrior[5]-generalPrior[4])*u_distribution(generator);
    for (int j = 3+nL; j < splDim[i]; j++) {
      spls[i][j] = generalPrior[6] + (generalPrior[7]-generalPrior[6])*u_distribution(generator);
    }
  }
}

void outProcToFile(const RealVector spls, const int ndim, const int
nspl, std::string fname) {

  FILE *myfile  = fopen(fname.c_str(),"w") ;
  for (int j = 0; j < nspl; j++) {
    for (int i = 0; i < ndim; i++)
      fprintf(myfile,"%24.18e ",spls[j*ndim+i]);
    fprintf(myfile,"\n");
  }
  fclose(myfile);

  return ;

}

void outProcToFile(const RealMatrix& spls, std::string fname) {

  FILE *myfile  = fopen(fname.c_str(),"w") ;
  // int offset = 0;
  for (int j = 0; j < spls.size(); j++) {
    for (int i = 0; i < spls[j].size(); i++)
      fprintf(myfile,"%24.18e ",spls[j][i]);
    fprintf(myfile,"\n");
    // offset += ndim[nspl];
  }
  fclose(myfile);

  return ;

}

void shuffle_spls(RealVector &spls, RealVector &llik, RealVector &lprior, std::default_random_engine &generator) {
  // Shuffle the samples randomly
  int nspl = llik.size();
  int ndim = spls.size()/nspl;

  IntVector idx(nspl);
  for (int j=0; j<nspl; j++) idx[j]=j;

  RealVector splsTmp(nspl*ndim), llikTmp(nspl), lpriorTmp(nspl);
  shuffle (idx.begin(), idx.end(), generator);
  for (int j = 0; j < nspl; j++) {
    llikTmp[j] = llik[idx[j]];
    lpriorTmp[j] = lprior[idx[j]];
    for (int i = 0; i < ndim; i++) splsTmp[j*ndim+i] = spls[idx[j]*ndim+i];
  }

  for (int j = 0; j < nspl; j++) {
    llik[j] = llikTmp[j];
    lprior[j] = lpriorTmp[j];
    for (int i = 0; i < ndim; i++) spls[j*ndim+i] = splsTmp[j*ndim+i];
  }

  return ;

}


double birthDeath(RealVector& splCand, double* spls, const char BorD, const int Zmin, const int Zmax, 
					double jumpStd, std::default_random_engine &generator, double lambda,
					std::normal_distribution<double> &n_distribution){
	// Dimensionality of previous sample
	int nDim = splCand.size();
	int nLayers = (nDim-1)/2;
	int nLayersNew;
	double* Zbeds = new double[nLayers-1];
    double Pjump = 0; // q(m|m*)/q(m*|m)
	Zbeds[0] = spls[nLayers+2];
	for(int i = 1; i < nLayers-1; i++){
		Zbeds[i] = Zbeds[i-1] + spls[nLayers+2 + i];
	}
	if(BorD == 'b'){
		nLayersNew = nLayers+1;
		int ZnewOrder = 0;
		splCand.resize(nDim+2);
		int Znew = rand()%(Zmax - Zmin);
		Znew += Zmin; 
		for(int i = 0; i < nLayers-1; i++){
			if(Znew > Zbeds[i])
				ZnewOrder += 1;
			else
				break;
		}
		int j = 0;
		splCand[0] = spls[0];
		splCand[1] = spls[1];
		for(int i = 0; i < nLayersNew; i++){
			if(i == ZnewOrder){
        double jumpMean = splCand[2+i];
        double tempP = n_distribution(generator);
				splCand[2+i] = spls[2+i]+jumpStd*tempP;
				j = -1;
        // Pjump = (7-double(nLayers))/(double(nLayers)+1);   
        char jumpType = 'G';
        RealVector jumpV(1);
        jumpV[0] = splCand[2+i];
        tempP = LogPrior::eval_dim(jumpV,&jumpMean,&jumpStd,&jumpType,lambda);
        Pjump = double(Zmax - Zmin)/(double(nLayers+1)*exp(tempP)); 
        // Pjump = Pjump/exp(tempP);
        // Pjump = 1./(6*exp(tempP)); 

        // if(Pjump < 0)
        //   std::cout << "Error: negative probability. \n";
			}
			else{
				splCand[2+i] = spls[2+i+j];
			}
		}
		j = 0;
		for(int i = 0; i < nLayersNew-1; i++){
			if(i == ZnewOrder){
				splCand[2+nLayersNew+i] = Znew;
				j = -1;
			}
			else
				splCand[2+nLayersNew+i] = Zbeds[i+j];
		}


	}
	else{
		nLayersNew = nLayers-1;
		splCand.resize(nDim-2);
		int Zremove = rand()%(nLayers-1);
		splCand[0] = spls[0];
		splCand[1] = spls[1];
		int j = 0;

		for(int i = 0; i < nLayersNew; i++){
			if(i == Zremove){
        double jumpMean = splCand[2+i];
        double tempP;
        char jumpType = 'G';
        RealVector jumpV(1);
        
				splCand[2+i] = (spls[2+i]+spls[2+i+1])/2;

				j = 1;
        jumpV[0] = splCand[2+i];
        tempP = LogPrior::eval_dim(jumpV,&jumpMean,&jumpStd,&jumpType,lambda);
        // Pjump = double(nLayers)*exp(tempP)/double(Zmax - Zmin);
        // Pjump = 6*exp(tempP);
        Pjump = double(nLayers)*exp(tempP)/double(Zmax - Zmin);
        // if(Pjump < 0)
        //   std::cout << "Error: negative probability. \n";
			}
			else
				splCand[2+i] = spls[2+i+j];
		}
		j = 0;
		for(int i = 0; i < nLayersNew-1; i++){
			if(i == Zremove){
				splCand[2+nLayersNew+i] = Zbeds[i+1];
				j = 1;
			}
			else{
				splCand[2+nLayersNew+i] = Zbeds[i+j];
			}
		}
	}
	if(nLayersNew > 2)
		for(int i = nLayersNew*2; i > nLayersNew+2; i--){
			splCand[i] = splCand[i] - splCand[i-1];
			if(splCand[i] < 0)
				std::cout <<"Error: negative thickness generated." <<std::endl;
		}

  return log(Pjump);
}

void sort_pair(RealMatrix& spls, IntVector& dims){
	IntVector V(dims.size());
	std::iota(V.begin(),V.end(),0);
	std::sort(V.begin(),V.end(), [&](int i,int j){return dims[i]<dims[j];} );
	RealMatrix tempSpls(spls.size());
	IntVector tempDims(spls.size());
	for(int i = 0; i < dims.size(); i++){
		tempSpls[i] = spls[V[i]];
		tempDims[i] = dims[V[i]];
	}
	dims = tempDims;
	spls = tempSpls;
}

void full_prior(int nLayers, double* generalPrior, std::vector<char> &Ptypes, 
	RealVector &PriorMin, RealVector &PriorMax){
	PriorMin.clear();
	PriorMax.clear();
	Ptypes.clear();
	int nP = nLayers*2 + 1;
	for(int i = 0 ; i < nP; i++)
        Ptypes.push_back('U');
    PriorMin.push_back(generalPrior[0]);
    PriorMin.push_back(generalPrior[0]);
    PriorMax.push_back(generalPrior[1]);
    PriorMax.push_back(generalPrior[1]);
    for(int i = 2; i < nLayers+2; i++){
      PriorMin.push_back(generalPrior[2]);
      PriorMax.push_back(generalPrior[3]);
    }
    PriorMin.push_back(generalPrior[4]);
    PriorMax.push_back(generalPrior[5]);
    for(int i = nLayers+3; i < nP; i++){
      PriorMin.push_back(generalPrior[6]);
      PriorMax.push_back(generalPrior[7]);
    }
}

int factorial(int n){
  if(n == 0)
    return 0;
  int result = 1;
  for(int i = 1; i <= n; i++)
    result *= i;
  return result;
}