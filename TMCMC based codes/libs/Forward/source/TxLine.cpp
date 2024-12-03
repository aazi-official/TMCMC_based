#include "TxLine.h"
#include <iostream>
#include <algorithm>

CTxLine::CTxLine()
{
}

CTxLine::~CTxLine()
{
}

bool CTxLine::calcTxLineParas(	in_ const SFormationConfig & ForConfig, in_ const DataType & Omega,	in_ const DataType & krho, out_ STxLineDataType & TxLinePara )
{
	// calculate chracteristic impeadance for each layer
	for ( unsigned int i = 0; i <= ForConfig.s_iNum_Layers+1; i++ )
	{
		Complex kz_tmp = std::sqrt( ForConfig.s_vHWaveNum[i] * ForConfig.s_vHWaveNum[i] - krho * krho * ForConfig.s_vTMAnsiRatio[i] );

		// check the logic operator to be or / and
		if ( kz_tmp.real() < 0.0 || kz_tmp.imag() > 0.0 )
		{
			kz_tmp = - kz_tmp;
		}
		TxLinePara.s_vWaveNumber_z[0][i] = kz_tmp;
		TxLinePara.s_vCharImpdZ[0][i] = kz_tmp/Omega/ForConfig.s_vHCompEpsilon[i];

		kz_tmp = std::sqrt( ForConfig.s_vHWaveNum[i] * ForConfig.s_vHWaveNum[i] - krho * krho * ForConfig.s_vTEAnsiRatio[i] );

		if ( kz_tmp.real() < 0.0 || kz_tmp.imag() > 0.0 )
		{
			kz_tmp = - kz_tmp;
		}
		TxLinePara.s_vWaveNumber_z[1][i] = kz_tmp;
		TxLinePara.s_vCharImpdZ[1][i] = Omega * ForConfig.s_vHCompMu[i]/kz_tmp;
	}
	//std::cout << TxLinePara.s_vCharImpdZ[1][28] << std::endl;

	for ( unsigned int iMode = 0; iMode < 2; iMode++ )
	{
		TxLinePara.s_vPhase[iMode][0] = CPLX_0;
		TxLinePara.s_vPhase[iMode][0] = CPLX_0;
		TxLinePara.s_vPhase2[iMode][1] = CPLX_0;
		TxLinePara.s_vPhase2[iMode][1] = CPLX_0;

		for ( unsigned int i = 2; i <= ForConfig.s_iNum_Layers; i++ )
		{
			TxLinePara.s_vPhase[iMode][i] = exp( - CPLX_J * TxLinePara.s_vWaveNumber_z[iMode][i]* ForConfig.s_vLayerThick[i] );
			TxLinePara.s_vPhase2[iMode][i] = TxLinePara.s_vPhase[iMode][i] * TxLinePara.s_vPhase[iMode][i];
		}
		
		TxLinePara.s_vPhase[iMode][ForConfig.s_iNum_Layers] = CPLX_0;
		TxLinePara.s_vPhase[iMode][ForConfig.s_iNum_Layers] = CPLX_0;
		TxLinePara.s_vPhase2[iMode][ForConfig.s_iNum_Layers+1] = CPLX_0;
		TxLinePara.s_vPhase2[iMode][ForConfig.s_iNum_Layers+1] = CPLX_0;

		// calculate the Fresnel Reflection Coefficient
		for ( unsigned int i = 0; i <= ForConfig.s_iNum_Layers; i ++ )
		{
			TxLinePara.s_vFrensnelRefCoef[iMode][i]  = ( TxLinePara.s_vCharImpdZ[iMode][i+1] - TxLinePara.s_vCharImpdZ[iMode][i] )
													 / ( TxLinePara.s_vCharImpdZ[iMode][i+1] + TxLinePara.s_vCharImpdZ[iMode][i] ) ;
		}
		//std::cout << TxLinePara.s_vCharImpdZ[iMode][28] << std::endl;
		//std::cout << TxLinePara.s_vCharImpdZ[iMode][27] << std::endl;
		//std::cout << TxLinePara.s_vFrensnelRefCoef[iMode][27] << std::endl;
		// calculate the Generalized Reflection Coefficient
		//TxLinePara.s_vRefCoefUp[iMode][i] = TLParameters.s_FrensnelRefCoef[iMode][iNumofLayers];
		for ( int i = 2; i <= ForConfig.s_iNum_Layers; i++  )
		{
			TxLinePara.s_vRefCoefUp[iMode][i] = ( - TxLinePara.s_vFrensnelRefCoef[iMode][i-1] + TxLinePara.s_vRefCoefUp[iMode][i-1] * TxLinePara.s_vPhase2[iMode][i-1] )
											  / ( CPLX_ONE - TxLinePara.s_vRefCoefUp[iMode][i-1] * TxLinePara.s_vFrensnelRefCoef[iMode][i-1] * TxLinePara.s_vPhase2[iMode][i-1]);
		}

		//TLParameters.s_RefCoefDown[iMode][1] = - TLParameters.s_FrensnelRefCoef[iMode][0];
		for ( int i = ForConfig.s_iNum_Layers-1; i >= 1; i-- )
		{
			TxLinePara.s_vRefCoefDown[iMode][i] = ( TxLinePara.s_vFrensnelRefCoef[iMode][i] + TxLinePara.s_vRefCoefDown[iMode][i+1] * TxLinePara.s_vPhase2[iMode][i+1] ) 
												/ ( CPLX_ONE + TxLinePara.s_vRefCoefDown[iMode][i+1] * TxLinePara.s_vFrensnelRefCoef[iMode][i] * TxLinePara.s_vPhase2[iMode][i+1]);
			//std::cout << TxLinePara.s_vFrensnelRefCoef[iMode][i] << std::endl;
			//std::cout << TxLinePara.s_vRefCoefDown[iMode][i+1] << std::endl;
			//std::cout << TxLinePara.s_vPhase2[iMode][i+1] << std::endl;
			//std::cout << TxLinePara.s_vRefCoefDown[iMode][i+1] << std::endl;
			//std::cout << TxLinePara.s_vFrensnelRefCoef[iMode][i] << std::endl;
			//std::cout << TxLinePara.s_vPhase2[iMode][i+1] << std::endl;
		}
		for ( unsigned int i = 2; i <= ForConfig.s_iNum_Layers-1; i++ )
		{
			TxLinePara.s_vTransCoefDown[iMode][i] = ( CPLX_ONE + TxLinePara.s_vRefCoefDown[iMode][i] ) * TxLinePara.s_vPhase[iMode][i]
												  / ( CPLX_ONE + TxLinePara.s_vRefCoefDown[iMode][i] * TxLinePara.s_vPhase2[iMode][i] );
		}
	}
	//std::cout << TxLinePara.s_vCharImpdZ[1][10] << std::endl;
	//std::cout << TxLinePara.s_vPhase[1][10] << std::endl;
	//std::cout << TxLinePara.s_vRefCoefUp[1][10] << std::endl;
	//std::cout << TxLinePara.s_vRefCoefDown[1][10] << std::endl;
	//std::cout << TxLinePara.s_vTransCoefDown[1][10] << std::endl;
	return true;
}

bool CTxLine::calcTxLineVIs(	in_ const SFormationConfig & ForConfig,	in_ const STxLineDataType & TxLinePara,
								in_ const SObsSrcPair OSPair,			in_ const DataType & krho,
								out_ STxLineVI & VI )
{
	Complex phase1[2] = {0.0};
	Complex phase2Term[2] = {0.0};
	Complex TransDown[2] = {0.0};
	Complex VCoeff[2] = {0.0};
	Complex ICoeff[2] = {0.0};
	Complex Tmp;
	STxLineVI VISrc;

	if ( OSPair.s_iSrcLayer == OSPair.s_iObsLayer )
	{
		calTLVIatSrcLayer( ForConfig, TxLinePara, OSPair.s_iSrcLayer, OSPair.s_Obs_z, OSPair.s_Src_zPrim, VI );

		//std::cout << VI.s_Iv[0] << VI.s_Ii[1] << VI.s_Iv[1] << VI.s_Vi[1] << VI.s_Vv[1] << std::endl;
	}
	else
	{
		DataType maxZ = std::max( OSPair.s_Obs_z, OSPair.s_Src_zPrim );
		DataType minZ = std::min( OSPair.s_Obs_z, OSPair.s_Src_zPrim );

		unsigned int maxLayer = std::max( OSPair.s_iObsLayer, OSPair.s_iSrcLayer );
		unsigned int minLayer = std::min( OSPair.s_iObsLayer, OSPair.s_iSrcLayer );

		calTLVIatSrcLayer( ForConfig, TxLinePara, minLayer, ForConfig.s_vLayerBoundary[minLayer], maxZ,	VISrc );
		
		//std::cout << VI.s_Iv[0] << VI.s_Ii[1] << VI.s_Iv[1] << VI.s_Vi[1] << VI.s_Vv[1] << std::endl;
		
		for ( unsigned int iMode = 0; iMode < 2; iMode ++ )
		{
			TransDown[iMode] = CPLX_ONE;
			for ( unsigned int i = minLayer+1; i <= maxLayer-1; i++ )
			{
				TransDown[iMode] = TransDown[iMode] * TxLinePara.s_vTransCoefDown[iMode][i];
			}
			
			phase1[iMode] = exp( - CPLX_J * TxLinePara.s_vWaveNumber_z[iMode][maxLayer] * ( ForConfig.s_vLayerBoundary[maxLayer-1] - minZ ));
			phase2Term[iMode] = exp( - DataType(2.0) * CPLX_J * TxLinePara.s_vWaveNumber_z[iMode][maxLayer] * ( minZ - ForConfig.s_vLayerBoundary[maxLayer] ))
							  * TxLinePara.s_vRefCoefDown[iMode][maxLayer];

			TransDown[iMode] = TransDown[iMode]/(CPLX_ONE + TxLinePara.s_vRefCoefDown[iMode][maxLayer]*TxLinePara.s_vPhase2[iMode][maxLayer] ) * phase1[iMode] ;

			ICoeff[iMode] = - TransDown[iMode] * ( CPLX_ONE - phase2Term[iMode] ) / TxLinePara.s_vCharImpdZ[iMode][maxLayer];

			if( iMode == 0 )
			{
				VI.s_Iv[iMode] = VISrc.s_Vv[iMode] * ICoeff[iMode];
			}
			else
			{
				VCoeff[iMode] = TransDown[iMode] * ( CPLX_ONE + phase2Term[iMode] );
				VI.s_Vi[iMode] = VISrc.s_Vi[iMode] * VCoeff[iMode];
				VI.s_Iv[iMode] = VISrc.s_Vv[iMode] * ICoeff[iMode];
				VI.s_Vv[iMode] = VISrc.s_Vv[iMode] * VCoeff[iMode];
				VI.s_Ii[iMode] = VISrc.s_Vi[iMode] * ICoeff[iMode];
			}
			
			//std::cout << VI.s_Vi[iMode] << TxLinePara.s_vCharImpdZ[iMode][maxLayer] << TxLinePara.s_vRefCoefDown[iMode][maxLayer] << std::endl;
		}
		
		if ( OSPair.s_iObsLayer < OSPair.s_iSrcLayer )
		{
			Tmp = VI.s_Ii[1];
			VI.s_Ii[1] = - VI.s_Vv[1];
			VI.s_Vv[1] = - Tmp;
		}
		//std::cout << VI.s_Iv[0] << VI.s_Ii[1] << VI.s_Iv[1] << VI.s_Vi[1] << VI.s_Vv[1] << std::endl;
	}

	return true;
}

bool CTxLine::calcAsympTxLineVI(in_ const SFormationConfig & ForConfig,	in_ const STxLineDataType & TxLinePara,
								in_ const SObsSrcPair OSPair,			in_ const DataType & krho,
								out_ STxLineVI & AsympVI )
{
	Complex phase1[2];

	double Sign = SIGN( OSPair.s_Obs_z, OSPair.s_Src_zPrim );
	/* If the source and observation points are in the same layer */
	if( OSPair.s_iSrcLayer == OSPair.s_iObsLayer )
	{
		
		for ( unsigned int iMode = 0; iMode < 2; iMode ++ )
		{
			Complex Gamma = CPLX_0;

			if( OSPair.onTopFlag )
			{
				Gamma = ForConfig.s_vStaticFresnelDn[iMode][OSPair.s_iSrcLayer];
			}

			Gamma = ForConfig.s_vStaticFresnelDn[iMode][OSPair.s_iSrcLayer];

			//TM & TE
			AsympVI.s_Iv[iMode] = (0.5/TxLinePara.s_vCharImpdZ[iMode][OSPair.s_iSrcLayer])*( CPLX_ONE - Gamma );

			if(iMode == 1)//TE
			{
				AsympVI.s_Vi[iMode] = (0.5*TxLinePara.s_vCharImpdZ[iMode][OSPair.s_iSrcLayer])*( CPLX_ONE + Gamma );;
				AsympVI.s_Vv[iMode] = 0.5*( Sign - Gamma );
				AsympVI.s_Ii[iMode] = 0.5*( Sign + Gamma );
			}
		}
	}
	else
	{
		std::cout<< "[Warning] : The Asmptotic Transmission Line is NOT Supported for NON-Same Layer!" <<std::endl;
	}
	return true;
}

bool CTxLine::calcExtractTerm(	in_ const SFormationConfig & ForConfig,	in_ const SObsSrcPair & OSPair, in_ const DataType & Omega,
								out_ std::vector<Complex> & rsltVec )
{
	rsltVec.resize(6);
	std::vector<Complex> SRI[2];
	
	double Sign = SIGN( OSPair.s_Obs_z, OSPair.s_Src_zPrim );

	/* calculate the Sommerfeld and Related identities analytically */
	analyticSRITerm( ForConfig.s_vHWaveNum[OSPair.s_iSrcLayer],	ForConfig.s_vTMAnsiRatio[OSPair.s_iSrcLayer], OSPair.s_rho,	SRI[0] ); //TM
	analyticSRITerm( ForConfig.s_vHWaveNum[OSPair.s_iSrcLayer],	ForConfig.s_vTEAnsiRatio[OSPair.s_iSrcLayer], OSPair.s_rho,	SRI[1] ); //TE

	Complex Gamma_TE = CPLX_0;
	Complex Gamma_TM = CPLX_0;

	if( OSPair.onTopFlag )
	{
		Gamma_TM = ForConfig.s_vStaticFresnelDn[0][OSPair.s_iSrcLayer];
		Gamma_TE = ForConfig.s_vStaticFresnelDn[1][OSPair.s_iSrcLayer];
	}

	Gamma_TM = ForConfig.s_vStaticFresnelDn[0][OSPair.s_iSrcLayer];
	Gamma_TE = ForConfig.s_vStaticFresnelDn[1][OSPair.s_iSrcLayer];
	
	//rslt[0]=S0{Iv_e} //TM
	//rslt[1]=S0{Iv_h} //TE
	//rslt[2]=S0{Vi_h*krho*krho}
	//rslt[3]=S1{Ii_h*krho}
	//rslt[4]=S1{Vv_h*krho}
	//rslt[5]={Iv_h,Iv_e}
					
	Complex jOmegaEpst = CPLX_J * Omega * ForConfig.s_vHCompEpsilon[OSPair.s_iSrcLayer];// *EPSILON_0;
	Complex jOmegaMut = CPLX_J * Omega * ForConfig.s_vHCompMu[OSPair.s_iSrcLayer];// *MU_0;
	rsltVec[0] = jOmegaEpst * SRI[0][0] * ( CPLX_ONE - Gamma_TM );
	rsltVec[1] = CPLX_ONE / jOmegaMut * SRI[1][1] * ( CPLX_ONE - Gamma_TE );
	rsltVec[2] = jOmegaMut * (SRI[1][2]) * ( CPLX_ONE + Gamma_TE );
	rsltVec[3] = (SRI[1][3]) * ( Sign + Gamma_TE );
	rsltVec[4] = (SRI[1][3]) * ( Sign - Gamma_TE );
	rsltVec[5] = CPLX_ONE/jOmegaMut * SRI[1][4] * ( CPLX_ONE - Gamma_TE ) - jOmegaEpst * SRI[0][5] * ( CPLX_ONE - Gamma_TM );

	return true;
}

bool CTxLine::analyticSRITerm(	in_ const Complex & wavenumber,	in_ const Complex & nu,	
								in_ const DataType & rho, out_ std::vector<Complex> & SRI )
{
	SRI.resize(6);
	Complex sqrt_nu = sqrt(nu);
	Complex Rp = rho/sqrt_nu;
	Complex OnePlusJKR = CPLX_ONE + CPLX_J * wavenumber * Rp;
	Complex exp_mJKR = std::exp( - CPLX_J * wavenumber * Rp );

	//S0{Iv_e} //TM  //Ia
	SRI[0] = exp_mJKR / 4.0 / PI / Rp / nu;
	
	//S0{Iv_h} //TE  //Ib
	SRI[1] = - SRI[0] * OnePlusJKR/Rp/Rp;
	
	//S0{Vi_h*krho*krho} //TE //Ic
	SRI[2] = (wavenumber * wavenumber * SRI[0] + SRI[1])/nu ;

	//rslt[3]=S1{Ii_h*krho} // TE // =0
	//rslt[4]=S1{Vv_h*krho} // TE // =0
	SRI[3] = CPLX_0 ;
	
	// Id //S1{Iv_h/krho} //TE
	SRI[4] = ( SRI[0] * nu + CPLX_J * wavenumber / 4.0 / PI ) / rho;
	
	// Ie //S1{Iv_e/krho} //TM
	SRI[5] = (CPLX_ONE - exp_mJKR) / 4.0 / PI / CPLX_J / wavenumber / rho;

	//rslt[5]={Iv_h,Iv_e}
	//Intgrand = ( m_vTxLineVIs[iWght].s_Iv[1] - m_vTxLineVIs[iWght].s_Iv[0] );
	//rsltVec[5] = CPLX_ONE/jOmegaMut * SRI[1][4] * ( CPLX_ONE - Gamma_TE ) - jOmegaEpst * SRI[0][5] * ( CPLX_ONE - Gamma_TM );
	
	return true;
}

bool CTxLine::calTLVIatSrcLayer(in_ const SFormationConfig & ForConfig,	in_ const STxLineDataType & TxLinePara,
								in_ const unsigned int & iSrcLayerNum,	in_ const DataType & Obs_z,	in_ const DataType & Src_zPrim, 
								out_ STxLineVI & VI )
{
	Complex phase1[2] = {0.0};
	Complex Dn[2] = {0.0};
	DataType Zeta_n[4] = {0.0};
	Complex Rn[2][4] = {0.0};
	Complex phase[2][4] = {0.0};
	Complex Temp1[2][4] = {0.0};
	double Sgn = SIGN(Obs_z, Src_zPrim);
	
	Zeta_n[0] = Obs_z + Src_zPrim - 2.0 * ForConfig.s_vLayerBoundary[iSrcLayerNum];
	Zeta_n[1] = 2.0 * ForConfig.s_vLayerBoundary[iSrcLayerNum-1] - Obs_z - Src_zPrim;
	Zeta_n[2] = 2.0 * ForConfig.s_vLayerThick[iSrcLayerNum] + ( Obs_z - Src_zPrim );
	Zeta_n[3] = 2.0 * ForConfig.s_vLayerThick[iSrcLayerNum] - ( Obs_z - Src_zPrim );

	for ( unsigned int iMode = 0; iMode < 2; iMode ++ )
	{
		phase1[iMode] = std::exp( - CPLX_J * TxLinePara.s_vWaveNumber_z[iMode][iSrcLayerNum] * abs( Obs_z - Src_zPrim ));

		Rn[iMode][0] = TxLinePara.s_vRefCoefDown[iMode][iSrcLayerNum];
		Rn[iMode][1] = TxLinePara.s_vRefCoefUp[iMode][iSrcLayerNum];
		Rn[iMode][2] = Rn[iMode][1] * Rn[iMode][0];
		Rn[iMode][3] = Rn[iMode][2];

		for( unsigned int index = 0; index < 4; index++ )
		{
			phase[iMode][index] = exp( - CPLX_J * TxLinePara.s_vWaveNumber_z[iMode][iSrcLayerNum]*Zeta_n[index] );
			Temp1[iMode][index] = Rn[iMode][index] * phase[iMode][index];
		}

		Dn[iMode] = CPLX_ONE - Rn[iMode][0] * Rn[iMode][1] * TxLinePara.s_vPhase2[iMode][iSrcLayerNum];
	}

	VI.s_Iv[0] = DataType(0.5) / TxLinePara.s_vCharImpdZ[0][iSrcLayerNum] 
			   * ( phase1[0] + CPLX_ONE/Dn[0] * ( - Temp1[0][0] - Temp1[0][1] + Temp1[0][2] + Temp1[0][3]) );
	
	VI.s_Iv[1] = DataType(0.5) / TxLinePara.s_vCharImpdZ[1][iSrcLayerNum] 
			   * ( phase1[1] + CPLX_ONE/Dn[1] * ( - Temp1[1][0] - Temp1[1][1] + Temp1[1][2] + Temp1[1][3]) );

	VI.s_Vi[1] = DataType(0.5) * TxLinePara.s_vCharImpdZ[1][iSrcLayerNum] 
			   * ( phase1[1] + CPLX_ONE/Dn[1] * (   Temp1[1][0] + Temp1[1][1] + Temp1[1][2] + Temp1[1][3]) );

	VI.s_Vv[0] = DataType(0.5) * ( Sgn * phase1[0] - CPLX_ONE / Dn[0] * (   Temp1[0][0] - Temp1[0][1] - Temp1[0][2] + Temp1[0][3]) );
	VI.s_Vv[1] = DataType(0.5) * ( Sgn * phase1[1] - CPLX_ONE / Dn[1] * (   Temp1[1][0] - Temp1[1][1] - Temp1[1][2] + Temp1[1][3]) );
	VI.s_Ii[1] = DataType(0.5) * ( Sgn * phase1[1] - CPLX_ONE / Dn[1] * ( - Temp1[1][0] + Temp1[1][1] - Temp1[1][2] + Temp1[1][3]) );

	return true;
}