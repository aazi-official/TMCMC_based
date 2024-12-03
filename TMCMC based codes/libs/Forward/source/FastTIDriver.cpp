#include "FastTIDriver.h"

CFastTIDriver::CFastTIDriver( )
{
}

CFastTIDriver::~CFastTIDriver( )
{
}

bool CFastTIDriver::initToolConfig( in_ const std::string & szToolFile )
{
	m_fTOOLIn.open( szToolFile );

	if ( !m_fTOOLIn )
	{
		std::cout<< "Unable to Open TOOL configuration File: " << szToolFile << std::endl;
		return false;
	}

	m_fTOOLIn >> m_dFreq;
	m_dOmega = m_dFreq * TWOPI;
	const unsigned int MAX=160;
	char szComment[MAX]; 
	m_fTOOLIn.getline(szComment,MAX); 
	m_fTOOLIn.getline(szComment,MAX); 

	// Number of Transmitter and Receivers
	m_fTOOLIn >> m_sToolConf.s_iNum_Trans >> m_sToolConf.s_iNum_Recvs;

	m_sToolConf.s_iNum_TransTurns.resize(m_sToolConf.s_iNum_Trans);
	m_sToolConf.s_vDistTrans.resize(m_sToolConf.s_iNum_Trans);
	m_sToolConf.s_iNum_RecvsTurns.resize(m_sToolConf.s_iNum_Recvs);
	m_sToolConf.s_vDistRecvs.resize(m_sToolConf.s_iNum_Recvs);

	// Input the unit for tool
	std::string szUnit;
	m_fTOOLIn >> szUnit;
	convertUnit( szUnit, m_sToolConf.s_FACTOR_T );

	// Input the transmitters and receivers info
	for( unsigned int i = 0; i < m_sToolConf.s_iNum_Trans; i++ )
	{
		m_fTOOLIn >> m_sToolConf.s_iNum_TransTurns[i] >> m_sToolConf.s_vDistTrans[i];
		m_sToolConf.s_vDistTrans[i] = - m_sToolConf.s_vDistTrans[i] * m_sToolConf.s_FACTOR_T;
	}

	for( unsigned int i = 0; i < m_sToolConf.s_iNum_Recvs; i++ )
	{
		m_fTOOLIn >> m_sToolConf.s_iNum_RecvsTurns[i] >> m_sToolConf.s_vDistRecvs[i];
		m_sToolConf.s_vDistRecvs[i] = - m_sToolConf.s_vDistRecvs[i] * m_sToolConf.s_FACTOR_T;
	}

	DataType TRdist = 0;
	for( unsigned int itran = 0; itran < m_sToolConf.s_iNum_Trans; itran++ )
	{
		for( unsigned int jrecv = 0; jrecv < m_sToolConf.s_iNum_Recvs; jrecv++ )
		{
			TRdist = abs( m_sToolConf.s_vDistRecvs[jrecv] - m_sToolConf.s_vDistTrans[itran] );
			m_sToolConf.s_dMaxDist = max( TRdist, m_sToolConf.s_dMaxDist );
		}
	}

	m_fTOOLIn.close();
	return true;
}

bool CFastTIDriver::initFormationPara( in_ const std::string & szFormationFile )
{
	m_fFormIn.open( szFormationFile );

	if ( !m_fFormIn )
	{
		std::cout<< "Unable to Open TOOL configuration File: " << szFormationFile << std::endl;
		return false;
	}

	DataType Alfa, Beta, Gama;
	m_fFormIn >> Alfa >> Beta >> Gama;

	// determine the well type

	Alfa = fmod(Alfa, 360.0);

	if( Alfa == 0.0 || Alfa == 180.0)
	{
		m_WellType = ENUM_VERTICAL;
	}
	else if (Alfa == 90.0 || Alfa == 270.0)
	{
		m_WellType = ENUM_HORIZONTAL;
	}
	else
	{
		m_WellType = ENUM_DEVIATED;

		//if( Alfa > 90.0 )
		//{
		//	Alfa = Alfa - 90.0;
		//}
	}

	// convert from degree to radius
	Alfa = Alfa * DEG_TO_RAD;
	Beta = Beta * DEG_TO_RAD;
	Gama = Gama * DEG_TO_RAD;
	// calculate rotation matrix R
	m_sFormConf.s_dyRotMatrixR.xx =  cos(Alfa)*cos(Beta)*cos(Gama) - sin(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixR.xy = -cos(Alfa)*cos(Beta)*sin(Gama) - sin(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixR.xz =  sin(Alfa)*cos(Beta);
	m_sFormConf.s_dyRotMatrixR.yx =  cos(Alfa)*sin(Beta)*cos(Gama) + cos(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixR.yy = -cos(Alfa)*sin(Beta)*sin(Gama) + cos(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixR.yz =  sin(Alfa)*sin(Beta);
	m_sFormConf.s_dyRotMatrixR.zx = -sin(Alfa)*cos(Gama);
	m_sFormConf.s_dyRotMatrixR.zy =  sin(Alfa)*sin(Gama);
	m_sFormConf.s_dyRotMatrixR.zz =  cos(Alfa);

	// calculate inverse rotation matrix R
	m_sFormConf.s_dyRotMatrixInvR.xx = cos(Alfa)*cos(Beta)*cos(Gama)-sin(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixInvR.xy = cos(Alfa)*sin(Beta)*cos(Gama)+cos(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixInvR.xz =-sin(Alfa)*cos(Gama);
	m_sFormConf.s_dyRotMatrixInvR.yx =-cos(Alfa)*cos(Beta)*sin(Gama)-sin(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixInvR.yy =-cos(Alfa)*sin(Beta)*sin(Gama)+cos(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixInvR.yz = sin(Alfa)*sin(Gama);
	m_sFormConf.s_dyRotMatrixInvR.zx = sin(Alfa)*cos(Beta);
	m_sFormConf.s_dyRotMatrixInvR.zy = sin(Alfa)*sin(Beta);
	m_sFormConf.s_dyRotMatrixInvR.zz = cos(Alfa);

	m_sFormConf.s_dAlfa = Alfa;
	m_sFormConf.s_dBeta = Beta;
	m_sFormConf.s_dGama = Gama;

	// Input the unit for formation
	std::string szUnit;
	m_fFormIn >> szUnit;
	m_fFormIn >> szUnit;
	convertUnit( szUnit, m_sFormConf.s_FACTOR_T );

	// 
	int sign;
	m_fFormIn >> sign;
	m_sFormConf.s_bSign = ( sign == 0 ? false : true );

	const unsigned int MAX=160;
	char szComment[MAX]; 

	m_fFormIn.getline(szComment,MAX); 

	m_fFormIn >> m_sFormConf.s_iNum_Layers;
	m_sFormConf.s_vHResistivity.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVResistivity.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHRelEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVRelEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHRelMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVRelMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHCompEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVCompEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHCompMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVCompMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vjOmegaVMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vTMAnsiRatio.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vTEAnsiRatio.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHWaveNum.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vStaticFresnelDn[0].resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vStaticFresnelDn[1].resize(m_sFormConf.s_iNum_Layers + 2);

	m_sFormConf.s_K0 = m_dOmega / VLIGHT_0; 

	// Need Check ?????
	m_sFormConf.s_vLayerBoundary.resize(m_sFormConf.s_iNum_Layers+1);
	m_sFormConf.s_vLayerThick.resize(m_sFormConf.s_iNum_Layers+1);

	for( unsigned int i = 1; i <= m_sFormConf.s_iNum_Layers; i++ )
	{
		if( i != m_sFormConf.s_iNum_Layers )
		{
			m_fFormIn >> m_sFormConf.s_vLayerBoundary[i] >> m_sFormConf.s_vHResistivity[i] >> m_sFormConf.s_vVResistivity[i]
					  >> m_sFormConf.s_vHRelEpsilon[i] >> m_sFormConf.s_vVRelEpsilon[i] >> m_sFormConf.s_vHRelMu[i];

			m_sFormConf.s_vVRelMu[i] = m_sFormConf.s_vHRelMu[i];
			// Need Check ?????
			m_sFormConf.s_vLayerBoundary[i] = - m_sFormConf.s_vLayerBoundary[i];
		}
		else
		{
			m_fFormIn >> m_sFormConf.s_vHResistivity[i] >> m_sFormConf.s_vVResistivity[i]
					  >> m_sFormConf.s_vHRelEpsilon[i] >> m_sFormConf.s_vVRelEpsilon[i] >> m_sFormConf.s_vHRelMu[i];
		}

		m_sFormConf.s_vHCompEpsilon[i] = m_sFormConf.s_vHRelEpsilon[i] * EPSILON_0 - CPLX_J * (CPLX_ONE/m_sFormConf.s_vHResistivity[i] ) / m_dOmega;
		m_sFormConf.s_vVCompEpsilon[i] = m_sFormConf.s_vVRelEpsilon[i] * EPSILON_0 - CPLX_J * (CPLX_ONE/m_sFormConf.s_vVResistivity[i] ) / m_dOmega;
	
		m_sFormConf.s_vHCompMu[i] = m_sFormConf.s_vHRelMu[i] * MU_0;
		m_sFormConf.s_vVRelMu[i] = m_sFormConf.s_vHRelMu[i];
		m_sFormConf.s_vVCompMu[i] = m_sFormConf.s_vVRelMu[i] * MU_0;
		m_sFormConf.s_vjOmegaVMu[i] = CPLX_J * m_dOmega * m_sFormConf.s_vVCompMu[i];

		m_sFormConf.s_vTMAnsiRatio[i] = m_sFormConf.s_vHCompEpsilon[i]/m_sFormConf.s_vVCompEpsilon[i];
		m_sFormConf.s_vTEAnsiRatio[i] = m_sFormConf.s_vHCompMu[i]/m_sFormConf.s_vVCompMu[i];
		m_sFormConf.s_vHWaveNum[i] = m_dOmega * sqrt( m_sFormConf.s_vHCompEpsilon[i] * m_sFormConf.s_vHCompMu[i] );

		m_sFormConf.s_vLayerBoundary[i] = m_sFormConf.s_vLayerBoundary[i] * m_sFormConf.s_FACTOR_T * cos(m_sFormConf.s_bSign*m_sFormConf.s_dAlfa);
		
	}

	m_sFormConf.s_vHResistivity[0] = m_sFormConf.s_vHResistivity[1];
	m_sFormConf.s_vHResistivity[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHResistivity[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVResistivity[0] = m_sFormConf.s_vVResistivity[1];
	m_sFormConf.s_vVResistivity[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVResistivity[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHRelEpsilon[0] = m_sFormConf.s_vHRelEpsilon[1];
	m_sFormConf.s_vHRelEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHRelEpsilon[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVRelEpsilon[0] = m_sFormConf.s_vVRelEpsilon[1];
	m_sFormConf.s_vVRelEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVRelEpsilon[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHRelMu[0] = m_sFormConf.s_vHRelMu[1];
	m_sFormConf.s_vHRelMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHRelMu[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVRelMu[0] = m_sFormConf.s_vVRelMu[1];
	m_sFormConf.s_vVRelMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVRelMu[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHCompEpsilon[0] = m_sFormConf.s_vHCompEpsilon[1];
	m_sFormConf.s_vHCompEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHCompEpsilon[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVCompEpsilon[0] = m_sFormConf.s_vVCompEpsilon[1];
	m_sFormConf.s_vVCompEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVCompEpsilon[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHCompMu[0] = m_sFormConf.s_vHCompMu[1];
	m_sFormConf.s_vHCompMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHCompMu[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVCompMu[0] = m_sFormConf.s_vVCompMu[1];
	m_sFormConf.s_vVCompMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVCompMu[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vjOmegaVMu[0] = m_sFormConf.s_vjOmegaVMu[1];
	m_sFormConf.s_vjOmegaVMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vjOmegaVMu[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vTMAnsiRatio[0] = m_sFormConf.s_vTMAnsiRatio[1];
	m_sFormConf.s_vTMAnsiRatio[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vTMAnsiRatio[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vTEAnsiRatio[0] = m_sFormConf.s_vTEAnsiRatio[1];
	m_sFormConf.s_vTEAnsiRatio[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vTEAnsiRatio[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHWaveNum[0] = m_sFormConf.s_vHWaveNum[1];
	m_sFormConf.s_vHWaveNum[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHWaveNum[m_sFormConf.s_iNum_Layers];

	for( unsigned int i = 1; i <= m_sFormConf.s_iNum_Layers; i++ )
	{
		// for TM mode;
		//DataType TM_Tmp_L = sqrt( m_sFormConf.s_vHRelEpsilon[i] * m_sFormConf.s_vVRelEpsilon[i] );
		//DataType TM_Tmp_LP1 = sqrt( m_sFormConf.s_vHRelEpsilon[i+1] * m_sFormConf.s_vVRelEpsilon[i+1] );
		
		Complex TM_Tmp_L = sqrt(m_sFormConf.s_vHCompEpsilon[i] * m_sFormConf.s_vVCompEpsilon[i]);
		Complex TM_Tmp_LP1 = sqrt(m_sFormConf.s_vHCompEpsilon[i + 1] * m_sFormConf.s_vVCompEpsilon[i + 1]);
		m_sFormConf.s_vStaticFresnelDn[0][i] = (TM_Tmp_L - TM_Tmp_LP1) / (TM_Tmp_L + TM_Tmp_LP1);

		// for TE mode;
		//DataType TE_Tmp_L = sqrt( m_sFormConf.s_vHRelMu[i] * m_sFormConf.s_vVRelMu[i] );
		//DataType TE_Tmp_LP1 = sqrt( m_sFormConf.s_vHRelMu[i+1] * m_sFormConf.s_vVRelMu[i+1] );
		Complex TE_Tmp_L = sqrt(m_sFormConf.s_vHCompMu[i] * m_sFormConf.s_vVCompMu[i]);
		Complex TE_Tmp_LP1 = sqrt(m_sFormConf.s_vHCompMu[i + 1] * m_sFormConf.s_vVCompMu[i + 1]);
		m_sFormConf.s_vStaticFresnelDn[1][i] = (TE_Tmp_LP1 - TE_Tmp_L) / (TE_Tmp_L + TE_Tmp_LP1);
	}

	//m_sFormConf.s_vLayerBoundary[0] = m_sFormConf.s_vLayerBoundary[1];
	//m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers] = m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers-1];

	m_fFormIn >> m_dDepthSt >> m_dStepSz >> m_dDepthEd;
	
	m_NumLgPts = (m_dDepthEd - m_dDepthSt)/m_dStepSz + 1;

	m_dDepthSt = - m_dDepthSt * m_sFormConf.s_FACTOR_T;
	m_dStepSz = - m_dStepSz * m_sFormConf.s_FACTOR_T;
	m_dDepthEd = - m_dDepthEd * m_sFormConf.s_FACTOR_T;

	DataType maxLogz = max(m_dDepthSt, m_dDepthEd);
	DataType minLogz = min(m_dDepthSt, m_dDepthEd);
	DataType buff = 10.0;

	if( maxLogz < m_sFormConf.s_vLayerBoundary[1] )
		m_sFormConf.s_vLayerBoundary[0] = m_sFormConf.s_vLayerBoundary[1] + buff;
	else
		m_sFormConf.s_vLayerBoundary[0] = maxLogz + m_sToolConf.s_dMaxDist + buff;

	// [D.W.] Modify to calculate s_vLayerBoundary[m_sFormConf.s_iNum_Layers] correctly

	if( minLogz > m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers-1] )
		m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers] = m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers-1] - buff;
	else
		m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers] = minLogz - m_sToolConf.s_dMaxDist - buff;
	
	for( unsigned int i = 1; i <= m_sFormConf.s_iNum_Layers; i++ )
		m_sFormConf.s_vLayerThick[i] = m_sFormConf.s_vLayerBoundary[i-1] - m_sFormConf.s_vLayerBoundary[i];

	m_fFormIn.close();
	return true;
}

bool CFastTIDriver::calcLogs()
{
	SObsSrcPair OSPair;
	Complex rslt[6] = {CPLX_0};
	m_vTxLineParas.resize(HTRANS_PTS);

	for( unsigned int i = 0; i < HTRANS_PTS; i++ )
	{
		m_vTxLineParas[i].Init(m_sFormConf.s_iNum_Layers);
	}

	m_vTxLineVIs.resize(HTRANS_PTS);

	m_vHfieldDyd.resize(m_sToolConf.s_iNum_Trans*m_sToolConf.s_iNum_Recvs*m_NumLgPts);

	//omp_init_lock( &m_TxParalock );
	//omp_init_lock( &m_TxVIlock );

	//
	
	//
	std::complex<DataType> S2;
	Complex jOmegaMuz = CPLX_0;
	Complex jOmegaMuzPrim = CPLX_0;
	unsigned int iCnt = 0;
	for( unsigned int iTran = 0; iTran < m_sToolConf.s_iNum_Trans; iTran++ )
	{
		for( unsigned int iRecv = 0; iRecv < m_sToolConf.s_iNum_Recvs; iRecv++ )
		{
			Complex Norm_Factor = m_sToolConf.s_iNum_TransTurns[iTran] * m_sToolConf.s_iNum_RecvsTurns[iRecv] * CPLX_J * m_dOmega * MU_0;
			//Complex Norm_Factor = CPLX_J * m_dOmega * MU_0;
			// for each transmitter and reciever pair, rho remains the same

			// [D.W.] Modify to calculate s_dBeta correctly

			DataType Numerator = 0.0;
			DataType Denominator = (m_sToolConf.s_vDistTrans[iTran] - m_sToolConf.s_vDistRecvs[iRecv])* sin(m_sFormConf.s_dAlfa);
			m_sFormConf.s_dBeta = atan2(Numerator, Denominator);

			DataType cosBeta = cos(m_sFormConf.s_dBeta);
			DataType sinBeta = sin(m_sFormConf.s_dBeta);
			DataType cos2Beta = cos(2 * m_sFormConf.s_dBeta);
			DataType sin2Beta = sin(2 * m_sFormConf.s_dBeta);

			OSPair.s_rho = abs( (m_sToolConf.s_vDistTrans[iTran] - m_sToolConf.s_vDistRecvs[iRecv])* sin( m_sFormConf.s_dAlfa ) );

			m_iHanksMin = 60;
			m_iHanksMax = 180;
			m_vHanksMin.resize(HTRANS_NUM);
			m_vHanksMax.resize(HTRANS_NUM);

			m_iHanksMin = 0;
			m_iHanksMax = 282;

			calcAllocKrhoVec(m_iHanksMin, m_iHanksMax, OSPair.s_rho);

			for( unsigned int iLgPt = 1; iLgPt <= m_NumLgPts; iLgPt++ )
			{
				OSPair.s_Src_zPrim = ( m_dDepthSt + ( iLgPt - 1 ) * m_dStepSz ) * cos( m_sFormConf.s_bSign * m_sFormConf.s_dAlfa ) 
								   + m_sToolConf.s_vDistTrans[iTran] * cos( m_sFormConf.s_dAlfa );

				bool onTopFlag = false;

				locateLayerNum( m_sFormConf, OSPair.s_Src_zPrim, OSPair.s_iSrcLayer, onTopFlag);
				if (m_WellType == ENUM_HORIZONTAL)
				{
					OSPair.s_Obs_z = OSPair.s_Src_zPrim;
					OSPair.s_iObsLayer = OSPair.s_iSrcLayer;
					OSPair.onTopFlag = onTopFlag;
				}
				else
				{
					OSPair.s_Obs_z = (m_dDepthSt + (iLgPt - 1) * m_dStepSz) * cos(m_sFormConf.s_bSign * m_sFormConf.s_dAlfa)
						+ m_sToolConf.s_vDistRecvs[iRecv] * cos(m_sFormConf.s_dAlfa);
					locateLayerNum(m_sFormConf, OSPair.s_Obs_z, OSPair.s_iObsLayer, onTopFlag);
				}
				
				if( m_WellType == ENUM_VERTICAL )
				{
					// TO DO: Handle vertical well;
					CQuadRule DERule(m_TxLineObj, m_sFormConf, m_dOmega, 1E-6);
					CQuadRule::iCount = 0;
					Complex s1[2],s2[2];
					DataType a = 1.e-6*m_sFormConf.s_K0;
					DataType b = m_sFormConf.s_K0;
					//DERule.tanhsinhQuad( a, b, OSPair, s1 );
					DERule.mixedQuad( a, OSPair, s2 );
					rslt[0] = (s1[0] + s2[0])/TWOPI*Norm_Factor;
					rslt[2] = (s1[1] + s2[1])/TWOPI*Norm_Factor;
					S2 = std::complex<DataType>(0.0, 0.0);
					rslt[0].imag(-rslt[0].imag());
					rslt[2].imag(-rslt[2].imag());
					//std::cout << "DE Quadarture Rule: " << CQuadRule::iCount << "Points Employed!" << std::endl;

				}
				else if ( m_WellType == ENUM_DEVIATED || m_WellType == ENUM_HORIZONTAL )
				{
					precalcTxLineVIs(m_iHanksMin, m_iHanksMax, OSPair);
					//#pragma omp parallel num_threads(8)
					{
						//#pragma omp for
#pragma omp parallel for
						for( int iHXform = 0; iHXform <= 5; iHXform++ )
						{
							calcHANKTrans( iHXform, OSPair, rslt[iHXform] );
						
							rslt[iHXform] = rslt[iHXform]/TWOPI*Norm_Factor; 
							rslt[iHXform].imag( -rslt[iHXform].imag() );
						}
					}

					// If horizontal well is considered, add back the extracted terms.
					if( m_WellType == ENUM_HORIZONTAL )
					{
						std::vector<Complex> ExtTerms;
						m_TxLineObj.calcExtractTerm( m_sFormConf, OSPair, m_dOmega, ExtTerms );
						for( int iHXform = 0; iHXform <= 5; iHXform++ )
						{
							ExtTerms[iHXform] = ExtTerms[iHXform] * Norm_Factor;
							ExtTerms[iHXform].imag(-ExtTerms[iHXform].imag());
							rslt[iHXform] = rslt[iHXform] + ExtTerms[iHXform];
							//rslt[iHXform].imag(-rslt[iHXform].imag());
						}
					}
					S2 = DataType(2.0) * rslt[5] / OSPair.s_rho - (rslt[1] - rslt[0]);
				}
				else
				{
					std::cout<< "[Warning]: Well TYPE is NOT Supported Yet!" << std::endl;
				}
				//std::cout<< "Min Spec: " << m_iHanksMin << ", Max Spec: " << m_iHanksMax << std::endl;

				m_iHanksMin = *std::min_element(m_vHanksMin.begin(), m_vHanksMin.end() );
				m_iHanksMax = *std::max_element(m_vHanksMax.begin(), m_vHanksMax.end() );

				jOmegaMuz = m_sFormConf.s_vjOmegaVMu[OSPair.s_iObsLayer];
				jOmegaMuzPrim = m_sFormConf.s_vjOmegaVMu[OSPair.s_iSrcLayer];
				// save data to dyadic, Note the convention here, H_ij is the H_j field due to i-comp M source 
				
				//
				//rslt[0]=S0{Iv_e}
				//rslt[1]=S0{Iv_h}
				//rslt[2]=S0{Vi_h*krho*krho}
				//rslt[3]=S1{Ii_h*krho}
				//rslt[4]=S1{Vv_h*krho}
				//rslt[5]={Iv_h,Iv_e}
				//S2=I5=S2{Iv_h-Iv_e};
				//
				
				m_vHfieldDyd[iCnt].xx = -DataType(0.5) * (rslt[0] + rslt[1] - cos2Beta * S2);
				m_vHfieldDyd[iCnt].xy = DataType(0.5) * (sin2Beta * S2);
				m_vHfieldDyd[iCnt].zx = cosBeta * rslt[3] / jOmegaMuzPrim;
				m_vHfieldDyd[iCnt].yx = m_vHfieldDyd[iCnt].xy;
				m_vHfieldDyd[iCnt].yy = -DataType(0.5) * (rslt[0] + rslt[1] + cos2Beta * S2);
				m_vHfieldDyd[iCnt].zy = sinBeta * rslt[3] / jOmegaMuzPrim;
				m_vHfieldDyd[iCnt].xz = cosBeta * rslt[4] / jOmegaMuz;
				m_vHfieldDyd[iCnt].yz = sinBeta * rslt[4] / jOmegaMuz;
				m_vHfieldDyd[iCnt].zz = rslt[2] / jOmegaMuzPrim / jOmegaMuz;
				
				// Apply rotation matrix operation
				m_vHfieldDyd[iCnt] = m_sFormConf.s_dyRotMatrixInvR * ( m_vHfieldDyd[iCnt] * m_sFormConf.s_dyRotMatrixR );

				iCnt++;
				m_VIsFlg.reset();
			}
			m_SpcFlg.reset();
		}
	}
	
	//omp_destroy_lock( &m_TxParalock );
	//omp_destroy_lock( &m_TxVIlock );

	return true;
}

void CFastTIDriver::calcHANKTrans( in_ const unsigned int iHTrans, in_ const SObsSrcPair OSPair, out_ Complex & rslt )
{
	//Complex CMAX1 = CPLX_0;
	rslt = CPLX_0;
	Complex T = CPLX_0;
	int NONE = 0;
	Complex TMAX = CPLX_0;

	double Y1 = HANKX_Y1/double(OSPair.s_rho);
	unsigned int I = 130;
	double Y = Y1 * HANKX_E;
	unsigned int M = 110;
	goto L200;
	
L110:

	TMAX.real(std::max( abs( T.real() ), TMAX.real() ));
	TMAX.imag(std::max( abs( T.imag() ), TMAX.imag() ));
	I++;
	Y *= HANKX_E;

	if( I > 148 )
	{
		if( TMAX.real() == 0.0 && TMAX.imag() == 0.0 )
		{
			NONE = 1;
		}
		TMAX =TOL*TMAX;
		M = 120;
	}
	goto L200;

L120:

	if( abs(T.real()) <= TMAX.real() && abs(T.imag()) <= TMAX.imag() )
	{
		//std::cout<< I << std::endl;
		m_vHanksMax[iHTrans] = I;
		Y = Y1;
		M = 140;
		I = 129;
	}
	else
	{
		I++;
		Y *= HANKX_E;
		if( I > 282 )
		{
			m_vHanksMax[iHTrans] = I;
			Y = Y1;
			M = 140;
			I = 129;
		}
	}
	goto L200;

L140 :

	if( abs(T.real()) <= TMAX.real() && abs(T.imag()) <= TMAX.imag() && NONE == 0  )
	{
		rslt = rslt/OSPair.s_rho;
		//std::cout<< I << std::endl;
		m_vHanksMin[iHTrans] = I;
		return;
	}
	else
	{
		I--;
		Y = Y*HANKX_ER;
		if( I <= 0 )
		{
			rslt = rslt/OSPair.s_rho;
			//std::cout<< I << std::endl;
			m_vHanksMin[iHTrans] = I;
			return;
		}
	}

L200:
	DataType G = DataType(Y);

    getIntgrand( iHTrans, I, G, OSPair, T );
	if( iHTrans < 3 )  // S0:J0
	{
		T=T*HTRANS_W0[I];
	}
	else			   // S1:J1	
	{
		T=T*HTRANS_W1[I];
	}
		
    rslt=rslt+T;
    switch( M )
	{
	case 110:
		goto L110;
		break;
	case 120:
		goto L120;
		break;
	case 140:
		goto L140;
		break;
	default:
		break;
	}
	return;
}

void CFastTIDriver::getIntgrand( in_ const unsigned int iHTrans, in_ const unsigned int iWght, in_ const DataType krho, in_ const SObsSrcPair OSPair, out_ Complex & Intgrand )
{
	if( ! m_SpcFlg.test(iWght) )
	{
		STxLineDataType TxLinePara;
		for( unsigned int iMode = 0; iMode < 2; iMode ++ )
		{
			TxLinePara.Init(m_sFormConf.s_iNum_Layers);
		}
		m_TxLineObj.calcTxLineParas( m_sFormConf, m_dOmega,	krho, TxLinePara );
		//while ( ! omp_test_lock(&m_TxParalock) )
		//{
		//}
		m_vTxLineParas[iWght] = TxLinePara;
		m_SpcFlg.set(iWght);
		//omp_unset_lock(&m_TxParalock);
	}
	
	if( ! m_VIsFlg.test(iWght) )
	{	
		STxLineVI TxLineVI;
		m_TxLineObj.calcTxLineVIs( m_sFormConf, m_vTxLineParas[iWght], OSPair, krho, TxLineVI );

		if( m_WellType == ENUM_HORIZONTAL )
		{
			STxLineVI AsympTxLineVI;
			m_TxLineObj.calcAsympTxLineVI( m_sFormConf, m_vTxLineParas[iWght], OSPair, krho, AsympTxLineVI );
			TxLineVI -= AsympTxLineVI;
		}
		
		//while ( ! omp_test_lock(&m_TxVIlock) )
		//{	
		//}
		m_vTxLineVIs[iWght] = TxLineVI;
		m_VIsFlg.set(iWght);
		//omp_unset_lock(&m_TxVIlock);
	}
	
	
	
	//
	//rslt[0]=S0{Iv_e}
	//rslt[1]=S0{Iv_h}
	//rslt[2]=S0{Vi_h*krho*krho}
	//rslt[3]=S1{Ii_h*krho}
	//rslt[4]=S1{Vv_h*krho}
	//rslt[5]={Iv_h,Iv_e}
	//S2=I5=S2{Iv_h-Iv_e};
	//
				
	switch(iHTrans)
	{
	case 0:
		Intgrand = m_vTxLineVIs[iWght].s_Iv[0]*krho;
		break;
	case 1:
		Intgrand = m_vTxLineVIs[iWght].s_Iv[1]*krho;
		break;
	case 2:
		Intgrand = m_vTxLineVIs[iWght].s_Vi[1]*krho*krho*krho;
		break;
	case 3:
		Intgrand = m_vTxLineVIs[iWght].s_Ii[1]*krho*krho;
		break;
	case 4:
		Intgrand = m_vTxLineVIs[iWght].s_Vv[1]*krho*krho;
		break;
	case 5:
		Intgrand = ( m_vTxLineVIs[iWght].s_Iv[1] - m_vTxLineVIs[iWght].s_Iv[0] );
		break;
	default:
		break;
	}
	return;
}

void CFastTIDriver::calcAllocKrhoVec(in_ const unsigned int iStart, in_ const unsigned int iEnd, in_ const DataType & rho)
{
	m_vSpckrho.resize(HTRANS_PTS);

	double Y1 = HANKX_Y1 / double(rho);

	unsigned int I = 129;
	m_vSpckrho[I] = Y1;

	for (int I = 129; I < HTRANS_PTS - 1; I++)
	{
		m_vSpckrho[I + 1] = m_vSpckrho[I] * HANKX_E;
	}
	for (int I = 129; I > 0; I--)
	{
		m_vSpckrho[I - 1] = m_vSpckrho[I] * HANKX_ER;
	}

//#pragma omp parallel num_threads(8)
	{
//#pragma omp for
#pragma omp parallel for
		for (int iWght = iStart; iWght <= iEnd; iWght++)
		{
			m_TxLineObj.calcTxLineParas(m_sFormConf, m_dOmega, m_vSpckrho[iWght], m_vTxLineParas[iWght]);
			m_SpcFlg.set(iWght);
		}
	}
}

void CFastTIDriver::precalcTxLineVIs( in_ const unsigned int iStart, in_ const unsigned int iEnd, in_ const SObsSrcPair OSPair )
{

//#pragma omp parallel num_threads(8)
	{
#pragma omp parallel for
		//#pragma omp for
		for (int iWght = iStart; iWght < iEnd; iWght++)
		//for (int iWght = iEnd - 1; iWght >= iStart; iWght--)
		{
			//m_TxLineObj.calcTxLineParas(m_sFormConf, m_dOmega, m_vSpckrho[iWght], m_vTxLineParas[iWght]);
			//m_SpcFlg.set(iWght);

			STxLineVI TxLineVI;
			m_TxLineObj.calcTxLineVIs(m_sFormConf, m_vTxLineParas[iWght], OSPair, m_vSpckrho[iWght], TxLineVI);

			if (m_WellType == ENUM_HORIZONTAL)
			{
				STxLineVI AsympTxLineVI;
				m_TxLineObj.calcAsympTxLineVI(m_sFormConf, m_vTxLineParas[iWght], OSPair, m_vSpckrho[iWght], AsympTxLineVI);
				TxLineVI -= AsympTxLineVI;
			}
			m_vTxLineVIs[iWght] = TxLineVI;
			m_VIsFlg.set(iWght);
		}
	}
	
}

void CFastTIDriver::locateLayerNum( in_ const SFormationConfig & FormConf, in_ const DataType & z, out_ unsigned int & iLayer, out_ bool & onTopFlag  )
{
	onTopFlag = false;

	//if( z >= FormConf.s_vLayerBoundary[1] )
	//{
	//	iLayer = 1;
	//	if( z == FormConf.s_vLayerBoundary[1] )
	//	{
	//		onTopFlag == true;
	//	}
	//}
	//else if( z < FormConf.s_vLayerBoundary[FormConf.s_iNum_Layers-1] )
	//{
	//	iLayer = FormConf.s_iNum_Layers;
	//}
	//else
	//{
	//	for( unsigned int i = 2; i < FormConf.s_iNum_Layers; i++ )
	//	{
	//		if( z >= FormConf.s_vLayerBoundary[i] && z < FormConf.s_vLayerBoundary[i-1] )
	//		{
	//			iLayer = i;
	//			if( z == FormConf.s_vLayerBoundary[i] )
	//			{
	//				onTopFlag == true;
	//			}
	//			break;
	//		}
	//	}
	//}

	//int iLayerID = 0;

	//for (int i = 0; i < FormConf.s_iNum_Layers; i++)
	//{
	//	if (z == FormConf.s_vLayerBoundary[i])
	//	{
	//		onTopFlag == true;
	//	}

	//	if (z>FormConf.s_vLayerBoundary[i])
	//	{
	//		break;
	//	}

	//	iLayerID++;
	//}

	//iLayer = iLayerID;

	// [D.W.] Correct Layer Number Calculation

	if (z > FormConf.s_vLayerBoundary[1])
	{
		iLayer = 1;
	}
	else if (z <= FormConf.s_vLayerBoundary[FormConf.s_iNum_Layers - 1])
	{
		iLayer = FormConf.s_iNum_Layers;
		if (z == FormConf.s_vLayerBoundary[FormConf.s_iNum_Layers - 1])
		{
			onTopFlag == true;
		}

	}
	else
	{
		for (unsigned int i = 2; i < FormConf.s_iNum_Layers; i++)
		{
			if (z > FormConf.s_vLayerBoundary[i] && z <= FormConf.s_vLayerBoundary[i - 1])
			{
				iLayer = i;
				if (z == FormConf.s_vLayerBoundary[i-1])
				{
					onTopFlag == true;
				}
				break;
			}
		}
	}

	return;
}

void CFastTIDriver::convertUnit( in_ const string & szUnit, out_ DataType & FACTOR_T )
{
	if ( szUnit == "M" || szUnit == "m" )
	{
		FACTOR_T = M_TO_METER;
	}
	else if ( szUnit == "DM" || szUnit == "dm" )
	{
		FACTOR_T = DM_TO_METER;
	}
	else if ( szUnit == "CM" || szUnit == "cm" )
	{
		FACTOR_T = CM_TO_METER;
	}
	else if ( szUnit == "MM" || szUnit == "mm" )
	{
		FACTOR_T = MM_TO_METER;
	}
	else if ( szUnit == "INCH" || szUnit == "inch" )
	{
		FACTOR_T = INCH_TO_METER;
	}
	else if ( szUnit == "FOOT" || szUnit == "foot" )
	{
		FACTOR_T = FOOT_TO_METER;
	}
	else
	{
		std::cout<< "Unknown Unit " << szUnit << " Is NOT Supported!" << std::endl;
	}
	return;
}

void CFastTIDriver::writeOutput()
{
	// Write H_log.DAT file
	//m_fOutFile.open( "H_log.DAT" );
	//if ( !m_fOutFile )
	//{
	//	std::cout<< "Unable to Open output File: " << "H_log.DAT" << std::endl;
	//	return;
	//}

	//// Write Aresis.DAT file
	//std::ofstream fResistOut;
	//fResistOut.open( "Aresis.DAT" );
	//if ( !fResistOut )
	//{
	//	std::cout<< "Unable to Open output File: " << "Aresis.DAT" << std::endl;

	//	return;
	//}
	//
	//m_fOutFile << setw(15)<< "DEPTH" 
	//		   << setw(15) << "REAL(Hxx)" << setw(15) << "IMAG(Hxx)" << setw(15) << "REAL(Hxy)" << setw(15) << "IMAG(Hxy)" << setw(15) << "REAL(Hxz)" << setw(15) << "IMAG(Hxz)"
	//		   << setw(15) << "REAL(Hyx)" << setw(15) << "IMAG(Hyx)" << setw(15) << "REAL(Hyy)" << setw(15) << "IMAG(Hyy)" << setw(15) << "REAL(Hyz)" << setw(15) << "IMAG(Hyz)"
	//		   << setw(15) << "REAL(Hzx)" << setw(15) << "IMAG(Hzx)" << setw(15) << "REAL(Hzy)" << setw(15) << "IMAG(Hzy)" << setw(15) << "REAL(Hzz)" << setw(15) << "IMAG(Hzz)" << std::endl;

	//fResistOut << setw(15)<< "DEPTH" 
	//		   << setw(15) << "REAL(Rxx)" << setw(15) << "IMAG(Rxx)" << setw(15) << "REAL(Ryy)" << setw(15) << "IMAG(Ryy)" << setw(15) << "REAL(Rzz)" << setw(15) << "IMAG(Rzz)" << std::endl;

	//char buffer [33];
	//unsigned int iCnt = 0;
	//DataType AK[9] = {0.0};
	//DataType BK[9] = {0.0};
	//DataType AK_T[9] = {0.0};
	//DataType BK_T[9] = {0.0};
	//for( unsigned int iTran = 0; iTran < m_sToolConf.s_iNum_Trans; iTran++ )
	//{
	//	for( unsigned int iRecv = 0; iRecv < m_sToolConf.s_iNum_Recvs; iRecv++ )
	//	{
	//		DataType AL = abs( m_sToolConf.s_vDistRecvs[iRecv] - m_sToolConf.s_vDistTrans[iTran] );
	//		AK[8] = -( m_dOmega * MU_0 ) * ( m_dOmega * MU_0 )/( 4.0 * PI * AL );
	//		AK[0] = 0.5 * AK[8];			AK[1] = 0.5 * AK[8];			AK[3] = 0.5 * AK[8];			AK[4] = 0.5 * AK[8];
	//		AK[2] = 0.25 * AK[8];			AK[5] = 0.25 * AK[8];			AK[6] = 0.25 * AK[8];			AK[7] = 0.25 * AK[8];

	//		BK[8] = m_dOmega * MU_0 / ( 2.0 * PI * AL * AL * AL );
	//		BK[0] = - 0.5 * BK[8];			BK[4] = - 0.5 * BK[8];

	//		for( unsigned int it = 0; it < 9; it++ )
	//		{
	//			AK_T[it] = AK_T[it] + AK[it] * m_sToolConf.s_iNum_TransTurns[iTran] * m_sToolConf.s_iNum_RecvsTurns[iRecv];
	//			BK_T[it] = BK_T[it] + BK[it] * m_sToolConf.s_iNum_TransTurns[iTran] * m_sToolConf.s_iNum_RecvsTurns[iRecv];
	//		}
	//	}
	//}

	//for( unsigned int iLgPt = 0; iLgPt < m_NumLgPts; iLgPt++ )
	//{
	//	CDyadic<Complex> Hfield;
	//	DataType zPos = -( m_dDepthSt + iLgPt * m_dStepSz ) / m_sFormConf.s_FACTOR_T;
	//	for( unsigned int iTran = 0; iTran < m_sToolConf.s_iNum_Trans; iTran++ )
	//	{
	//		for( unsigned int iRecv = 0; iRecv < m_sToolConf.s_iNum_Recvs; iRecv++ )
	//		{
	//			iCnt = ( iTran * m_sToolConf.s_iNum_Recvs  + iRecv ) * m_NumLgPts + iLgPt;
	//			Hfield = m_vHfieldDyd[iCnt] + Hfield;
	//		}
	//	}

	//	m_fOutFile << std::setiosflags( ios::scientific | ios::right ) << setw(15)<< zPos << setw(15)
	//				<< setw(15) << Hfield.xx.real() << setw(15) << Hfield.xx.imag()
	//				<< setw(15) << Hfield.xy.real() << setw(15) << Hfield.xy.imag()
	//				<< setw(15) << Hfield.xz.real() << setw(15) << Hfield.xz.imag()
	//				<< setw(15) << Hfield.yx.real() << setw(15) << Hfield.yx.imag()
	//				<< setw(15) << Hfield.yy.real() << setw(15) << Hfield.yy.imag()
	//				<< setw(15) << Hfield.yz.real() << setw(15) << Hfield.yz.imag()
	//				<< setw(15) << Hfield.zx.real() << setw(15) << Hfield.zx.imag()
	//				<< setw(15) << Hfield.zy.real() << setw(15) << Hfield.zy.imag()
	//				<< setw(15) << Hfield.zz.real() << setw(15) << Hfield.zz.imag() << std::endl;

	//	DataType AI_xx = m_dOmega * MU_0 * Hfield.xx.real() - BK_T[0];
	//	DataType AI_yy = m_dOmega * MU_0 * Hfield.yy.real() - BK_T[4];
	//	DataType AI_zz = m_dOmega * MU_0 * Hfield.zz.real() - BK_T[8];
	//	DataType AR_xx = - m_dOmega * MU_0 * Hfield.xx.imag() / AK_T[0];
	//	DataType AR_yy = - m_dOmega * MU_0 * Hfield.yy.imag() / AK_T[4];
	//	DataType AR_zz = - m_dOmega * MU_0 * Hfield.zz.imag() / AK_T[8];

	//	fResistOut << std::setiosflags( ios::scientific | ios::right ) << setw(15)<< zPos << setw(15)
	//				<< 1.0/AR_xx << setw(15) << AK_T[0]/AI_xx << setw(15) 
	//				<< 1.0/AR_yy << setw(15) << AK_T[4]/AI_yy << setw(15) 
	//				<< 1.0/AR_zz << setw(15) << AK_T[8]/AI_zz << setw(15) << std::endl;
	//}
	//m_fOutFile.close();
	//fResistOut.close();

	// Write HLog_ij.DAT file
	for( unsigned int iTran = 0; iTran < m_sToolConf.s_iNum_Trans; iTran++ )
	{
		for( unsigned int iRecv = 0; iRecv < m_sToolConf.s_iNum_Recvs; iRecv++ )
		{
			string szFile = "Fast1D2017_Hlog_" + std::to_string(iTran + 1) + std::to_string(iRecv + 1) + ".DAT";
			m_fOutFile.open( szFile );
			if ( !m_fOutFile )
			{
				std::cout<< "Unable to Open output File: " << szFile << std::endl;
				return;
			}
			m_fOutFile << setw(15)<< "ZR" 
				<< setw(15) << "REAL(Hxx)" << setw(15) << "IMAG(Hxx)" << setw(15) << "REAL(Hxy)" << setw(15) << "IMAG(Hxy)" << setw(15) << "REAL(Hxz)" << setw(15) << "IMAG(Hxz)"
				<< setw(15) << "REAL(Hyx)" << setw(15) << "IMAG(Hyx)" << setw(15) << "REAL(Hyy)" << setw(15) << "IMAG(Hyy)" << setw(15) << "REAL(Hyz)" << setw(15) << "IMAG(Hyz)"
				<< setw(15) << "REAL(Hzx)" << setw(15) << "IMAG(Hzx)" << setw(15) << "REAL(Hzy)" << setw(15) << "IMAG(Hzy)" << setw(15) << "REAL(Hzz)" << setw(15) << "IMAG(Hzz)" << std::endl;

			for( unsigned int iLgPt = 0; iLgPt < m_NumLgPts; iLgPt++ )
			{
				DataType zPos = -( m_dDepthSt + iLgPt * m_dStepSz ) / m_sFormConf.s_FACTOR_T;

				int iCnt = (iTran * m_sToolConf.s_iNum_Recvs + iRecv) * m_NumLgPts + iLgPt;
				m_fOutFile << std::setiosflags( ios::scientific | ios::right ) << setw(15) << zPos << setw(15)
					<< setw(15) << m_vHfieldDyd[iCnt].xx.real() << setw(15) << m_vHfieldDyd[iCnt].xx.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].xy.real() << setw(15) << m_vHfieldDyd[iCnt].xy.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].xz.real() << setw(15) << m_vHfieldDyd[iCnt].xz.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].yx.real() << setw(15) << m_vHfieldDyd[iCnt].yx.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].yy.real() << setw(15) << m_vHfieldDyd[iCnt].yy.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].yz.real() << setw(15) << m_vHfieldDyd[iCnt].yz.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].zx.real() << setw(15) << m_vHfieldDyd[iCnt].zx.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].zy.real() << setw(15) << m_vHfieldDyd[iCnt].zy.imag()
					<< setw(15) << m_vHfieldDyd[iCnt].zz.real() << setw(15) << m_vHfieldDyd[iCnt].zz.imag() << std::endl;
			}
			m_fOutFile.close();
		}
	}
	
}


// modified initialization functions
bool CFastTIDriver::initToolConfig( const DataType freq, const DataType locTrans, const DataType locRecvs )
{
	m_dFreq = freq;
	m_dOmega = m_dFreq * TWOPI;

	// Number of Transmitter and Receivers
	m_sToolConf.s_iNum_Trans = 1;
	m_sToolConf.s_iNum_Recvs = 1;

	m_sToolConf.s_iNum_TransTurns.resize(m_sToolConf.s_iNum_Trans);
	m_sToolConf.s_vDistTrans.resize(m_sToolConf.s_iNum_Trans);
	m_sToolConf.s_iNum_RecvsTurns.resize(m_sToolConf.s_iNum_Recvs);
	m_sToolConf.s_vDistRecvs.resize(m_sToolConf.s_iNum_Recvs);

	// Input the unit for tool
	std::string szUnit;
	szUnit = "inch";
	convertUnit( szUnit, m_sToolConf.s_FACTOR_T );

	// Input the transmitters and receivers info
	for( unsigned int i = 0; i < m_sToolConf.s_iNum_Trans; i++ )
	{
		m_sToolConf.s_iNum_TransTurns[i] = 1;
		m_sToolConf.s_vDistTrans[i] = locTrans;
		m_sToolConf.s_vDistTrans[i] = - m_sToolConf.s_vDistTrans[i] * m_sToolConf.s_FACTOR_T;
	}

	for( unsigned int i = 0; i < m_sToolConf.s_iNum_Recvs; i++ )
	{
		m_sToolConf.s_iNum_RecvsTurns[i] = 1;
		m_sToolConf.s_vDistRecvs[i] = locRecvs;
		m_sToolConf.s_vDistRecvs[i] = - m_sToolConf.s_vDistRecvs[i] * m_sToolConf.s_FACTOR_T;
	}

	DataType TRdist = 0;
	for( unsigned int itran = 0; itran < m_sToolConf.s_iNum_Trans; itran++ )
	{
		for( unsigned int jrecv = 0; jrecv < m_sToolConf.s_iNum_Recvs; jrecv++ )
		{
			TRdist = abs( m_sToolConf.s_vDistRecvs[jrecv] - m_sToolConf.s_vDistTrans[itran] );
			m_sToolConf.s_dMaxDist = max( TRdist, m_sToolConf.s_dMaxDist );
		}
	}
	return true;
}

bool CFastTIDriver::initFormationPara( const DataType Dip, const DataType Azimuth, const unsigned int numLayers, const DataType* ZBed, const DataType* Rh, const DataType* Rv, const DataType* traj )
{
	DataType Alfa, Beta, Gama;
	Alfa = Dip;
	Beta = Azimuth; 
	Gama = 0;

	// determine the well type

	Alfa = fmod(Alfa, 360.0);

	if( Alfa == 0.0 || Alfa == 180.0)
	{
		m_WellType = ENUM_VERTICAL;
	}
	else if (Alfa == 90.0 || Alfa == 270.0)
	{
		m_WellType = ENUM_HORIZONTAL;
	}
	else
	{
		m_WellType = ENUM_DEVIATED;
	}

	// convert from degree to radius
	Alfa = Alfa * DEG_TO_RAD;
	Beta = Beta * DEG_TO_RAD;
	Gama = Gama * DEG_TO_RAD;
	// calculate rotation matrix R
	m_sFormConf.s_dyRotMatrixR.xx =  cos(Alfa)*cos(Beta)*cos(Gama) - sin(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixR.xy = -cos(Alfa)*cos(Beta)*sin(Gama) - sin(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixR.xz =  sin(Alfa)*cos(Beta);
	m_sFormConf.s_dyRotMatrixR.yx =  cos(Alfa)*sin(Beta)*cos(Gama) + cos(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixR.yy = -cos(Alfa)*sin(Beta)*sin(Gama) + cos(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixR.yz =  sin(Alfa)*sin(Beta);
	m_sFormConf.s_dyRotMatrixR.zx = -sin(Alfa)*cos(Gama);
	m_sFormConf.s_dyRotMatrixR.zy =  sin(Alfa)*sin(Gama);
	m_sFormConf.s_dyRotMatrixR.zz =  cos(Alfa);

	// calculate inverse rotation matrix R
	m_sFormConf.s_dyRotMatrixInvR.xx = cos(Alfa)*cos(Beta)*cos(Gama)-sin(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixInvR.xy = cos(Alfa)*sin(Beta)*cos(Gama)+cos(Beta)*sin(Gama);
	m_sFormConf.s_dyRotMatrixInvR.xz =-sin(Alfa)*cos(Gama);
	m_sFormConf.s_dyRotMatrixInvR.yx =-cos(Alfa)*cos(Beta)*sin(Gama)-sin(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixInvR.yy =-cos(Alfa)*sin(Beta)*sin(Gama)+cos(Beta)*cos(Gama);
	m_sFormConf.s_dyRotMatrixInvR.yz = sin(Alfa)*sin(Gama);
	m_sFormConf.s_dyRotMatrixInvR.zx = sin(Alfa)*cos(Beta);
	m_sFormConf.s_dyRotMatrixInvR.zy = sin(Alfa)*sin(Beta);
	m_sFormConf.s_dyRotMatrixInvR.zz = cos(Alfa);

	m_sFormConf.s_dAlfa = Alfa;
	m_sFormConf.s_dBeta = Beta;
	m_sFormConf.s_dGama = Gama;

	// Input the unit for formation
	std::string szUnit;
	szUnit = "foot";
	convertUnit( szUnit, m_sFormConf.s_FACTOR_T );

	int sign;
	sign = 0;
	m_sFormConf.s_bSign = ( sign == 0 ? false : true );

	m_sFormConf.s_iNum_Layers = numLayers;
	m_sFormConf.s_vHResistivity.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVResistivity.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHRelEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVRelEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHRelMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVRelMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHCompEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVCompEpsilon.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHCompMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vVCompMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vjOmegaVMu.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vTMAnsiRatio.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vTEAnsiRatio.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vHWaveNum.resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vStaticFresnelDn[0].resize(m_sFormConf.s_iNum_Layers+2);
	m_sFormConf.s_vStaticFresnelDn[1].resize(m_sFormConf.s_iNum_Layers + 2);

	m_sFormConf.s_K0 = m_dOmega / VLIGHT_0; 

	// Need Check ?????
	m_sFormConf.s_vLayerBoundary.resize(m_sFormConf.s_iNum_Layers+1);
	m_sFormConf.s_vLayerThick.resize(m_sFormConf.s_iNum_Layers+1);

	for( unsigned int i = 1; i <= m_sFormConf.s_iNum_Layers; i++ )
	{
		if( i != m_sFormConf.s_iNum_Layers )
		{
			m_sFormConf.s_vLayerBoundary[i] = ZBed[i-1];
			m_sFormConf.s_vHResistivity[i] = Rh[i-1];
			m_sFormConf.s_vVResistivity[i] = Rv[i-1];
			m_sFormConf.s_vHRelEpsilon[i] = 1;
			m_sFormConf.s_vVRelEpsilon[i] = 1;
			m_sFormConf.s_vHRelMu[i] = 1;

			m_sFormConf.s_vVRelMu[i] = m_sFormConf.s_vHRelMu[i];
			// Need Check ?????
			m_sFormConf.s_vLayerBoundary[i] = - m_sFormConf.s_vLayerBoundary[i];
		}
		else
		{
			m_sFormConf.s_vHResistivity[i] = Rh[m_sFormConf.s_iNum_Layers-1];
			m_sFormConf.s_vVResistivity[i] = Rv[m_sFormConf.s_iNum_Layers-1];
			m_sFormConf.s_vHRelEpsilon[i] = 1;
			m_sFormConf.s_vVRelEpsilon[i] = 1;
			m_sFormConf.s_vHRelMu[i] = 1;
		}

		m_sFormConf.s_vHCompEpsilon[i] = m_sFormConf.s_vHRelEpsilon[i] * EPSILON_0 - CPLX_J * (CPLX_ONE/m_sFormConf.s_vHResistivity[i] ) / m_dOmega;
		m_sFormConf.s_vVCompEpsilon[i] = m_sFormConf.s_vVRelEpsilon[i] * EPSILON_0 - CPLX_J * (CPLX_ONE/m_sFormConf.s_vVResistivity[i] ) / m_dOmega;
	
		m_sFormConf.s_vHCompMu[i] = m_sFormConf.s_vHRelMu[i] * MU_0;
		m_sFormConf.s_vVRelMu[i] = m_sFormConf.s_vHRelMu[i];
		m_sFormConf.s_vVCompMu[i] = m_sFormConf.s_vVRelMu[i] * MU_0;
		m_sFormConf.s_vjOmegaVMu[i] = CPLX_J * m_dOmega * m_sFormConf.s_vVCompMu[i];

		m_sFormConf.s_vTMAnsiRatio[i] = m_sFormConf.s_vHCompEpsilon[i]/m_sFormConf.s_vVCompEpsilon[i];
		m_sFormConf.s_vTEAnsiRatio[i] = m_sFormConf.s_vHCompMu[i]/m_sFormConf.s_vVCompMu[i];
		m_sFormConf.s_vHWaveNum[i] = m_dOmega * sqrt( m_sFormConf.s_vHCompEpsilon[i] * m_sFormConf.s_vHCompMu[i] );

		m_sFormConf.s_vLayerBoundary[i] = m_sFormConf.s_vLayerBoundary[i] * m_sFormConf.s_FACTOR_T * cos(m_sFormConf.s_bSign*m_sFormConf.s_dAlfa);
		
	}

	m_sFormConf.s_vHResistivity[0] = m_sFormConf.s_vHResistivity[1];
	m_sFormConf.s_vHResistivity[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHResistivity[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVResistivity[0] = m_sFormConf.s_vVResistivity[1];
	m_sFormConf.s_vVResistivity[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVResistivity[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHRelEpsilon[0] = m_sFormConf.s_vHRelEpsilon[1];
	m_sFormConf.s_vHRelEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHRelEpsilon[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVRelEpsilon[0] = m_sFormConf.s_vVRelEpsilon[1];
	m_sFormConf.s_vVRelEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVRelEpsilon[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHRelMu[0] = m_sFormConf.s_vHRelMu[1];
	m_sFormConf.s_vHRelMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHRelMu[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVRelMu[0] = m_sFormConf.s_vVRelMu[1];
	m_sFormConf.s_vVRelMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVRelMu[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHCompEpsilon[0] = m_sFormConf.s_vHCompEpsilon[1];
	m_sFormConf.s_vHCompEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHCompEpsilon[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVCompEpsilon[0] = m_sFormConf.s_vVCompEpsilon[1];
	m_sFormConf.s_vVCompEpsilon[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVCompEpsilon[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHCompMu[0] = m_sFormConf.s_vHCompMu[1];
	m_sFormConf.s_vHCompMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHCompMu[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vVCompMu[0] = m_sFormConf.s_vVCompMu[1];
	m_sFormConf.s_vVCompMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vVCompMu[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vjOmegaVMu[0] = m_sFormConf.s_vjOmegaVMu[1];
	m_sFormConf.s_vjOmegaVMu[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vjOmegaVMu[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vTMAnsiRatio[0] = m_sFormConf.s_vTMAnsiRatio[1];
	m_sFormConf.s_vTMAnsiRatio[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vTMAnsiRatio[m_sFormConf.s_iNum_Layers];
	m_sFormConf.s_vTEAnsiRatio[0] = m_sFormConf.s_vTEAnsiRatio[1];
	m_sFormConf.s_vTEAnsiRatio[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vTEAnsiRatio[m_sFormConf.s_iNum_Layers];

	m_sFormConf.s_vHWaveNum[0] = m_sFormConf.s_vHWaveNum[1];
	m_sFormConf.s_vHWaveNum[m_sFormConf.s_iNum_Layers+1] = m_sFormConf.s_vHWaveNum[m_sFormConf.s_iNum_Layers];

	for( unsigned int i = 1; i <= m_sFormConf.s_iNum_Layers; i++ )
	{
		// for TM mode;
		//DataType TM_Tmp_L = sqrt( m_sFormConf.s_vHRelEpsilon[i] * m_sFormConf.s_vVRelEpsilon[i] );
		//DataType TM_Tmp_LP1 = sqrt( m_sFormConf.s_vHRelEpsilon[i+1] * m_sFormConf.s_vVRelEpsilon[i+1] );
		
		Complex TM_Tmp_L = sqrt(m_sFormConf.s_vHCompEpsilon[i] * m_sFormConf.s_vVCompEpsilon[i]);
		Complex TM_Tmp_LP1 = sqrt(m_sFormConf.s_vHCompEpsilon[i + 1] * m_sFormConf.s_vVCompEpsilon[i + 1]);
		m_sFormConf.s_vStaticFresnelDn[0][i] = (TM_Tmp_L - TM_Tmp_LP1) / (TM_Tmp_L + TM_Tmp_LP1);

		// for TE mode;
		//DataType TE_Tmp_L = sqrt( m_sFormConf.s_vHRelMu[i] * m_sFormConf.s_vVRelMu[i] );
		//DataType TE_Tmp_LP1 = sqrt( m_sFormConf.s_vHRelMu[i+1] * m_sFormConf.s_vVRelMu[i+1] );
		Complex TE_Tmp_L = sqrt(m_sFormConf.s_vHCompMu[i] * m_sFormConf.s_vVCompMu[i]);
		Complex TE_Tmp_LP1 = sqrt(m_sFormConf.s_vHCompMu[i + 1] * m_sFormConf.s_vVCompMu[i + 1]);
		m_sFormConf.s_vStaticFresnelDn[1][i] = (TE_Tmp_LP1 - TE_Tmp_L) / (TE_Tmp_L + TE_Tmp_LP1);
	}

	//m_sFormConf.s_vLayerBoundary[0] = m_sFormConf.s_vLayerBoundary[1];
	//m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers] = m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers-1];

	m_dDepthSt = traj[0];
	m_dStepSz = traj[1];
	m_dDepthEd = traj[2];
	
	m_NumLgPts = (m_dDepthEd - m_dDepthSt)/m_dStepSz + 1;

	m_dDepthSt = - m_dDepthSt * m_sFormConf.s_FACTOR_T;
	m_dStepSz = - m_dStepSz * m_sFormConf.s_FACTOR_T;
	m_dDepthEd = - m_dDepthEd * m_sFormConf.s_FACTOR_T;

	DataType maxLogz = max(m_dDepthSt, m_dDepthEd);
	DataType minLogz = min(m_dDepthSt, m_dDepthEd);
	DataType buff = 10.0;

	if( maxLogz < m_sFormConf.s_vLayerBoundary[1] )
		m_sFormConf.s_vLayerBoundary[0] = m_sFormConf.s_vLayerBoundary[1] + buff;
	else
		m_sFormConf.s_vLayerBoundary[0] = maxLogz + m_sToolConf.s_dMaxDist + buff;

	// [D.W.] Modify to calculate s_vLayerBoundary[m_sFormConf.s_iNum_Layers] correctly

	if( minLogz > m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers-1] )
		m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers] = m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers-1] - buff;
	else
		m_sFormConf.s_vLayerBoundary[m_sFormConf.s_iNum_Layers] = minLogz - m_sToolConf.s_dMaxDist - buff;
	
	for( unsigned int i = 1; i <= m_sFormConf.s_iNum_Layers; i++ )
		m_sFormConf.s_vLayerThick[i] = m_sFormConf.s_vLayerBoundary[i-1] - m_sFormConf.s_vLayerBoundary[i];

	return true;
}