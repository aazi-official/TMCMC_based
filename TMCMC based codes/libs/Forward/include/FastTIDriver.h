#ifndef FASTTIDRIVER_H
#define FASTTIDRIVER_H

#include <string>
#include <iostream>
#include <fstream>
#include <bitset>
#include <iomanip>
#include <algorithm> 
#include <math.h> 

//#include <omp.h>

#include "CommUtil.h"
#include "TxLine.h"
#include "QuadRule.h"

class CFastTIDriver
{
public:

	CFastTIDriver();

	~CFastTIDriver();

	bool initToolConfig( in_ const std::string & szToolFile );

	bool initFormationPara( in_ const std::string & szFormationFile );

	bool calcLogs();

	void writeOutput();

	// modified parts
	bool initToolConfig( const DataType freq, const DataType locTrans, const DataType locRecvs );
	bool initFormationPara( const DataType Dip, const DataType Azimuth, const unsigned int numLayers, const DataType* ZBed, const DataType* Rh, const DataType* Rv, const DataType* traj );
	vector<CDyadic<Complex> > m_vHfieldDyd;

private:

	void convertUnit( in_ const string & szUnit, out_ DataType & FACTOR_T );

	void calcHANKTrans( in_ const unsigned int iHTrans, in_ const SObsSrcPair OSPair, out_ Complex & rslt );

	void getIntgrand( in_ const unsigned int iHTrans, in_ const unsigned int iWght, in_ const DataType krho, in_ const SObsSrcPair OSPair, out_ Complex & rslt );

	void locateLayerNum( in_ const SFormationConfig & FormConf, in_ const DataType & z, out_ unsigned int & iLayer, out_ bool & onTopFlag );

	void calcAllocKrhoVec(in_ const unsigned int iStart, in_ const unsigned int iEnd, const DataType & rho);

	void precalcTxLineVIs( in_ const unsigned int iStart, in_ const unsigned int iEnd, in_ const SObsSrcPair OSPair );

private:

	std::ifstream m_fTOOLIn;
	std::ifstream m_fFormIn;
	std::ofstream m_fOutFile;

	DataType m_dFreq;
	DataType m_dOmega;

	SToolConfig m_sToolConf;
	SFormationConfig m_sFormConf;

	EWellType m_WellType;
	DataType m_dDepthSt;
	DataType m_dDepthEd;
	DataType m_dStepSz;
	unsigned int m_NumLgPts;

	CTxLine m_TxLineObj;
	vector<DataType> m_vSpckrho;

	std::bitset<HTRANS_PTS> m_SpcFlg;
	vector<STxLineDataType> m_vTxLineParas;
	omp_lock_t m_TxParalock;

	std::bitset<HTRANS_PTS> m_VIsFlg;
	vector<STxLineVI> m_vTxLineVIs;
	omp_lock_t m_TxVIlock;


	unsigned int m_iHanksMin;
	unsigned int m_iHanksMax;
	vector<unsigned int> m_vHanksMin;
	vector<unsigned int> m_vHanksMax;

};

#endif