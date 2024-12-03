#ifndef QUADRULE_H
#define QUADRULE_H

#include "CommUtil.h"
#include "TxLine.h"

class CQuadRule
{

public:

	CQuadRule(CTxLine & TxLine, const SFormationConfig & ForConfig, DataType Omega, DataType EPS = 1E-7, DataType KAI = 1E-15, int MAXLEVEL = 5,int NMAX = 24, DataType Eta = 1);

	~CQuadRule();

	bool tanhsinhQuad( const DataType & a, const DataType & b, const SObsSrcPair& OSPair, Complex val[] );

	bool mixedQuad( const DataType & a, const SObsSrcPair& OSPair, Complex val[] );

	static int iCount;

private:

	void calcIntgrnd(const DataType & Gama, const SObsSrcPair& OSPair, Complex I[]);

	void TruncIndex( const SObsSrcPair& OSPair, const DataType & eh, int & n, Complex s[], const bool & DEFlag );

	void PartSum( const SObsSrcPair& OSPair, const DataType & eh, const DataType & e2h, const int n, Complex s[], const bool & DEFlag );

	void Term( const SObsSrcPair& OSPair, const DataType & ekh, Complex t[], const bool & DEFlag );

private:

	DataType m_DE_a;
	DataType m_DE_b;
	DataType m_DE_sgm;

	DataType m_MX_a;

	int m_MAXLEVEL;

	int m_NMAX;
	DataType m_KAI;
	DataType m_EPS;

	DataType m_Eta;

	CTxLine & m_TxLine;

	const SFormationConfig & m_ForConfig;

	DataType m_Omega;

	STxLineDataType m_TxLineDT;
	STxLineVI m_TxLineVI;

};

#endif