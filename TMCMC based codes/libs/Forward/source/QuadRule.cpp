#include "QuadRule.h"

int CQuadRule::iCount = 0;

CQuadRule::CQuadRule(CTxLine & TxLine, const SFormationConfig & ForConfig, DataType Omega, DataType EPS, DataType KAI, int MAXLEVEL,int NMAX, DataType Eta)
	: m_TxLine(TxLine), m_ForConfig(ForConfig), m_Omega(Omega),m_EPS(EPS), m_KAI(KAI), m_MAXLEVEL(MAXLEVEL), m_NMAX(NMAX),m_Eta(Eta)
{
	m_TxLineDT.Init(m_ForConfig.s_iNum_Layers);
}

CQuadRule::~CQuadRule()
{
}


bool CQuadRule::tanhsinhQuad( const DataType & a, const DataType & b, const SObsSrcPair& OSPair, Complex val[] )
{
	m_DE_a = a;
	m_DE_b = b;

	DataType sigma = ( b - a )/2.0;
	m_DE_sgm = sigma;

	DataType Gama = ( b + a )/2.0;

	int m = 0;
	DataType h = 1.5;

	DataType eh = exp(h);

	// call the function;
	Complex Integrand[2];
	calcIntgrnd( Gama, OSPair, Integrand );
	Complex s[2];

	s[0] = m_Eta * Integrand[0];
	s[1] = m_Eta * Integrand[1];
	int n = 0;
	TruncIndex( OSPair, eh, n, s, true );
	Complex old[2];
	old[0] = sigma * h * s[0];
	old[1] = sigma * h * s[1];

	while( m < m_MAXLEVEL )
	{
		m++;
		DataType e2h = eh;
		h = h/2.0;
		eh = exp(h);
		PartSum( OSPair, eh, e2h, n, s, true );
		val[0] = old[0]/2.0 + sigma * h * s[0];
		val[1] = old[1]/2.0 + sigma * h * s[1];

		if( abs(val[0] - old[0]) < m_EPS*abs(val[0]) && abs(val[1] - old[1]) < m_EPS*abs(val[1]) )
		{
			return true;
		}
		old[0] = val[0];
		old[1] = val[1];
		n = 2 * n;
	}

	return true;
}

bool CQuadRule::mixedQuad( const DataType & a, const SObsSrcPair& OSPair, Complex val[] )
{
	m_MX_a = a;

	int m = 0;
	DataType h = 1.0;

	DataType eh = exp(h);
	DataType delta = exp(-1.0);
	DataType wght = 2.0 * delta;

	// call the function;
	Complex s[2];
	calcIntgrnd( a + delta, OSPair, s );
	s[0] = wght * s[0];
	s[1] = wght * s[1];

	int n1 = 0;
	int n2 = 0;
	TruncIndex( OSPair, eh, n1, s, false );
	TruncIndex( OSPair, 1.0/eh, n2, s, false );
	Complex old[2];
	old[0] = h * s[0];
	old[1] = h * s[1];

	while( m < m_MAXLEVEL )
	{
		m++;
		DataType e2h = eh;
		h = h/2.0;
		eh = exp(h);
		Complex s1[2], s2[2];
		PartSum( OSPair, eh, e2h, n1, s1, false );
		PartSum( OSPair, 1.0/eh, 1.0/e2h, n2, s2, false );

		val[0] = old[0]/2.0 + h * (s1[0] + s2[0]);
		val[1] = old[1]/2.0 + h * (s1[1] + s2[1]);

		if( abs(val[0] - old[0]) < m_EPS*abs(val[0]) && abs(val[1] - old[1]) < m_EPS*abs(val[1]) )
		{
			return true;
		}
		old[0] = val[0];
		old[1] = val[1];
		n1 = 2 * n1;
		n2 = 2 * n2;
	}
	return true;
}


void CQuadRule::calcIntgrnd(const DataType & Gama, const SObsSrcPair& OSPair, Complex I[])
{
	iCount = iCount + 1;
	
	m_TxLine.calcTxLineParas(m_ForConfig, m_Omega, Gama, m_TxLineDT);
	m_TxLine.calcTxLineVIs(m_ForConfig, m_TxLineDT, OSPair, Gama, m_TxLineVI);

	I[0] = (m_TxLineVI.s_Iv[0] + m_TxLineVI.s_Iv[1])*Gama;
	I[1] = (m_TxLineVI.s_Vi[1])*Gama*Gama*Gama;
}

void CQuadRule::TruncIndex( const SObsSrcPair& OSPair, const DataType & eh, int & n, Complex s[], const bool & DEFlag )
{
	DataType ekh = eh;
	Complex t[2];
	n = 0;
	while( n < m_NMAX )
	{
		n++;
		Term(OSPair, ekh, t, DEFlag);
		s[0] = s[0] + t[0];
		s[1] = s[1] + t[1];
		if( abs(t[0]) < m_KAI*abs(s[0]) && abs(t[1]) < m_KAI*abs(s[1]) )
		{
			return;
		}
		ekh = ekh * eh;
	}
}

void CQuadRule::PartSum( const SObsSrcPair& OSPair, const DataType & eh, const DataType & e2h, const int n, Complex s[], const bool & DEFlag )
{
	DataType ekh = eh;
	Complex t[2];
	Term(OSPair, ekh, s, DEFlag);
	for( int k = 2; k<= n; k++ )
	{
		ekh = ekh * e2h;
		Term(OSPair, ekh, t, DEFlag);
		s[0] = s[0] + t[0];
		s[1] = s[1] + t[1];
	}
}

void CQuadRule::Term( const SObsSrcPair& OSPair, const DataType & ekh, Complex t[], const bool & DEFlag )
{
		
	if( DEFlag )
	{
		DataType q = exp(- m_Eta * ( ekh - 1./ekh ) );
		DataType delta = 2.0 * q /( 1.0 + q);
		DataType wght = m_Eta * ( ekh + 1.0/ekh ) * delta/( 1.0 + q);
		//call the function
		DataType Gama1 = m_DE_a + delta * m_DE_sgm;
		Complex t1[2];
		calcIntgrnd( Gama1, OSPair, t1);

		DataType Gama2 = m_DE_b - delta * m_DE_sgm;
		Complex t2[2];
		calcIntgrnd( Gama2, OSPair, t2);

		t[0] = wght * (t1[0] + t2[0]);
		t[1] = wght * (t1[1] + t2[1]);
	}
	else
	{
		DataType delta = exp(- ekh )/ekh;
		DataType wght = ( 1 + ekh)*delta;
		// call the function
		DataType Gama = m_MX_a + delta;
		calcIntgrnd( Gama, OSPair, t);
		t[0] = wght * t[0];
		t[1] = wght * t[1];
	}
}