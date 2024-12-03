#ifndef TXLINE_H
#define TXLINE_H

#include "CommUtil.h"

class CTxLine
{
public:

	CTxLine();

	~CTxLine();

	bool calcTxLineVIs(		in_ const SFormationConfig & ForConfig,	in_ const STxLineDataType & TxLinePara,
							in_ const SObsSrcPair OSPair,			in_ const DataType & krho,
							out_ STxLineVI & VI );

	bool calcAsympTxLineVI( in_ const SFormationConfig & ForConfig,	in_ const STxLineDataType & TxLinePara,
							in_ const SObsSrcPair OSPair,			in_ const DataType & krho,
							out_ STxLineVI & AsympVI );

	bool calcExtractTerm(	in_ const SFormationConfig & ForConfig,	in_ const SObsSrcPair & OSPair, in_ const DataType & Omega,
							out_ std::vector<Complex> & rsltVec );

public:

	bool calcTxLineParas(	in_ const SFormationConfig & ForConfig,	in_ const DataType & Omega,	in_ const DataType & krho,  out_ STxLineDataType & TxLinePara );

	bool calTLVIatSrcLayer(	in_ const SFormationConfig & ForConfig,	in_ const STxLineDataType & TxLinePara,
								in_ const unsigned int & iSrcLayerNum,	in_ const DataType & Obs_z,	in_ const DataType & Src_zPrim, 
								out_ STxLineVI & VI );

	bool analyticSRITerm(	in_ const Complex & wavenumber,	in_ const Complex & nu,	
							in_ const DataType & rho, out_ std::vector<Complex> & SRI );

};

#endif