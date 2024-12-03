#ifndef CURVES_H
#define CURVES_H

#include "CommUtil.h"
#include "FastTIDriver.h"

void dipole1D(out_ complex<DataType> *response, in_ CFastTIDriver &fastTIDriverObj, in_ unsigned int numLayers, 
				in_ const DataType *Rh, in_ const DataType *Rv, in_ const DataType *ZBed, in_ const DataType Dip, in_ const DataType Azimuth, 
				in_ const DataType TVD, in_ const DataType freq, in_ DataType Trans, in_ DataType Recvs);

void directional_curve(DataType* curves, CFastTIDriver &fastTIDriverObj, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed,
	const DataType TVD, const DataType Dip, const DataType Azimuth, const DataType freq, const DataType sp);

void forward_geosphere(out_ DataType* resp, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed, const DataType TVD, const DataType Dip);

void forward_chevron(out_ DataType* resp, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed, const DataType TVD, const DataType Dip);

void forward(out_ DataType* resp, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed, const DataType TVD, const DataType Dip, const int freq);

#endif