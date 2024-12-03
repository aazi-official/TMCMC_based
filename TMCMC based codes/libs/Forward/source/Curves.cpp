#include "Curves.h"

void dipole1D(out_ complex<DataType> *response, in_ CFastTIDriver &fastTIDriverObj, in_ unsigned int numLayers, in_ const DataType *Rh, in_ const DataType *Rv, 
	in_ const DataType *ZBed, in_ const DataType TVD, in_ const DataType Dip, in_ const DataType Azimuth, in_ const DataType freq, in_ DataType Trans, in_ DataType Recvs) {
	DataType traj[3] = {TVD,1,TVD};

	fastTIDriverObj.initToolConfig(freq, Trans, Recvs);
	fastTIDriverObj.initFormationPara(Dip, Azimuth, numLayers, ZBed, Rh, Rv, traj);
	fastTIDriverObj.calcLogs();

	response[0] = fastTIDriverObj.m_vHfieldDyd[0].xx;
	response[1] = fastTIDriverObj.m_vHfieldDyd[0].xy;
	response[2] = fastTIDriverObj.m_vHfieldDyd[0].xz;
	response[3] = fastTIDriverObj.m_vHfieldDyd[0].yx;
	response[4] = fastTIDriverObj.m_vHfieldDyd[0].yy;
	response[5] = fastTIDriverObj.m_vHfieldDyd[0].yz;
	response[6] = fastTIDriverObj.m_vHfieldDyd[0].zx;
	response[7] = fastTIDriverObj.m_vHfieldDyd[0].zy;
	response[8] = fastTIDriverObj.m_vHfieldDyd[0].zz;
}

// The Periscope type of curves including Symmetrized, Anti-symmetrized, Harmonic resistivity, and Resistivity anisotropy
void directional_curve(DataType* curves, CFastTIDriver &fastTIDriverObj, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed,
	const DataType TVD, const DataType Dip, const DataType Azimuth, const DataType freq, const DataType sp) {
	
	complex<DataType> xx, xz, yy, zx, zz;
	DataType att_symm, ps_symm, att_anti, ps_anti, att_bulk, ps_bulk, att_anis, ps_anis;
	DataType traj[3] = {TVD,1,TVD};

	fastTIDriverObj.initToolConfig(freq, -sp/2, sp/2);
	fastTIDriverObj.initFormationPara(Dip, Azimuth, numLayers, ZBed, Rh, Rv, traj);
	fastTIDriverObj.calcLogs();

	xx = fastTIDriverObj.m_vHfieldDyd[0].xx;
	xz = fastTIDriverObj.m_vHfieldDyd[0].xz;
	yy = fastTIDriverObj.m_vHfieldDyd[0].yy;
	zx = fastTIDriverObj.m_vHfieldDyd[0].zx;
	zz = fastTIDriverObj.m_vHfieldDyd[0].zz;

	att_symm = -20 * log10( abs( ((zz + zx)*(zz - xz)) / ((zz - zx)*(zz + xz)) ) );
	ps_symm = 180 / PI*arg( ((zz + zx)*(zz - xz)) / ((zz - zx)*(zz + xz)) );
	att_anti = 20 * log10( abs( ((zz + zx)*(zz + xz)) / ((zz - zx)*(zz - xz)) ) );
	ps_anti = -180 / PI*arg( ((zz + zx)*(zz + xz)) / ((zz - zx)*(zz - xz)) );
	att_bulk = 20 * log10( abs( (xx + yy) / (2. * zz) ) );
	ps_bulk = 180 / PI*arg( -(xx + yy) / (2. * zz) );
	att_anis = 20 * log10( abs(xx / yy) );
	ps_anis = 180 / PI*arg( xx / yy );

	curves[0] = ps_symm; curves[1] = att_symm; curves[2] = ps_anti; curves[3] = att_anti;
	curves[4] = ps_bulk; curves[5] = att_bulk; curves[6] = ps_anis; curves[7] = att_anis;
}

// Conventional resistivity measurement (non-directional), adopted from Guidewave curve construction
void resistivity_curve(DataType* curves, CFastTIDriver &fastTIDriverObj, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed,
	const DataType TVD, const DataType Dip, const DataType Azimuth, const DataType freq, const DataType R1, const DataType R2) {
	
	complex<DataType> vr1, vr2;
	DataType att, ps;
	DataType traj[3] = {TVD,1,TVD};

	fastTIDriverObj.initToolConfig(freq, 0, R1);
	fastTIDriverObj.initFormationPara(Dip, Azimuth, numLayers, ZBed, Rh, Rv, traj);
	fastTIDriverObj.calcLogs();
	vr1 = fastTIDriverObj.m_vHfieldDyd[0].zz;
	
	fastTIDriverObj.initToolConfig(freq, 0, R2);
	fastTIDriverObj.calcLogs();
	vr2 = fastTIDriverObj.m_vHfieldDyd[0].zz;

	att = 20 * log10( abs(vr1) / abs(vr2) );
	ps = 180 / PI * ( arg(vr2) - arg(vr1) );

	curves[0] = ps; curves[1] = att;
}

// Lab simulation: 155/75/45ft, 2/6/12/24kHz, Periscope type; 24inch, 2MHz Guidewave type conventional resistivity curves. 42 curves in total.
void forward_chevron(out_ DataType* resp, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed, const DataType TVD, const DataType Dip) {
	DataType curves[8];
	CFastTIDriver fastTIDriverObj;
	CFastTIDriver &refObj = fastTIDriverObj;

	// Group 1: 155ft spacing, 2kHz frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e3, 1860);
	for (int i = 0; i < 8; i++)
		resp[i] = curves[i];
	// Group 2: 75ft spacing, 6/12kHz frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 6e3, 900);
	for (int i = 0; i < 8; i++)
		resp[i+8] = curves[i];
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 12e3, 900);
	for (int i = 0; i < 8; i++)
		resp[i+16] = curves[i];
	// Group 3: 45ft spacing, 12/24kHz frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 12e3, 540);
	for (int i = 0; i < 8; i++)
		resp[i+24] = curves[i];
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 540);
	for (int i = 0; i < 8; i++)
		resp[i+32] = curves[i];
	// Group 4: 19/25in spacing, 2MHz frequency, conventional
	resistivity_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e6, 19, 25);
	resp[40] = curves[0];
	resp[41] = curves[1];
}

// Real geosphere tool: one Tx, two Rx. Spacing: 80/34ft; frequency: 2/6/12kHz. Curve ordered by official manual. 48 curves in total.
void forward_geosphere(out_ DataType* resp, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed, const DataType TVD, const DataType Dip) {
	DataType curves[8];
	CFastTIDriver fastTIDriverObj;
	CFastTIDriver &refObj = fastTIDriverObj;
	DataType R1, R2, F1, F2, F3;
	R1 = 408.66; R2 = 954.33;
	F1 = 2e3; F2 = 6e3; F3 = 12e3;

	// Group 1: R1 spacing, F1 frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, F1, R1);
	resp[42] = curves[0]; resp[36] = curves[1];
	resp[6] = curves[2]; resp[0] = curves[3];
	resp[30] = curves[4]; resp[24] = curves[5];
	resp[18] = curves[6]; resp[12] = curves[7];

	// Group 2: R2 spacing, F1 frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, F1, R2);
	resp[43] = curves[0]; resp[37] = curves[1];
	resp[7] = curves[2]; resp[1] = curves[3];
	resp[31] = curves[4]; resp[25] = curves[5];
	resp[19] = curves[6]; resp[13] = curves[7];

	// Group 3: R1 spacing, F2 frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, F2, R1);
	resp[44] = curves[0]; resp[38] = curves[1];
	resp[8] = curves[2]; resp[2] = curves[3];
	resp[32] = curves[4]; resp[26] = curves[5];
	resp[20] = curves[6]; resp[14] = curves[7];

	// Group 4: R2 spacing, F2 frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, F2, R2);
	resp[45] = curves[0]; resp[39] = curves[1];
	resp[9] = curves[2]; resp[3] = curves[3];
	resp[33] = curves[4]; resp[27] = curves[5];
	resp[21] = curves[6]; resp[15] = curves[7];

	// Group 5: R1 spacing, F3 frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, F3, R1);
	resp[46] = curves[0]; resp[40] = curves[1];
	resp[10] = curves[2]; resp[4] = curves[3];
	resp[34] = curves[4]; resp[28] = curves[5];
	resp[22] = curves[6]; resp[16] = curves[7];

	// Group 6: R2 spacing, F3 frequency
	directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, F3, R2);
	resp[47] = curves[0]; resp[41] = curves[1];
	resp[11] = curves[2]; resp[5] = curves[3];
	resp[35] = curves[4]; resp[29] = curves[5];
	resp[23] = curves[6]; resp[17] = curves[7];
}

// Lab simulation: 150/80/40/20/10/4ft, 2/6/12/24/48/96/400/2000kHz, Periscope type; curves at each frequency can be selected.
void forward(out_ DataType* resp, const unsigned int numLayers, const DataType *Rh, DataType const *Rv, const DataType *ZBed, const DataType TVD, const DataType Dip, const int freq) {
	DataType curves[8];
	CFastTIDriver fastTIDriverObj;
	CFastTIDriver &refObj = fastTIDriverObj;

	switch (freq) {
	case 2:		// Frequency: 2kHz; Spacing: 150/80/40ft
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e3, 150 * 12);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e3, 80 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e3, 40 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	case 6:		// Frequency: 6kHz; Spacing: 80/40/28ft
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 6e3, 80 * 12);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 6e3, 40 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 6e3, 28 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	case 12:		// Frequency: 12kHz; Spacing: 40/28/20ft
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 12e3, 40 * 12);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 12e3, 28 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 12e3, 20 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	case 24:		// Frequency: 24kHz; Spacing: 28/20/16ft
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 28 * 12);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 20 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 16 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	case 48:		// Frequency: 48kHz; Spacing: 20/16/8ft
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 20 * 12);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 16 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 24e3, 8 * 12);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	case 100:		// Frequency: 100kHz; Spacing: 180/96/84inch
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 100e3, 180);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 100e3, 96);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 100e3, 84);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	case 400:		// Frequency: 400kHz; Spacing: 84/34/22inch
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 400e3, 84);
		for (int i = 0; i < 8; i++)
			resp[i] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 400e3, 34);
		for (int i = 0; i < 8; i++)
			resp[i + 8] = curves[i];
		directional_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 400e3, 22);
		for (int i = 0; i < 8; i++)
			resp[i + 16] = curves[i];
		break;
	default:	// resistivity measurement; frequency: 400kHz/2MHz
		resistivity_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 400e3, 37, 43);
		resp[0] = curves[0];
		resp[1] = curves[1];
		resistivity_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 400e3, 31, 37);
		resp[2] = curves[0];
		resp[3] = curves[1];
		resistivity_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e6, 25, 31);
		resp[4] = curves[0];
		resp[5] = curves[1];
		resistivity_curve(curves, refObj, numLayers, Rh, Rv, ZBed, TVD, Dip, 0, 2e6, 19, 25);
		resp[6] = curves[0];
		resp[7] = curves[1];
	}
}