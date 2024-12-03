//#pragma once
#ifndef COMMUTIL_H
#define COMMUTIL_H

#include <complex>
#include <vector>
#include <math.h>
#include <omp.h>

#include "Vector.h"
using namespace std;

// 
#define in_
#define out_
#define inout_

#ifndef NULL
#define NULL (0)
#endif

typedef double DataType;
typedef std::complex<DataType> Complex;

// typedef std::vector< DataType > DataType1D;
// typedef std::vector< std::vector< DataType > > DataType2D;

static const DataType PI = 3.141592653589793238462643;
static const DataType TWOPI = PI * 2.0;
static const DataType DEG_TO_RAD = PI/180.;

// defined light speed in vacuum [m/s]
static const DataType VLIGHT_0 = 2.99792458e+8; 

//permeability of free space [H/m]
static const DataType MU_0 = PI*4.0e-07;

//permittivity of free space [F/m]
static const DataType EPSILON_0 = 1.0/(VLIGHT_0*VLIGHT_0*MU_0);

// define unit conversion constant
static const DataType M_TO_METER = 1.0;
static const DataType DM_TO_METER = 0.1;
static const DataType CM_TO_METER = 0.01;
static const DataType MM_TO_METER = 0.001;
static const DataType INCH_TO_METER = 25.4/1000.;
static const DataType FOOT_TO_METER = 12.*25.4/1000;

// static const DataType ETA_0 = std::sqrt(MU_0/EPSILON_0);

static const double HANKX_E = 1.221402758160169834E0;
static const double HANKX_ER = 0.818730753077981859E0;
static const double HANKX_Y1 = 0.7358852661479794460E0;

static const Complex CPLX_ONE = Complex(1,0);
static const Complex CPLX_0 = Complex(0,0);
static const Complex CPLX_J = Complex(0,1);

typedef struct tagSToolConfig
{
	// factor from current unit to meter
	DataType s_FACTOR_T;

	// number of transmitters
	unsigned int s_iNum_Trans;
	
	// number of receivers
	unsigned int s_iNum_Recvs;

	// number of turns for transmitters
	std::vector<DataType> s_iNum_TransTurns;

	// number of turns for receivers
	std::vector<DataType> s_iNum_RecvsTurns;
	
	// distance to tool center for transmitters
	std::vector<DataType> s_vDistTrans;
	
	// distance to tool center for receivers
	std::vector<DataType> s_vDistRecvs;

	// maximum distance from transmitter to receiver
	DataType s_dMaxDist;

}SToolConfig;

typedef struct tagSFormationConfig
{
	// Dip angle
	DataType s_dAlfa;

	// Azimuthal angle
	DataType s_dBeta;

	// Orientiation angle
	DataType s_dGama;

	CDyadic<DataType> s_dyRotMatrixR;

	CDyadic<DataType> s_dyRotMatrixInvR;

	// factor from current unit to meter
	DataType s_FACTOR_T;

	// true for apparent depth, false for true depth
	bool s_bSign;

	// number of layers
	unsigned int s_iNum_Layers;

	// Horizontal Resistivity for each layer
	std::vector<DataType> s_vHResistivity;

	// Vertical Resistivity for each layer
	std::vector<DataType> s_vVResistivity;

	// Horizontal Relative Permittivity for each layer
	std::vector<DataType> s_vHRelEpsilon;

	// Vertical Relative Permittivity for each layer
	std::vector<DataType> s_vVRelEpsilon;

	// Relative Permeability for each layer
	std::vector<DataType> s_vHRelMu;
	std::vector<DataType> s_vVRelMu;

	// Relative Complex Horizontal Permittivity
	std::vector<Complex> s_vHCompEpsilon;

	// Relative Complex Vertical Permittivity
	std::vector<Complex> s_vVCompEpsilon;
	
	// Relative Complex Horizontal Permittivity
	std::vector<Complex> s_vHCompMu;

	// Relative Complex Vertical Permittivity
	std::vector<Complex> s_vVCompMu;

	std::vector<Complex> s_vTMAnsiRatio;
	std::vector<Complex> s_vTEAnsiRatio;

	std::vector<Complex> s_vjOmegaVMu;

	std::vector<Complex> s_vHWaveNum;

	// Layer Boundaries: s_vLayerBoundary[i] lower boundary of i-th layer
	std::vector<DataType> s_vLayerBoundary;

	// Layer Boundaries: s_vLayerDepth[i] depth of i-th layer
	std::vector<DataType> s_vLayerThick;

	// Static Fresnel Reflection Coefficient looking down;
	std::vector<Complex> s_vStaticFresnelDn[2];

	DataType s_K0;

}SFormationConfig;

typedef struct tagTxLineDataType
{
	/* i=0,1  i=0 for TM (e), i=1 for TE (h)  */

	// Wavenumber Kz for each layer
	std::vector<Complex> s_vWaveNumber_z[2];

	std::vector<Complex> s_vPhase[2];

	std::vector<Complex> s_vPhase2[2];

	// Characteristic impedance
	std::vector<Complex> s_vCharImpdZ[2];

	// Frensnel Reflection
	std::vector<Complex> s_vFrensnelRefCoef[2];

	// Reflection Looking UP
	std::vector<Complex> s_vRefCoefUp[2];

	// Reflection Looking DOWN
	std::vector<Complex> s_vRefCoefDown[2];

	std::vector<Complex> s_vTransCoefDown[2];

	inline void Init(const int NumLayer)
	{
		for (int iMode = 0; iMode < 2; iMode++)
		{
			this->s_vCharImpdZ[iMode].resize(NumLayer+2);
			this->s_vFrensnelRefCoef[iMode].resize(NumLayer+2);
			this->s_vPhase[iMode].resize(NumLayer+2);
			this->s_vPhase2[iMode].resize(NumLayer+2);
			this->s_vRefCoefDown[iMode].resize(NumLayer+2);
			this->s_vRefCoefUp[iMode].resize(NumLayer+2);
			this->s_vTransCoefDown[iMode].resize(NumLayer+2);
			this->s_vWaveNumber_z[iMode].resize(NumLayer+2);
		}
	}
	

}STxLineDataType;

typedef enum tagWellType
{ 

ENUM_VERTICAL = 0,

ENUM_DEVIATED = 1,

ENUM_HORIZONTAL = 2

}EWellType;

typedef struct tagSTxLineVI
{

	/* i=0,1  i=0 for TM (e), i=1 for TE (h)  */

	// Voltage due to current source
	Complex s_Vi[2];

	// Current due to current source
	Complex s_Ii[2];

	// Voltage due to voltage source
	Complex s_Vv[2];

	// Current due to voltage source
	Complex s_Iv[2];

	inline tagSTxLineVI& operator-=(const tagSTxLineVI &rhs)
	{
		for ( int i = 0; i<2; i++ )
		{
			this->s_Vi[i] = this->s_Vi[i] - rhs.s_Vi[i];
			this->s_Ii[i] = this->s_Ii[i] - rhs.s_Ii[i];
			this->s_Vv[i] = this->s_Vv[i] - rhs.s_Vv[i];
			this->s_Iv[i] = this->s_Iv[i] - rhs.s_Iv[i];
		}
		return (*this);
	}

}STxLineVI;

typedef struct tagSObsSrcPair
{

	DataType s_Obs_z;

	unsigned int s_iObsLayer;

	DataType s_Src_zPrim;

	unsigned int s_iSrcLayer;

	DataType s_rho;

	bool onTopFlag;

}SObsSrcPair;


static const DataType TOL = 1.E-7;

#define HTRANS_PTS 283
#define HTRANS_NUM 6

static const DataType HTRANS_W0[HTRANS_PTS] = 

	{ 2.1969101E-11,4.1201161E-09, -6.1322980E-09,7.2479291E-09,
     -7.9821627E-09,8.5778983E-09, -9.1157294E-09,9.6615250E-09,
     -1.0207546E-08,1.0796633E-08, -1.1393033E-08,1.2049873E-08,
     -1.2708789E-08,1.3446466E-08, -1.4174300E-08,1.5005577E-08,
     -1.5807160E-08,1.6747136E-08, -1.7625961E-08,1.8693427E-08,
     -1.9650840E-08,2.0869789E-08, -2.1903555E-08,2.3305308E-08,
     -2.4407377E-08,2.6033678E-08, -2.7186773E-08,2.9094334E-08,
     -3.0266804E-08,3.2534013E-08, -3.3672072E-08,3.6408936E-08,
     -3.7425022E-08,4.0787921E-08, -4.1543242E-08,4.5756842E-08,
     -4.6035233E-08,5.1425075E-08, -5.0893896E-08,5.7934897E-08,
     -5.6086570E-08,6.5475248E-08, -6.1539913E-08,7.4301996E-08,
     -6.7117043E-08,8.4767837E-08, -7.2583120E-08,9.7366568E-08,
     -7.7553611E-08,1.1279873E-07, -8.1416723E-08,1.3206914E-07,
     -8.3217217E-08,1.5663185E-07, -8.1482581E-08,1.8860593E-07,  
     -7.3963141E-08,2.3109673E-07, -5.7243707E-08,2.8867452E-07,
     -2.6163525E-08,3.6808773E-07,  2.7049871E-08,4.7932617E-07,
      1.1407365E-07,6.3720626E-07,  2.5241961E-07,8.6373487E-07,
      4.6831433E-07,1.1916346E-06,  8.0099716E-07,1.6696015E-06,
      1.3091334E-06,2.3701475E-06,  2.0803829E-06,3.4012978E-06,

      3.2456774E-06, 4.9240402E-06, 5.0005198E-06, 7.1783540E-06,
      7.6367633E-06, 1.0522038E-05, 1.1590021E-05, 1.5488635E-05,
      1.7510398E-05, 2.2873836E-05, 2.6368006E-05, 3.3864387E-05,
      3.9610390E-05, 5.0230379E-05, 5.9397373E-05, 7.4612122E-05,
      8.8951409E-05, 1.1094809E-04, 1.3308026E-04, 1.6511335E-04,
      1.9895671E-04, 2.4587195E-04, 2.9728181E-04, 3.6629770E-04,
      4.4402013E-04, 5.4589361E-04, 6.6298832E-04, 8.1375348E-04,
      9.8971624E-04, 1.2132772E-03, 1.4772052E-03, 1.8092022E-03,
      2.2045122E-03, 2.6980811E-03, 3.2895354E-03, 4.0238764E-03,
      4.9080203E-03, 6.0010999E-03, 7.3216878E-03, 8.9489225E-03,
      1.0919448E-02, 1.3340696E-02, 1.6276399E-02, 1.9873311E-02,
      2.4233627E-02, 2.9555699E-02, 3.5990069E-02, 4.3791529E-02,
      5.3150319E-02, 6.4341372E-02, 7.7506720E-02, 9.2749987E-02,
      1.0980561E-01, 1.2791555E-01, 1.4525830E-01, 1.5820085E-01,
      1.6058576E-01, 1.4196085E-01, 8.9781222E-02,-1.0238278E-02,
     -1.5083434E-01,-2.9059573E-01,-2.9105437E-01,-3.7973244E-02,
      3.8273717E-01, 2.2014118E-01,-4.7342635E-01, 1.9331133E-01,
      5.3839527E-02,-1.1909845E-01, 9.9317051E-02,-6.6152628E-02,
      4.0703241E-02,-2.4358316E-02, 1.4476533E-02,-8.6198067E-03,

      5.1597053E-03,-3.1074602E-03,1.8822342E-03,-1.1456545E-03,
      7.0004347E-04,-4.2904226E-04,2.6354444E-04,-1.6215439E-04,
      9.9891279E-05,-6.1589037E-05,3.7996921E-05,-2.3452250E-05,
      1.4479572E-05,-8.9417427E-06,5.5227518E-06,-3.4114252E-06,
      2.1074101E-06,-1.3019229E-06,8.0433617E-07,-4.9693681E-07,
      3.0702417E-07,-1.8969219E-07,1.1720069E-07,-7.2412496E-08,
      4.4740283E-08,-2.7643004E-08,1.7079403E-08,-1.0552634E-08,
      6.5200311E-09,-4.0284597E-09,2.4890232E-09,-1.5378695E-09,
      9.5019040E-10,-5.8708696E-10,3.6273937E-10,-2.2412348E-10,
      1.3847792E-10,-8.5560821E-11,5.2865474E-11,-3.2664392E-11,
      2.0182948E-11,-1.2470979E-11,7.7057678E-12,-4.7611713E-12,
      2.9415274E-12,-1.8170081E-12,1.1221034E-12,-6.9271067E-13,
      4.2739744E-13,-2.6344388E-13,1.6197105E-13,-9.9147443E-14,
      6.0487998E-14,-3.6973097E-14,2.2817964E-14,-1.4315547E-14,
      9.1574735E-15,-5.9567236E-15,3.9209969E-15,-2.5911739E-15,
      1.6406939E-15,-8.8248590E-16,3.0195409E-16, 2.2622634E-17,
     -8.0942556E-17,-3.7172363E-17,1.9299542E-16,-3.3388160E-16,
      4.6174116E-16,-5.8627358E-16,7.2227767E-16,-8.7972941E-16,
      1.0211793E-15,-1.0940039E-15,1.0789555E-15,-9.7089714E-16,

      7.4110927E-16,-4.1700094E-16, 8.5977184E-17, 1.3396469E-16,
     -1.7838410E-16, 4.8975421E-17, 1.9398153E-16,-5.0046989E-16,
      8.3280985E-16,-1.1544640E-15, 1.4401527E-15,-1.6637066E-15,
      1.7777129E-15,-1.7322187E-15, 1.5247247E-15,-1.1771155E-15,
      6.9747910E-16,-1.2088956E-16,-4.8382957E-16, 1.0408292E-15,
     -1.5220450E-15, 1.9541597E-15,-2.4107448E-15, 2.9241438E-15,
     -3.5176475E-15, 4.2276125E-15,-5.0977851E-15, 6.1428456E-15,
     -7.3949962E-15, 8.8597601E-15,-1.0515959E-14, 1.2264584E-14,
     -1.3949870E-14, 1.5332490E-14,-1.6146782E-14, 1.6084121E-14,
     -1.4962523E-14, 1.2794804E-14,-9.9286701E-15, 6.8825809E-15,
     -4.0056107E-15, 1.5965079E-15,-7.2732961E-18,-4.0433218E-16,
     -6.5679655E-16, 3.3011866E-15,-7.3545910E-15, 1.2394851E-14,
     -1.7947697E-14, 2.3774303E-14,-3.0279168E-14, 3.9252831E-14,
     -5.5510504E-14, 9.0505371E-14,-1.7064873E-13};


static const DataType HTRANS_W1[HTRANS_PTS] = {

	 -4.2129715E-16, 5.3667031E-15,-7.1183962E-15, 8.9478500E-15,
     -1.0767891E-14, 1.2362265E-14,-1.3371129E-14, 1.3284178E-14,
     -1.1714302E-14, 8.4134738E-15,-3.7726725E-15,-1.4263879E-15,
      6.1279163E-15,-9.1102765E-15, 9.9696405E-15,-9.3649955E-15,
      8.6009018E-15,-8.9749846E-15, 1.1153987E-14,-1.4914821E-14,
      1.9314024E-14,-2.3172388E-14, 2.5605477E-14,-2.6217555E-14,
      2.5057768E-14,-2.2485539E-14, 1.9022752E-14,-1.5198084E-14,
      1.1422464E-14,-7.9323958E-15, 4.8421406E-15,-2.1875032E-15,
     -3.2177842E-17, 1.8637565E-15,-3.3683643E-15, 4.6132219E-15,
     -5.6209538E-15, 6.4192841E-15,-6.8959928E-15, 6.9895792E-15,
     -6.5355935E-15, 5.6125163E-15,-4.1453931E-15, 2.6358827E-15,
     -9.5104370E-16, 1.4600474E-16, 5.6166519E-16, 8.2899246E-17,
      5.0032100E-16, 4.3752205E-16, 2.1052293E-15,-9.5451973E-16,
      6.4004437E-15,-2.1926177E-15, 1.1651003E-14, 5.8415433E-16,
      1.8044664E-14, 1.0755745E-14, 3.0159022E-14, 3.3506138E-14,
      5.8709354E-14, 8.1475200E-14, 1.2530006E-13, 1.8519112E-13,
      2.7641786E-13, 4.1330823E-13, 6.1506209E-13, 9.1921659E-13,
      1.3698462E-12, 2.0447427E-12, 3.0494477E-12, 4.5501001E-12,
      6.7870250E-12, 1.0126237E-11, 1.5104976E-11, 2.2536053E-11,

      3.3617368E-11, 5.0153839E-11, 7.4818173E-11, 1.1161804E-10,
      1.6651222E-10, 2.4840923E-10, 3.7058109E-10, 5.5284353E-10,
      8.2474468E-10, 1.2303750E-09, 1.8355034E-09, 2.7382502E-09,
      4.0849867E-09, 6.0940898E-09, 9.0913020E-09, 1.3562651E-08,
      2.0233058E-08, 3.0184244E-08, 4.5029477E-08, 6.7176304E-08,
      1.0021488E-07, 1.4950371E-07, 2.2303208E-07, 3.3272689E-07,
      4.9636623E-07, 7.4049804E-07, 1.1046805E-06, 1.6480103E-06,
      2.4585014E-06, 3.6677163E-06, 5.4714550E-06, 8.1626422E-06,
      1.2176782E-05, 1.8166179E-05, 2.7099223E-05, 4.0428804E-05,
      6.0307294E-05, 8.9971508E-05, 1.3420195E-04, 2.0021123E-04,
      2.9860417E-04, 4.4545291E-04, 6.6423156E-04, 9.9073275E-04,
      1.4767050E-03, 2.2016806E-03, 3.2788147E-03, 4.8837292E-03,
      7.2596811E-03, 1.0788355E-02, 1.5973323E-02, 2.3612041E-02,
      3.4655327E-02, 5.0608141E-02, 7.2827752E-02, 1.0337889E-01,
      1.4207357E-01, 1.8821315E-01, 2.2996815E-01, 2.5088500E-01,
      2.0334626E-01, 6.0665451E-02,-2.0275683E-01,-3.5772336E-01,
     -1.8280529E-01, 4.7014634E-01, 7.2991233E-03,-3.0614594E-01,
      2.4781735E-01,-1.1149185E-01, 2.5985386E-02, 1.0850279E-02,
     -2.2830217E-02, 2.4644647E-02,-2.2895284E-02, 2.0197032E-02,

     -1.7488968E-02,1.5057670E-02,-1.2953923E-02,1.1153254E-02,
     -9.6138436E-03,8.2952090E-03,-7.1628321E-03,6.1882910E-03,
     -5.3482055E-03,4.6232056E-03,-3.9970542E-03,3.4560118E-03,
     -2.9883670E-03,2.5840861E-03,-2.2345428E-03,1.9323046E-03,
     -1.6709583E-03,1.4449655E-03,-1.2495408E-03,1.0805480E-03,
     -9.3441130E-04,8.0803899E-04,-6.9875784E-04,6.0425624E-04,
     -5.2253532E-04,4.5186652E-04,-3.9075515E-04,3.3790861E-04,
     -2.9220916E-04,2.5269019E-04,-2.1851585E-04,1.8896332E-04,
     -1.6340753E-04,1.4130796E-04,-1.2219719E-04,1.0567099E-04,
     -9.1379828E-05,7.9021432E-05,-6.8334412E-05,5.9092726E-05,
     -5.1100905E-05,4.4189914E-05,-3.8213580E-05,3.3045496E-05,
     -2.8576356E-05,2.4711631E-05,-2.1369580E-05,1.8479514E-05,
     -1.5980307E-05,1.3819097E-05,-1.1950174E-05,1.0334008E-05,
     -8.9364160E-06,7.7278366E-06,-6.6827083E-06,5.7789251E-06,
     -4.9973715E-06,4.3215167E-06,-3.7370660E-06,3.2316575E-06,
     -2.7946015E-06,2.4166539E-06,-2.0898207E-06,1.8071890E-06,
     -1.5627811E-06,1.3514274E-06,-1.1686576E-06,1.0106059E-06,
     -8.7392952E-07,7.5573750E-07,-6.5353002E-07,5.6514528E-07,
     -4.8871388E-07,4.2261921E-07,-3.6546333E-07,3.1603732E-07

     -2.7329579E-07,2.3633470E-07,-2.0437231E-07,1.7673258E-07,
     -1.5283091E-07,1.3216174E-07,-1.1428792E-07,9.8831386E-08,
     -8.5465227E-08,7.3906734E-08,-6.3911437E-08,5.5267923E-08,
     -4.7793376E-08,4.1329702E-08,-3.5740189E-08,3.0906612E-08,
     -2.6726739E-08,2.3112160E-08,-1.9986424E-08,1.7283419E-08,
     -1.4945974E-08,1.2924650E-08,-1.1176694E-08,9.6651347E-09,
     -8.3580023E-09,7.2276490E-09,-6.2501673E-09,5.4048822E-09,
     -4.6739154E-09,4.0418061E-09,-3.4951847E-09,3.0224895E-09,
     -2.6137226E-09,2.2602382E-09,-1.9545596E-09,1.6902214E-09,
     -1.4616324E-09,1.2639577E-09,-1.0930164E-09,9.4519327E-10,
     -8.1736202E-10,7.0681930E-10,-6.1122713E-10,5.2856342E-10,
     -4.5707937E-10,3.9526267E-10,-3.4180569E-10,2.9557785E-10,
     -2.5560176E-10,2.2103233E-10,-1.9113891E-10,1.6528994E-10,
     -1.4294012E-10,1.2361991E-10,-8.2740936E-11};


inline DataType SIGN( DataType a, DataType b )
{
	DataType tmp = a-b;
	if (tmp >= 0)
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}

#endif