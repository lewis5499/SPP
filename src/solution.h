/*------------------------------------------------------------------------------
* solution.h : SPP software common functions
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references : [1]"RTK_Structs.h"
*              [2] rtklib_2.4.3_b34
*
* version : $Revision: 1.1 $ $Date: 2023/10/10 15:37:41 $
* history : 2023/10/10 1.0 new
*
*-----------------------------------------------------------------------------*/
#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "strclib.h"
#include <string>

/* Conversion function between universal time, GPS time and simplified Julian day------------------*/
bool CT2MJD(const COMMONTIME *CT, MJDTIME *MJDT);
bool MJD2CT(const MJDTIME *MJDT, COMMONTIME *CT);
bool GPST2MJD(const GPSTIME *GT, MJDTIME *MJDT);
bool MJD2GPST(const MJDTIME *MJDT, GPSTIME *GT);
bool CT2GPST(const COMMONTIME *CT, GPSTIME *GT);
bool GPST2CT(const GPSTIME *GT, COMMONTIME *CT);
double TimeDiff(const GPSTIME *GT2, const GPSTIME*GT1);
const char* time2str(const GPSTIME *GT);
std::string getCurrentTimeString();
std::string convertToStdString(const char* input);

/* Conversion between deg and rad------------------------------------------------------------------*/
void deg2dms(const double deg, double *dms, int ndec);
double dms2deg(const double *dms);

/* Conversions: ecef Pos, geodetic Pos, local coordinate Pos---------------------------------------*/
void blh2enuTransMat(const double *blhPos, double *e2nMat);
void xyz2blh(const double *xyzPos, double *blhPos);
void blh2xyz(const double *blhPos, double *xyzPos);
void xyz2enu(const double *xyzPos, const double *refxyzPos, double *enuPos);
void compSatElAz(const double *Xr, const double *Xs, double *Elev, double *Azim); 
void denuPos(const double *X0, const double *Xr, double *dNeu);  
double geodist(const double *rs, const double *rr, double *e);

/* satellite common functions----------------------------------------------------------------------*/
double svClkCorr(const double dt, const double *params);
double svClkdCorr(const double dt, const double *params);
double *earthRotCorr(const double *vec3, const double time, const double omga);
double getPIF(EPOCHOBS *obs, const int idx, const int sys);

/* outlier detection-------------------------------------------------------------------------------*/
int markOutlier(EPOCHOBS *obs);
int initMWGF(EPOCHOBS *obs, const int idx);
void detectOutlier(EPOCHOBS *obs, MWGF *lastEpoch);
void deepCopy(const MWGF *source, MWGF *dst);

/* buff deep copy----------------------------------------------------------------------------------*/
int bread(int8_t *oribuff, uint8_t *dstbuff, int copiedlength);

/* make the obs in order --------------------------------------------------------------------------*/
int inOrder(EPOCHOBS *obs);

/* Tropospheric correction ------------------------------------------------------------------------*/
double Hopfield(const double H, const double Elev);
void tropCorr_getEl(EPOCHOBS *obs, const double *Xr, double *BLHr, const int *locs, const int num);

/* compute satPos and satClk Corr -----------------------------------------------------------------*/
int compSatPosVel(GPSEPHREC *ephgps, BDSEPHREC *ephbds, EPOCHOBS *obs, const bool isSPP);
int gpsPosVel(GPSEPHREC *ephgps, EPOCHOBS *obs, const int prn_1, const bool isSPP, const double dt);
int bdsPosVel(BDSEPHREC *ephbds, EPOCHOBS *obs, const int prn_1, const bool isSPP, const double dt);
int findSatidx(EPOCHOBS *obs, const short prn, const int sys);

/* spp: precise SatClockOft/Sft, the satPos with Earth rotation Correction ------------------------*/
int satClk_EarthRotCorr(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, const double *Xr);
int gpsSatClk(EPOCHOBS *obs, GPSEPHREC *ephgps, const int idx, const double dt);
int bdsSatClk(EPOCHOBS *obs, BDSEPHREC *ephbds, const int idx, const double dt);

/* spp-spv ----------------------------------------------------------------------------------------*/
int SPP_SPV(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *pos);
int spp_LS_Epoch(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *pos);
int spv_LS_Epoch(EPOCHOBS *obs, POSRES *pos);
int spp_KF_Epoch(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *pos);
int sppfindValidSat_locs(EPOCHOBS *obs, int *locs, int &num, int &numgps, int &numbds);
int spvfindValidSat_locs(EPOCHOBS *obs, int *locs, int &num);

/* Least Squares ----------------------------------------------------------------------------------*/
int spp_LS(EPOCHOBS *obs, const int *locs, double *Xr, double *rcvClkoft, double *Qxx, double *Dxx,
	int necessNum, int num, double &sigmaPos, double &Pdop, double *v);
int spv_LS(EPOCHOBS *obs, const int *locs, const int num, const int necessNum, double *vel, double *Xr,
	double &sigmaVel, double &rcvClksft);

/* Extended Kalman Filter -------------------------------------------------------------------------*/
int smoother(const double *xf, const double *Qf, const double *xb, const double *Qb, int n, double *xs, double *Qs);
int oneStep_Prediction(const double *Xk_1, const double *Dxk_1, double *Xkk_1, double *Dxkk_1);
int measurement_Update(const double *Dxkk_1, const double *Xkk_1, double *Xk, double *Dxk, EPOCHOBS *obs, const int *locs, int num);

#endif

