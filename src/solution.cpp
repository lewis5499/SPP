#include "solution.h"
#include "consts.h"
#include "matrix.h"
#include <cmath>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

using namespace matrix;

/* Convert Common Time to MJD --------------------------------------------------
* convert Common Time to MJD
* args   : const COMMONTIME * CT    I   to be converted
*          MJDTIME * MJDT           O   target TIME
* return : bool				        O   success or failure
*-----------------------------------------------------------------------------*/
bool CT2MJD(const COMMONTIME *CT, MJDTIME *MJDT){
	int y, m = 0;
	//First determine whether CT is initialized/available.
	//If the month is January or February, subtract 1 from the year and add 12 to the month.
	if (CT->Month == 0 || CT->Day == 0){
		return false;
	}else if (CT->Month <= 2){
		y = CT->Year - 1;
		m = CT->Month + 12;
	}else{
		y = CT->Year;
		m = CT->Month;
	}                                                                        //30.6001:Avoid loss of computer accuracy
	MJDT->Days = (int)(365.25*y) + (int)(30.6001*(m + 1)) + int(CT->Day) - 679019;//1720981.5-2400000.5
	MJDT->FracDay = (CT->Hour + CT->Minute / 60.0 + CT->Second / 3600.0) / 24.0;
	return true;
}

/* Convert MJD to Common Time --------------------------------------------------
* convert MJD to Common Time
* args   : const MJDTIME * MJDT     I   to be converted
*          COMMONTIME * CT          O   target TIME
* return : bool				        O   success or failure
*-----------------------------------------------------------------------------*/
bool MJD2CT(const MJDTIME *MJDT, COMMONTIME *CT){
	if (MJDT->Days == 0 && abs(MJDT->FracDay) < EPSILON){
		return false;
	}            //2400000.5+0.5, a is the integer part
	int a = (int)(MJDT->Days + MJDT->FracDay + 2400001);
	int b = a + 1537;
	int c = (int)((b - 122.1) / 365.25);
	int d = (int)(365.25*c);
	int e = (int)((b - d) / 30.6001);

	CT->Day = b - d - (int)(30.6001*e);
	CT->Month = e - 1 - 12 * (int)(e / 14);
	CT->Year = c - 4715 - (int)((7 + CT->Month) / 10);

	CT->Hour = (int)(MJDT->FracDay * 24);
	CT->Minute = (int)((MJDT->FracDay * 24 - CT->Hour) * 60);
	CT->Second = ((MJDT->FracDay * 24 - CT->Hour) * 60 - CT->Minute) * 60.0;
	return true;
}

/* Convert GPST to MJD ---------------------------------------------------------
* convert GPS Time to MJD
* args   : const GPSTIME * GT       I   to be converted
*          MJDTIME * MJDT           O   target TIME
* return : bool				        O   success or failure
*-----------------------------------------------------------------------------*/
bool GPST2MJD(const GPSTIME *GT, MJDTIME *MJDT){
	if (GT->Week == 0 && abs(GT->SecOfWeek) < EPSILON){
		return false;
	}
	MJDT->Days = 44244 + GT->Week * 7 + int(GT->SecOfWeek / 86400.0);
	MJDT->FracDay = 44244 + GT->Week * 7 + (GT->SecOfWeek / 86400.0) - MJDT->Days;
	return true;
}

/* Convert MJD to GPST ---------------------------------------------------------
* convert MJD to GPS Time
* args   : const MJDTIME * MJDT     I   to be converted
*          GPSTIME * GT             O   target TIME
* return : bool				        O   success or failure
*-----------------------------------------------------------------------------*/
bool MJD2GPST(const MJDTIME *MJDT, GPSTIME *GT){
	if (MJDT->Days == 0 && abs(MJDT->FracDay) < EPSILON){
		return false;
	}
	GT->Week = (int)((MJDT->Days + MJDT->FracDay - 44244) / 7.0);
	GT->SecOfWeek = (MJDT->Days + MJDT->FracDay - 44244 - GT->Week * 7) * 86400;
	return true;
}

/* Convert Common Time to GPST -------------------------------------------------
* convert Common Time to GPS Time 
* args   : const COMMONTIME * CT    I   to be converted
*          GPSTIME * GT             O   target TIME
* return : bool				        O   success or failure
*-----------------------------------------------------------------------------*/
bool CT2GPST(const COMMONTIME *CT, GPSTIME *GT){
	if (CT->Month == 0 || CT->Day == 0){
		return false;
	}
	MJDTIME *MJD = new MJDTIME;
	CT2MJD(CT, MJD);
	MJD2GPST(MJD, GT);
	delete MJD;
	return true;
}

/* Convert GPST to Common Time -------------------------------------------------
* convert GPS Time to Common Time
* args   : const GPSTIME * GT    I   to be converted
*          COMMONTIME * CT       O   target TIME
* return : bool				     O   success or failure
*-----------------------------------------------------------------------------*/
bool GPST2CT(const GPSTIME *GT, COMMONTIME *CT){
	if (GT->Week == 0 && abs(GT->SecOfWeek) < EPSILON){
		return false;
	}
	MJDTIME *MJD = new MJDTIME;
	GPST2MJD(GT, MJD);
	MJD2CT(MJD, CT);
	delete MJD;
	return true;
}

/* get difference of two GPST time ---------------------------------------------
* convert degree to degree-minute-second
* args   : const GPSTIME * GT2    I   GPST: GT2 subtract
*          const GPSTIME * GT1    I   GPST: GT1 subtracted
* return : double				  O   DiffGPST(seconds):GT2-GT1
*-----------------------------------------------------------------------------*/
double TimeDiff(const GPSTIME *GT2, const GPSTIME *GT1){
	//RETURN=SECONDS(GT2-GT1)
	return (GT2->Week - GT1->Week) * 7 * 86400.0 + (GT2->SecOfWeek - GT1->SecOfWeek);
}

/* get const char* type of GPST time ---------------------------------------------
* convert GPST time to const char*
* args   : const GPSTIME * GT    I   GPST: GT
* return : const char*           O   the fisrt address of the GPST str
*-----------------------------------------------------------------------------*/
const char* time2str(const GPSTIME *GT)
{
	string *str = new string("GPS Week: " + to_string(GT->Week) + "  Tow: " + to_string(GT->SecOfWeek) + "\n");
	return (*str).c_str();
}

/* get Current Time String ------------------------------------------------------
* get Current Time
* args   : NONE
* return : string				  O   the Current Time formatted as"%Y-%m-%d-%H-%M-%S"
*-----------------------------------------------------------------------------*/
std::string getCurrentTimeString()
{
	auto now = std::chrono::system_clock::now();
	std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
	std::tm* localTime = std::localtime(&currentTime);
	std::stringstream ss;
	ss << std::put_time(localTime, "%Y-%m-%d-%H-%M-%S");
	return ss.str();
}

/* convertToStdString Function --------------------------------------------------
 * Converts a C-style string (const char*) to a C++ std::string.
 * args
 *     input: const char* - The C-style string to be converted to std::string.
 * return
 *     std::string - The converted string in std::string format.
 * Notes
 *     If the input pointer is nullptr, the function returns an empty std::string.
 *-----------------------------------------------------------------------------*/
std::string convertToStdString(const char* input) {
	if (input == nullptr) {
		return std::string(); 
	} else {
		return std::string(input);
	}
}

/* convert degree to deg-min-sec -----------------------------------------------
* convert degree to degree-minute-second
* args   : const double deg       I   degree
*          double *dms            O   degree-minute-second {deg,min,sec}
*          int    ndec            I   number of decimals of second
* return : none
*-----------------------------------------------------------------------------*/
void deg2dms(const double deg, double *dms, int ndec){
	double sign = deg < 0.0 ? -1.0 : 1.0, a = fabs(deg);
	double unit = pow(0.1, ndec);
	dms[0] = floor(a); a = (a - dms[0])*60.0;
	dms[1] = floor(a); a = (a - dms[1])*60.0;
	dms[2] = floor(a / unit + 0.5)*unit;
	if (dms[2] >= 60.0){
		dms[2] = 0.0;
		dms[1] += 1.0;
		if (dms[1] >= 60.0){
			dms[1] = 0.0;
			dms[0] += 1.0;
		}
	}
	dms[0] *= sign;
}

/* convert deg-min-sec to degree -----------------------------------------------
* convert degree-minute-second to degree
* args   : const double *dms      I   degree-minute-second {deg,min,sec}
* return : double                 O   degree
*-----------------------------------------------------------------------------*/
double dms2deg(const double *dms){
	double sign = dms[0] < 0.0 ? -1.0 : 1.0;
	return sign * (fabs(dms[0]) + dms[1] / 60.0 + dms[2] / 3600.0);
}

/* geodetic position to local coordinate transfromation matrix -----------------
* compute geodetic position to local coordinate transfromation matrix
* args   : const double *blhPos     I   geodetic position {lat,lon} (rad)
*          double *e2nMat           O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by row-major order
*-----------------------------------------------------------------------------*/
void blh2enuTransMat(const double *refblhPos, double *e2nMat){
	double sinb = sin(refblhPos[0]), cosb = cos(refblhPos[0]), sinl = sin(refblhPos[1]), cosl = cos(refblhPos[1]);

	e2nMat[0] = -sinl;        e2nMat[1] = cosl;         e2nMat[2] = 0.0;
	e2nMat[3] = -sinb * cosl; e2nMat[4] = -sinb * sinl; e2nMat[5] = cosb;
	e2nMat[6] = cosb * cosl;  e2nMat[7] = cosb * sinl;  e2nMat[8] = sinb;
}

/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : const double *xyzPos    I   ecef position {x,y,z} (m)
*          double *blhPos          O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void xyz2blh(const double *xyzPos, double *blhPos){
	double e2 = F_WGS84 * (2.0 - F_WGS84), z, zk, n = R_WGS84, sinb;
	double r2 = xyzPos[0] * xyzPos[0] + xyzPos[1] * xyzPos[1];
	for (z = xyzPos[2], zk = 0.0; fabs(z - zk) >= 1E-5;){
		zk = z;
		sinb = z / sqrt(r2 + z * z);
		n = R_WGS84 / sqrt(1.0 - e2 * sinb*sinb);
		z = xyzPos[2] + n * e2*sinb;
	}// latitude is -90 to 90 degrees, we use atan();
	 // longitude is -180 to 180 degrees, we use atan2();
	 // additionally, it is necessary to distinguish the North and South Poles here.
	blhPos[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (xyzPos[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
	blhPos[1] = r2 > 1E-12 ? atan2(xyzPos[1], xyzPos[0]) : 0.0;
	blhPos[2] = sqrt(r2 + z * z) - n;
}

/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : const double *blhPos    I   geodetic position {lat,lon,h} (rad,m)
*          double *xyzPos          O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void blh2xyz(const double *blhPos, double *xyzPos){
	double sinb = sin(blhPos[0]), cosb = cos(blhPos[0]), sinl = sin(blhPos[1]), cosl = cos(blhPos[1]);
	double e2 = F_WGS84 * (2.0 - F_WGS84), n = R_WGS84 / sqrt(1.0 - e2 * sinb*sinb);

	xyzPos[0] = (n + blhPos[2])*cosb*cosl;
	xyzPos[1] = (n + blhPos[2])*cosb*sinl;
	xyzPos[2] = (n*(1.0 - e2) + blhPos[2])*sinb;
}

/* transfrom ecef position to local coordinate ---------------------------------
* compute ecef position to local coordinate transfromation matrix
* args   : const double *xyzPos        I   ecef position {x,y,z} (m)
*          const double *refxyzPos     I   reference ecef position {x0,y0,z0} (m)
*          double       *enuPos        O   local coordinate position {e,n,u} (m)
* return : none
*-----------------------------------------------------------------------------*/
void xyz2enu(const double *xyzPos, const double *refxyzPos, double *enuPos){
	/*double refblhPos[3] = { 0 }; xyz2blh(refxyzPos, refblhPos);
	double sinL, cosL, sinB, cosB;
	double dx = xyzPos[0] - refxyzPos[0];
	double dy = xyzPos[1] - refxyzPos[1];
	double dz = xyzPos[2] - refxyzPos[2];
	sinB = sin(refblhPos[0]), cosB = cos(refblhPos[0]);
	sinL = sin(refblhPos[1]), cosL = cos(refblhPos[1]);
	
	enuPos[0] = -sinL * dx + cosL * dy;
	enuPos[1] = -sinB * cosL*dx - sinB * sinL*dy + cosB * dz;
	enuPos[2] = cosB * cosL*dx + cosB * sinL*dy + sinB * dz;*/

	double dxyzPos[3], refblhPos[3], e2nmat[9];
	xyz2blh(refxyzPos, refblhPos);
	matSub(xyzPos, refxyzPos, dxyzPos, 3, 1);
	blh2enuTransMat(refblhPos, e2nmat);
	matMul(e2nmat, 3, 3, dxyzPos, 3, 1, enuPos);
}

/* Compute Satellite Elevation angle and Azimuth angle -------------------------
* compute satellite elevation angle and ezimuth angle based on sat Pos and base Pos
* args   : const double *Xr   I   base Pos: ecef position {x,y,z} (m)
*          const double *Xs   I   sat Pos:  ecef position {x,y,z} (m)
*		   double *Elev       O   satellite elevation angle relative to the base station
*          double *Azim       O   satellite azimuth angle relative to the base station
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void compSatElAz(const double *Xr, const double *Xs, double *Elev, double *Azim){
	double XrblhPos[3] = { 0,0,0 }, dxyzMat[3] = { 0,0,0 }, denuMat[3] = { 0,0,0 };
	double Xr2enuMat[9] = { 0,0,0,0,0,0,0,0,0 };
	xyz2blh(Xr, XrblhPos);
	blh2enuTransMat(XrblhPos, Xr2enuMat);
	matSub(Xs, Xr, dxyzMat, 3, 1);
	matMul(Xr2enuMat, 3, 3, dxyzMat, 3, 1, denuMat);

	double r2 = denuMat[1] * denuMat[1] + denuMat[0] * denuMat[0];
	*Elev = r2 > 1E-12 ? atan(denuMat[2] / sqrt(r2)) : (denuMat[2] > 0.0 ?  PI / 2.0 : -PI / 2.0);
	*Azim = atan2(denuMat[0], denuMat[1]);
}

/* Compute the Position error in the local coordinate -------------------------
* compute the position error in the local coordinate: enu, based on the known precise base station Pos
* args   : const double *X0      I   precise base Pos: ecef position {x,y,z} (m)
*          const double *Xr      I   computed base Pos: ecef position {x,y,z} (m)
*		   double *denuErrorMat  O   the pos Error of base station: enu position {dE, dN, dU} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
void denuPos(const double *X0, const double *Xr, double *denu){
	double X0blhPos[3] = { 0 }, X02enuMat[9] = { 0 }, dxyzMat[3] = { 0 };
	xyz2blh(X0, X0blhPos);
	blh2enuTransMat(X0blhPos, X02enuMat);
	matSub(Xr, X0, dxyzMat, 3, 1);
	matMul(X02enuMat, 3, 3, dxyzMat, 3, 1, denu);
}

/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
double geodist(const double *rs, const double *rr, double *e) {
	double r;
	int i;

	if (norm(rs, 3) < R_WGS84) return -1.0;
	for (i = 0; i < 3; i++) { e[i] = rs[i] - rr[i]; }
	r = norm(e, 3);
	for (i = 0; i < 3; i++) { e[i] /= r; }
	return r + Omega_WGS * (rs[0] * rr[1] - rs[1] * rr[0]) / CLIGHT;
}

/* get sat Clock bias correction */
double svClkCorr(const double dt, const double *params) {
	return params[0] + params[1]*dt + params[2]*SQR(dt);
}

/* get sat Clock drift correction */
double svClkdCorr(const double dt, const double *params) {
	return params[1] + 2 * params[2] * dt;
}

/* get earth rotation correction matrix */
double *earthRotCorr(const double *vec3, const double time, const double omga)
{
	double rot[9];
	rot[0] =  cos(omga*time); rot[1] = sin(omga*time); rot[2] = 0;
	rot[3] = -sin(omga*time); rot[4] = cos(omga*time); rot[5] = 0;
	rot[6] = 0;				  rot[7] = 0;			   rot[8] = 1;
	double *ivec = new double[3];
	matMul(rot, 3, 3, vec3, 3, 1, ivec);
	return ivec;
}

/* get Hopfield model correction */
double Hopfield(const double H, const double Elev)
{
	if (H<-10000.0||H>8848.86) { return 0.0; }
	double E = Elev * Deg;
	double RH = RH0 * exp(-0.0006396*(H - H0));
	double p = p0 * pow(1 - 0.0000226*(H - H0), 5.225);
	double T = T0 - 0.0065*(H - H0);
	double e = RH * exp(-37.2465 + 0.213166*T - 0.000256908*SQR(T));
	double hw = 11000.0;
	double hd = 40136.0 + 148.72*(T0 - 273.16);
	double Kw = 155.2*(1e-7) * 4810.0 / (SQR(T))*e*(hw - H);
	double Kd = 155.2*(1e-7)*p / T * (hd - H);

	return Kd / sin(sqrt(SQR(E) + 6.25)*Rad) + Kw / sin(sqrt(SQR(E) + 2.25)*Rad);
}

/* buff deep copy */
int bread(int8_t *oribuff, uint8_t *dstbuff, int copiedlength)
{
	int idx;
	for (idx = 0; idx < copiedlength;idx++){
		dstbuff[idx] = oribuff[idx];
	}
	return 0;
}


/* smoother --------------------------------------------------------------------
* combine forward and backward filters by fixed-interval smoother as follows:
*
*   xs=Qs*(Qf^-1*xf+Qb^-1*xb), Qs=(Qf^-1+Qb^-1)^-1)
*
* args   : double *xf       I   forward solutions (n x 1)
* args   : double *Qf       I   forward solutions covariance matrix (n x n)
*          double *xb       I   backward solutions (n x 1)
*          double *Qb       I   backward solutions covariance matrix (n x n)
*          int    n         I   number of solutions
*          double *xs       O   smoothed solutions (n x 1)
*          double *Qs       O   smoothed solutions covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : see reference [4] 5.2
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
int smoother(const double *xf, const double *Qf, const double *xb,
	const double *Qb, int n, double *xs, double *Qs)
{
	double *invQf = mat(n, n), *invQb = mat(n, n), *xx = mat(n, 1);
	int i, info = -1;

	matcpy(invQf, Qf, n, n);
	matcpy(invQb, Qb, n, n);
	if (!matinv(invQf, n) && !matinv(invQb, n)) {
		for (i = 0; i < n*n; i++) Qs[i] = invQf[i] + invQb[i];
		if (!(info = matinv(Qs, n))) {
			matmul("NN", n, 1, n, 1.0, invQf, xf, 0.0, xx);
			matmul("NN", n, 1, n, 1.0, invQb, xb, 1.0, xx);
			matmul("NN", n, 1, n, 1.0, Qs, xx, 0.0, xs);
		}
	}
	free(invQf); free(invQb); free(xx);
	return info;
}

int oneStep_Prediction(const double *Xk_1, const double *Dxk_1, double *Xkk_1, double *Dxkk_1)
{
	const double dt = 1.0;

	/* state transfer matrix */
	double STM[36];
	eye(STM, 6);
	STM[23] = dt;
	STM[29] = dt;
	double Trans_STM[36];
	matTran(STM, 6, 6, Trans_STM);
	
	/* state noise matrix */
	double Dwk_1[36];
	for (int i = 0; i < 36; i++) { Dwk_1[i] = 0.0; }
	Dwk_1[0]  = SQR(sigma_x)*dt;
	Dwk_1[7]  = SQR(sigma_y)*dt;
	Dwk_1[14] = SQR(sigma_z)*dt;
	Dwk_1[21] = SQR(sigma_bG)*dt + 1.0 / 3 * SQR(sigma_d)*dt*dt*dt;
	Dwk_1[22] = 1.0 / 3 * SQR(sigma_d)*dt*dt*dt;
	Dwk_1[23] = 1.0 / 2 * SQR(sigma_d)*dt*dt;
	Dwk_1[27] = 1.0 / 3 * SQR(sigma_d)*dt*dt*dt;
	Dwk_1[28] = SQR(sigma_bC)*dt + 1.0 / 3 * SQR(sigma_d)*dt*dt*dt;
	Dwk_1[29] = 1.0 / 2 * SQR(sigma_d)*dt*dt;
	Dwk_1[33] = 1.0 / 2 * SQR(sigma_d)*dt*dt;
	Dwk_1[34] = 1.0 / 2 * SQR(sigma_d)*dt*dt;
	Dwk_1[35] = SQR(sigma_d)*dt;

	/* ahead prediction */
	matMul(STM, 6, 6, Xk_1, 6, 1, Xkk_1);
	double temp1[36], temp2[36];
	matMul(STM, 6, 6, Dxk_1, 6, 6, temp1);
	matMul(temp1, 6, 6, Trans_STM, 6, 6, temp2);

	matAdd(temp2, Dwk_1, Dxkk_1, 6, 6);

	return 0;
}

int measurement_Update(const double *Dxkk_1, const double *Xkk_1, double *Xk, double *Dxk, EPOCHOBS *obs, const int *locs, int num)
{
	const int MAXNUM = 50;

	/* measurement noise matrix */
	double Dr[MAXNUM * MAXNUM], Dr_copy[MAXNUM * MAXNUM];
	eye(Dr, num); 
	eye(Dr_copy, num);
	scalar_matMul(var_obs,Dr, Dr, num*num);
	scalar_matMul(var_obs,Dr_copy, Dr_copy, num*num);

	/* design matrix & predicted residual */
	double H[MAXNUM * 6];
	double d0, dX[3];
	int idx;
	double V[MAXNUM];
	double F;
	for (int i = 0; i < num; i++)
	{
		idx = locs[i];
		matSub(Xkk_1, obs->SatPVT[idx].SatPos, dX, 3, 1);
		d0 = norm(dX, 3);

		H[i * 6 + 0] = dX[0] / d0;
		H[i * 6 + 1] = dX[1] / d0;
		H[i * 6 + 2] = dX[2] / d0;
		if (obs->SatPVT[idx].Sys == GPS) {
			H[i * 6 + 3] = 1.0;
			H[i * 6 + 4] = 0.0;
			F = d0 + Xkk_1[3] - CLIGHT * obs->SatPVT[idx].SatClkOft + obs->SatPVT[idx].TropCorr;
			V[i] = getPIF(obs, idx, GPS) - F;
		}else{
			H[i * 6 + 3] = 0.0;
			H[i * 6 + 4] = 1.0;
			F = d0 + Xkk_1[4] - CLIGHT * obs->SatPVT[idx].SatClkOft + obs->SatPVT[idx].TropCorr;
			V[i] = getPIF(obs, idx, BDS) - F;
		}
		H[i * 6 + 5] = 0.0;
	}

	/* gain matrix */
	double K[6 * MAXNUM];
	double temp1[6 * MAXNUM];
	double Trans_H[6 * MAXNUM];
	matTran(H, num, 6, Trans_H);
	matMul(Dxkk_1, 6, 6, Trans_H, 6, num, temp1);

	double temp2[MAXNUM * 6];
	matMul(H, num, 6, Dxkk_1, 6, 6, temp2);
	double temp3[MAXNUM * MAXNUM];
	matMul(temp2, num, 6, Trans_H, 6, num, temp3);
	matAdd(temp3, Dr, Dr, num, num);
	matinv(Dr, num);

	matMul(temp1, 6, num, Dr, num, num, K);

	/* Xk */
	double KV[6];
	matMul(K, 6, num, V, num, 1, KV);

	matAdd(Xkk_1, KV, Xk, 6, 1);

	/* Dxk */
	double I[36], KH[36], T[36];
	eye(I, 6);
	matMul(K, 6, num, H, num, 6, KH);
	matSub(I, KH, I, 6, 6);
	matTran(I, 6, 6, T);
	double temp4[36], temp5[36];
	matMul(I, 6, 6, Dxkk_1, 6, 6, temp4);
	matMul(temp4, 6, 6, T, 6, 6, temp5);

	double temp6[6 * MAXNUM];
	double Trans_K[MAXNUM * 6];
	double temp7[36];
	matTran(K, 6, num, Trans_K);
	matMul(K, 6, num, Dr_copy, num, num, temp6);
	matMul(temp6, 6, num, Trans_K, num, 6, temp7);

	matAdd(temp5, temp7, Dxk, 6, 6);

	return 0;
}

/* Get the PIF of 'obs.SatObs[idx]' according to sys (TGD Correction included)------------------ */
double getPIF(EPOCHOBS *obs, const int idx, const int sys)
{
	double p = 0.0;
	if (sys == BDS) {
		double p1, p3;
		p1 = obs->SatObs[idx].P[0]; p3 = obs->SatObs[idx].P[1];
		p = (p3 - p1 * SQR(FC13R) + CLIGHT * SQR(FC13R)*obs->SatPVT[idx].Tgd1) / (1 - SQR(FC13R));
	} else if (sys == GPS) {
		p = obs->SatPVT[idx].PIF;
	}
	return p;
}

/* Tropospheric delay correction --------------------------------------------------------------- */
void tropCorr_getEl(EPOCHOBS *obs, const double *Xr, double *BLHr, const int *locs, const int num)
{
	int i, idx;
	double El, Az, TropCorr, Xs[3];

	xyz2blh(Xr, BLHr);
	for (i = 0; i < num; i++) {
		idx = locs[i];
		Xs[0] = obs->SatPVT[idx].SatPos[0], Xs[1] = obs->SatPVT[idx].SatPos[1], Xs[2] = obs->SatPVT[idx].SatPos[2];
		compSatElAz(Xr, Xs, &El, &Az);
		TropCorr = Hopfield(BLHr[2], El);
		obs->SatPVT[idx].Azimuth = Az;
		obs->SatPVT[idx].Elevation = El;
		obs->SatPVT[idx].TropCorr = TropCorr;
		if (El < THRES_El*Rad) { obs->SatPVT[idx].Valid = false; }
		//when satPVT: 'track' 'halfc' ?
	}
}

/* Mark the outlier of obs.SatPVT[...] --------------------------------------------------------- */
int markOutlier(EPOCHOBS *obs)
{
	int idx;
	MWGF lastEpoch[MAXCHANNUM];
	for (idx = 0; idx < MAXCHANNUM; idx++) {
		deepCopy(&(obs->ComObs[idx]), &(lastEpoch[idx]));
		obs->ComObs[idx].reset();
	}
	for (idx = 0; idx < obs->SatNum; idx++) {
		if (initMWGF(obs, idx) < 0) {
			obs->SatPVT[idx].Valid = false;
			continue;
		}
	}
	detectOutlier(obs, lastEpoch);
	return 0;
}

/* Initialize the obs.ComObs[...] -------------------------------------------------------------- */
int initMWGF(EPOCHOBS *obs, const int idx)
{
	double f1, f2, p1, p2, l1, l2;
	obs->ComObs[idx].Prn = obs->SatObs[idx].Prn;
	obs->ComObs[idx].Sys = obs->SatObs[idx].System;
	obs->ComObs[idx].n = 1;
	if (obs->ComObs[idx].Sys == GPS) {
		f1 = FG1_GPS; f2 = FG2_GPS;
	}
	else if (obs->ComObs[idx].Sys == BDS) {
		f1 = FG1_BDS; f2 = FG3_BDS;
	}
	else return -1;
	l1 = obs->SatObs[idx].L[0]; l2 = obs->SatObs[idx].L[1];
	p1 = obs->SatObs[idx].P[0]; p2 = obs->SatObs[idx].P[1];
	obs->SatObs[idx].validP[0] = !(fabs(p1) < EPSILON);
	obs->SatObs[idx].validP[1] = !(fabs(p2) < EPSILON);
	obs->SatObs[idx].validL[0] = !(fabs(l1) < EPSILON);
	obs->SatObs[idx].validL[0] = !(fabs(l2) < EPSILON);
	if (fabs(p1) < EPSILON || fabs(p2) < EPSILON || fabs(l1) < EPSILON || fabs(l2) < EPSILON) { return -1; }

	//double widelen = CLIGHT / (f1 - f2);
	obs->ComObs[idx].PIF = (SQR(f1)*p1 - SQR(f2)*p2) / (SQR(f1) - SQR(f2));
	obs->ComObs[idx].LIF = (SQR(f1)*l1 - SQR(f2)*l2) / (SQR(f1) - SQR(f2));
	obs->ComObs[idx].LMW = ((f1*l1 - f2 * l2) / (f1 - f2) - (f1*p1 + f2 * p2) / (f1 + f2));
	obs->ComObs[idx].LGF = l1 - l2;
	obs->ComObs[idx].PGF = p1 - p2;

	obs->SatPVT[idx].PIF = obs->ComObs[idx].PIF;
	obs->SatPVT[idx].LIF = obs->ComObs[idx].LIF;
	return 0;
}

/* Detect the outlier of obs.SatPVT[...] using dMW, GF ----------------------------------------- */
void detectOutlier(EPOCHOBS *obs, MWGF *lastEpoch)
{
	double dLGF, dLMW; bool isValid; int i;
	for (i = 0; i < obs->SatNum; i++)
	{
		//SatPVT[i]: Unsuccessful solution, ephemeris unavailable, detection skips
		if (!obs->SatPVT[i].Sys) {
			//printf("%d %f\nThere is no satellite PVT result for the current epoch, and the ephemeris is unavailable,sys=%d,prn=%d\n",
			//		obs->Time.Week, obs->Time.SecOfWeek, obs->SatObs[i].System, obs->SatObs[i].Prn);
			continue;
		}
		if (!obs->SatPVT[i].Valid) {//Incomplete dual-frequency observation of the current epoch, detection skips
			//printf("%d %f\nThe dual-frequency observations of the current epoch are incomplete,sys=%d,prn=%d\n",
			//		obs->Time.Week, obs->Time.SecOfWeek, obs->SatPVT[i].Sys, obs->SatPVT[i].prn);
			continue;
		}
		bool hasfound = false;
		for (int j = 0; j < MAXCHANNUM; j++)
		{
			//the current epoch: Complete dual-frequency observations of SatObs[i], and has found the corresponding satellite 
			if ((obs->ComObs[i].Prn == lastEpoch[j].Prn) && (obs->ComObs[i].Sys == lastEpoch[j].Sys))
			{
				if (fabs(lastEpoch[j].LMW) < EPSILON)
				{	//Incomplete dual-frequency observations from the previous epoch, complete in the current epoch.
					//Poor observing conditions are believed to render them unusable.
					obs->ComObs[i].n = 1;
					obs->SatPVT[i].Valid = false;
					//printf("%d %f\nThe dual-frequency observation values in the previous epoch are incomplete and the current epoch is complete,sys=%d,prn=%d\n", 
					//		obs->Time.Week, obs->Time.SecOfWeek, obs->SatPVT[i].Sys, obs->SatPVT[i].prn);
					//system("pause");
					hasfound = true; break;//break
				}
				//Dual-frequency observations from both the previous and current epochs are complete
				//Cycle slips detection is conducted
				dLGF = obs->ComObs[i].LGF - lastEpoch[j].LGF;
				dLMW = obs->ComObs[i].LMW - lastEpoch[j].LMW;
				isValid = (dLGF <= THRES_dLGF) && (dLMW <= THRES_dLMW);
				obs->SatPVT[i].Valid = isValid;
				if (isValid) {//No Cycle slip
					obs->ComObs[i].n = lastEpoch[j].n + 1;
					obs->ComObs[i].LMW = (obs->ComObs[i].LMW + lastEpoch[j].LMW*lastEpoch[j].n) / (lastEpoch[j].n + 1);
				}
				else {//Cycle slips are detected
					obs->ComObs[i].n = 1;
					//printf("%d %f\nThe satellite has a cycle slip in the current epoch,sys=%d,prn=%d\n",
					//		obs->Time.Week, obs->Time.SecOfWeek, obs->SatPVT[i].Sys, obs->SatPVT[i].prn);
					//system("pause");
				}
				hasfound = true; break;
			}
		}
		//Complete dual-frequency observations in the current epoch
		//However, the corresponding satellite from the previous epoch could not be identified
		//Initial tracking does not confirm the absence of cycle slips; marked as unusable
		if (!hasfound) {
			obs->ComObs[i].n = 1;
			//printf("%d %f\nThe satellite was observed for the first time in the current epoch,sys=%d,prn=%d\n", 
			//		obs->Time.Week, obs->Time.SecOfWeek, obs->SatPVT[i].Sys, obs->SatPVT[i].prn);
			obs->SatPVT[i].Valid = false;
			//system("pause");
		}
	}
}

/* deepCopy 'MWGF' type ------------------------------------------------------------------------ */
void deepCopy(const MWGF *source, MWGF *dst)
{
	dst->Sys = source->Sys;
	dst->Prn = source->Prn;
	dst->n = source->n;
	dst->PIF = source->PIF;
	dst->LIF = source->LIF;
	dst->PGF = source->PGF;
	dst->LGF = source->LGF;
	dst->LMW = source->LMW;
}

/* make the obsdata/satPVTdata in order/accordance --------------------------------------------- */
int inOrder(EPOCHOBS *obs)
{
	SATOBS	  *ordsatobs = new SATOBS[MAXCHANNUM];
	SATRES *ordsatpvt = new SATRES[MAXCHANNUM];
	int idx = 0, i;
	for (i = 0; i < MAXCHANNUM; i++) {
		if (obs->SatObs[i].System) {
			memcpy(&ordsatobs[idx], &(obs->SatObs[i]), sizeof(SATOBS));
			memcpy(&ordsatpvt[idx], &(obs->SatPVT[i]), sizeof(SATRES));
			idx++;
		}
	}//shallow copy
	memcpy(obs->SatObs, ordsatobs, sizeof(SATOBS)*MAXCHANNUM);
	memcpy(obs->SatPVT, ordsatpvt, sizeof(SATRES)*MAXCHANNUM);
	delete[] ordsatobs;	delete[] ordsatpvt;
	return 1;
}

/* broadcast ephemeris to satellite position and clock bias -------------------------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,bds)
* args   : GPSEPHREC *ephgps     I     GPS broadcast ephemeris
*          BDSEPHREC *ephbds     I     Beidou broadcast ephemeris
*          EPOCHOBS  *obs        I/O   obs data of current epoch
* return : none
* notes  : satellite clock includes relativity correction, and code bias (tgd or bgd)
*----------------------------------------------------------------------------------------------- */
int compSatPosVel(GPSEPHREC *ephgps, BDSEPHREC *ephbds, EPOCHOBS *obs, const bool isSPP)
{
	int prn_1;
	double dt = 0.0;
	for (prn_1 = 0; prn_1 < MAXGPSNUM; prn_1++) {
		if (gpsPosVel(ephgps, obs, prn_1, isSPP, dt) > -1) {
			//printf("GPS sat: prn=%d has computed PVT successfully\n", j + 1);
			continue;
		}
	}
	for (prn_1 = 0; prn_1 < MAXBDSNUM; prn_1++) {
		if (bdsPosVel(ephbds, obs, prn_1, isSPP, dt) > -1) {
			//printf("BDS sat: prn=%d has computed PVT successfully\n", k + 1);
			continue;
		}
	}
	return 1;
}

/* compute gps satellite pos and vel, prn_1=prn-1------------------------------------------------------------- */
int gpsPosVel(GPSEPHREC *ephgps, EPOCHOBS *obs, const int prn_1, const bool isSPP, const double dt)
{
	double tk, n0, M, E, Ek, sinE, cosE, v, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, sini, cosi, mu, omge, ClkParams[3];
	int n, idx;

	idx = findSatidx(obs, prn_1 + 1, GPS);
	mu = GM_WGS84; omge = Omega_WGS;
	tk = TimeDiff(&(obs->Time), &(ephgps[prn_1].TOE));
	n0 = sqrt(mu) / pow(ephgps[prn_1].SqrtA, 3);
	ClkParams[0] = ephgps[prn_1].ClkBias;
	ClkParams[1] = ephgps[prn_1].ClkDrift;
	ClkParams[2] = ephgps[prn_1].ClkDriftRate;

	if (ephgps[prn_1].Sys && (!ephgps[prn_1].SVHealth) && (!(fabs(tk) > THRES_tk)) && (idx != -1))
	{	//ephemeris expired/didn't exist, or sat unhealty
		if (!isSPP) {
			obs->SatPVT[idx].SatClkOft = svClkCorr(TimeDiff(&(obs->Time), &(ephgps[prn_1].TOC)), ClkParams);
			obs->SatPVT[idx].SatClkSft = svClkdCorr(TimeDiff(&(obs->Time), &(ephgps[prn_1].TOC)), ClkParams);
			tk -= obs->SatPVT[idx].SatClkOft;
		}
		else {
			tk -= dt;
		}
		M = ephgps[prn_1].M0 + (n0 + ephgps[prn_1].DeltaN)*tk;
		for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER&&n < MAX_ITER_KEPLER; n++) {
			Ek = E; E -= (E - ephgps[prn_1].e*sin(E) - M) / (1.0 - ephgps[prn_1].e*cos(E));
		}
		if (n >= MAX_ITER_KEPLER) {
			printf("eph2pos: kepler iteration overflow GPS_prn=%d\n", ephgps[prn_1].PRN);
			return -1;
		}
		sinE = sin(E); cosE = cos(E);

		v = atan2(sqrt(1.0 - SQR(ephgps[prn_1].e))*sinE, cosE - ephgps[prn_1].e);
		u = v + ephgps[prn_1].omega;
		r = SQR(ephgps[prn_1].SqrtA)*(1.0 - ephgps[prn_1].e*cosE);
		i = ephgps[prn_1].i0 + ephgps[prn_1].iDot*tk;
		sin2u = sin(2.0*u); cos2u = cos(2.0*u);
		u += ephgps[prn_1].Cus*sin2u + ephgps[prn_1].Cuc*cos2u;
		r += ephgps[prn_1].Crs*sin2u + ephgps[prn_1].Crc*cos2u;
		i += ephgps[prn_1].Cis*sin2u + ephgps[prn_1].Cic*cos2u;
		x = r * cos(u); y = r * sin(u); cosi = cos(i); sini = sin(i);
		O = ephgps[prn_1].OMEGA + (ephgps[prn_1].OMEGADot - omge)*tk - omge * ephgps[prn_1].TOE.SecOfWeek;
		sinO = sin(O); cosO = cos(O);

		/* compute sat pos */
		double xk, yk, zk;
		xk = x * cosO - y * cosi*sinO;
		yk = x * sinO + y * cosi*cosO;
		zk = y * sini;

		/* compute sat vel */
		double E_dot, u_dot, f_dot, r_dot, I_dot, O_dot, xk_dot, yk_dot;
		E_dot = (n0 + ephgps[prn_1].DeltaN) / (1 - ephgps[prn_1].e*cosE);
		f_dot = E_dot * sqrt(1 - SQR(ephgps[prn_1].e)) / (1 - ephgps[prn_1].e*cosE);
		u_dot = f_dot + 2 * f_dot*(ephgps[prn_1].Cus*cos2u - ephgps[prn_1].Cuc*sin2u);
		r_dot = SQR(ephgps[prn_1].SqrtA)*ephgps[prn_1].e*sinE*E_dot + 2 * f_dot*(ephgps[prn_1].Crs*cos2u - ephgps[prn_1].Crc*sin2u);
		I_dot = ephgps[prn_1].iDot + 2 * f_dot*(ephgps[prn_1].Cis*cos2u - ephgps[prn_1].Cic*sin2u);
		O_dot = (ephgps[prn_1].OMEGADot - omge);
		xk_dot = r_dot * cos(u) - r * u_dot*sin(u);
		yk_dot = r_dot * sin(u) + r * u_dot*cos(u);

		double dXk, dYk, dZk;
		dXk = -yk * O_dot - (yk_dot*cosi - zk * I_dot)*sin(O) + xk_dot * cos(O);
		dYk = xk * O_dot + (yk_dot*cosi - zk * I_dot)*cos(O) + xk_dot * sin(O);
		dZk = yk_dot * sini + y * I_dot*cosi;

		obs->SatPVT[idx].Sys = GPS;
		obs->SatPVT[idx].prn = ephgps[prn_1].PRN;
		obs->SatPVT[idx].Tgd1 = ephgps[prn_1].TGD;
		obs->SatPVT[idx].SatPos[0] = xk;
		obs->SatPVT[idx].SatPos[1] = yk;
		obs->SatPVT[idx].SatPos[2] = zk;
		obs->SatPVT[idx].SatVel[0] = dXk;
		obs->SatPVT[idx].SatVel[1] = dYk;
		obs->SatPVT[idx].SatVel[2] = dZk;
		if (!isSPP) { obs->SatPVT[idx].Valid = true; }

	}
	else {
		return -1;
	}
	return 0;
}

/* compute bds satellite pos and vel, prn_1=prn-1------------------------------------------------------------- */
int bdsPosVel(BDSEPHREC *ephbds, EPOCHOBS *obs, const int prn_1, const bool isSPP, const double dt)
{
	double tk, n0, M, E, Ek, sinE, cosE, v, u, r, i, O, sin2u, cos2u, x, y, sinO, cosO, sini, cosi, mu, omge;
	double xg, yg, zg, sino, coso, ClkParams[3];
	int n, idx;

	idx = findSatidx(obs, prn_1 + 1, BDS);
	mu = GM_CGCS; omge = Omega_CGCS;
	tk = TimeDiff(&(obs->Time), &(ephbds[prn_1].TOE));
	n0 = sqrt(mu) / pow(ephbds[prn_1].SqrtA, 3);
	GPSTIME toeBDST(ephbds[prn_1].TOE.Week - GPST_BDT_WEEKS, ephbds[prn_1].TOE.SecOfWeek - GPST_BDT);
	ClkParams[0] = ephbds[prn_1].ClkBias;
	ClkParams[1] = ephbds[prn_1].ClkDrift;
	ClkParams[2] = ephbds[prn_1].ClkDriftRate;

	if (ephbds[prn_1].Sys && (!ephbds[prn_1].SVHealth) && (!(fabs(tk) > THRES_tk)) && (idx != -1))
	{	//ephemeris expired/didn't exist, or sat unhealty
		if (!isSPP) {
			obs->SatPVT[idx].SatClkOft = svClkCorr(TimeDiff(&(obs->Time), &(ephbds[prn_1].TOC)), ClkParams);
			obs->SatPVT[idx].SatClkSft = svClkdCorr(TimeDiff(&(obs->Time), &(ephbds[prn_1].TOC)), ClkParams);
			tk -= obs->SatPVT[idx].SatClkOft;
		}
		else {
			tk -= dt;
		}
		M = ephbds[prn_1].M0 + (n0 + ephbds[prn_1].DeltaN)*tk;
		for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER&&n < MAX_ITER_KEPLER; n++) {
			Ek = E; E -= (E - ephbds[prn_1].e*sin(E) - M) / (1.0 - ephbds[prn_1].e*cos(E));
		}
		if (n >= MAX_ITER_KEPLER) {
			printf("eph2pos: kepler iteration overflow GPS_prn=%d\n", ephbds[prn_1].PRN);
			return -1;  //iteration failed
		}
		sinE = sin(E); cosE = cos(E);

		v = atan2(sqrt(1.0 - SQR(ephbds[prn_1].e))*sinE, cosE - ephbds[prn_1].e);
		u = v + ephbds[prn_1].omega;
		r = SQR(ephbds[prn_1].SqrtA)*(1.0 - ephbds[prn_1].e*cosE);
		i = ephbds[prn_1].i0 + ephbds[prn_1].iDot*tk;
		sin2u = sin(2.0*u); cos2u = cos(2.0*u);
		u += ephbds[prn_1].Cus*sin2u + ephbds[prn_1].Cuc*cos2u;
		r += ephbds[prn_1].Crs*sin2u + ephbds[prn_1].Crc*cos2u;
		i += ephbds[prn_1].Cis*sin2u + ephbds[prn_1].Cic*cos2u;
		x = r * cos(u); y = r * sin(u); cosi = cos(i); sini = sin(i);

		/* compute sat pos */
		double xk, yk, zk, dXk, dYk, dZk;
		if (ephbds[prn_1].PRN <= 5 || ephbds[prn_1].PRN >= 59) {/* beidou geo satellite */
			O = ephbds[prn_1].OMEGA + ephbds[prn_1].OMEGADot*tk - omge * toeBDST.SecOfWeek;
			sinO = sin(O); cosO = cos(O);
			xg = x * cosO - y * cosi*sinO;
			yg = x * sinO + y * cosi*cosO;
			zg = y * sini;
			sino = sin(omge*tk); coso = cos(omge*tk);
			xk = xg * coso + yg * sino*COS_5 + zg * sino*SIN_5;
			yk = -xg * sino + yg * coso*COS_5 + zg * coso*SIN_5;
			zk = -yg * SIN_5 + zg * COS_5;

			/* compute sat vel */
			double E_dot, u_dot, f_dot, r_dot, I_dot, O_dot, xk_dot, yk_dot;
			E_dot = (n0 + ephbds[prn_1].DeltaN) / (1 - ephbds[prn_1].e*cosE);
			f_dot = E_dot * sqrt(1 - SQR(ephbds[prn_1].e)) / (1 - ephbds[prn_1].e*cosE);
			u_dot = f_dot + 2 * f_dot*(ephbds[prn_1].Cus*cos2u - ephbds[prn_1].Cuc*sin2u);
			r_dot = SQR(ephbds[prn_1].SqrtA)*ephbds[prn_1].e*sinE*E_dot + 2 * f_dot*(ephbds[prn_1].Crs*cos2u - ephbds[prn_1].Crc*sin2u);
			I_dot = ephbds[prn_1].iDot + 2 * f_dot*(ephbds[prn_1].Cis*cos2u - ephbds[prn_1].Cic*sin2u);
			O_dot = (ephbds[prn_1].OMEGADot);
			xk_dot = r_dot * cos(u) - r * u_dot*sin(u);
			yk_dot = r_dot * sin(u) + r * u_dot*cos(u);

			double dXg, dYg, dZg, dXk1, dYk1, dZk1, dXk2, dYk2, dZk2;
			dXg = -yg * O_dot - (yk_dot*cosi - zg * I_dot)*sin(O) + xk_dot * cos(O);
			dYg = xg * O_dot + (yk_dot*cosi - zg * I_dot)*cos(O) + xk_dot * sin(O);
			dZg = yk_dot * sini + y * I_dot*cosi;
			//part1
			dXk1 = dXg * coso + dYg * sino*COS_5 + dZg * sino*SIN_5;
			dYk1 = -dXg * sino + dYg * coso*COS_5 + dZg * coso*SIN_5;
			dZk1 = -dYg * SIN_5 + dZg * COS_5;
			//part2
			dXk2 = omge * (-sin(omge*tk)*xg + yg * COS_5*cos(omge*tk) + zg * SIN_5*cos(omge*tk));
			dYk2 = omge * (-cos(omge*tk)*xg - yg * COS_5*sin(omge*tk) - zg * SIN_5*sin(omge*tk));
			dZk2 = 0.0;
			//add
			dXk = dXk1 + dXk2;		dYk = dYk1 + dYk2;		dZk = dZk1 + dZk2;
		}
		else {
			O = ephbds[prn_1].OMEGA + (ephbds[prn_1].OMEGADot - omge)*tk - omge * toeBDST.SecOfWeek;
			sinO = sin(O); cosO = cos(O);
			xk = x * cosO - y * cosi*sinO;
			yk = x * sinO + y * cosi*cosO;
			zk = y * sini;

			/* compute sat vel */
			double E_dot, u_dot, f_dot, r_dot, I_dot, O_dot, xk_dot, yk_dot;
			E_dot = (n0 + ephbds[prn_1].DeltaN) / (1 - ephbds[prn_1].e*cosE);
			f_dot = E_dot * sqrt(1 - SQR(ephbds[prn_1].e)) / (1 - ephbds[prn_1].e*cosE);
			u_dot = f_dot + 2 * f_dot*(ephbds[prn_1].Cus*cos2u - ephbds[prn_1].Cuc*sin2u);
			r_dot = SQR(ephbds[prn_1].SqrtA)*ephbds[prn_1].e*sinE*E_dot + 2 * f_dot*(ephbds[prn_1].Crs*cos2u - ephbds[prn_1].Crc*sin2u);
			I_dot = ephbds[prn_1].iDot + 2 * f_dot*(ephbds[prn_1].Cis*cos2u - ephbds[prn_1].Cic*sin2u);
			O_dot = (ephbds[prn_1].OMEGADot - omge);
			xk_dot = r_dot * cos(u) - r * u_dot*sin(u);
			yk_dot = r_dot * sin(u) + r * u_dot*cos(u);

			dXk = -yk * O_dot - (yk_dot*cosi - zk * I_dot)*sin(O) + xk_dot * cos(O);
			dYk = xk * O_dot + (yk_dot*cosi - zk * I_dot)*cos(O) + xk_dot * sin(O);
			dZk = yk_dot * sini + y * I_dot*cosi;
		}
		obs->SatPVT[idx].Sys = BDS;
		obs->SatPVT[idx].prn = ephbds[prn_1].PRN;
		obs->SatPVT[idx].Tgd1 = ephbds[prn_1].TGD1;
		obs->SatPVT[idx].Tgd2 = ephbds[prn_1].TGD2;
		obs->SatPVT[idx].SatPos[0] = xk;
		obs->SatPVT[idx].SatPos[1] = yk;
		obs->SatPVT[idx].SatPos[2] = zk;
		obs->SatPVT[idx].SatVel[0] = dXk;
		obs->SatPVT[idx].SatVel[1] = dYk;
		obs->SatPVT[idx].SatVel[2] = dZk;
		if (!isSPP) { obs->SatPVT[idx].Valid = true; }

	}
	else {
		return -1;
	}
	return 0;
}

/* find the idx in obs.Satobs[idx] according to prn and sys ------------------------------------ */
int findSatidx(EPOCHOBS *obs, const short prn, const int sys) {
	int i;
	for (i = 0; i < MAXCHANNUM; i++) {
		if (obs->SatObs[i].Prn == prn && obs->SatObs[i].System == sys) { return i; }
	}
	return -1;
}

/* compute signal transmission time, Sat Clk/dClk, Sat POS/VEL with earth rotation correction -- */
int satClk_EarthRotCorr(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, const double *Xr)
{
	int prn, i, sys = UNKS;
	double e[3], dt, omge;

	for (i = 0; i < obs->SatNum; i++) {

		/* compute the Signal Transmission Time */
		if (!obs->SatPVT[i].Valid) { continue; }
		prn = obs->SatPVT[i].prn;
		sys = obs->SatPVT[i].Sys;
		dt = geodist(obs->SatPVT[i].SatPos, Xr, e) / CLIGHT;

		/* reCompute the Sat Clk/dClk */
		switch (sys) {
		case GPS:	gpsSatClk(obs, ephgps, prn - 1, dt); break;
		case BDS:	bdsSatClk(obs, ephbds, prn - 1, dt); break;
		default:
			printf("oem7 unknown system: sys=%d\n", sys);
		}

		/* reCompute the Sat POS/VEL without earth rotation correction */
		if (sys == GPS) {
			gpsPosVel(ephgps, obs, prn - 1, true, dt);
		}
		else if (sys == BDS) {
			bdsPosVel(ephbds, obs, prn - 1, true, dt);
		}

		/* earth rotation correction */
		omge = (sys == GPS) ? Omega_WGS : Omega_CGCS;
		double *pos = earthRotCorr(obs->SatPVT[i].SatPos, dt, omge);
		double *vel = earthRotCorr(obs->SatPVT[i].SatVel, dt, omge);
		obs->SatPVT[i].SatPos[0] = pos[0];
		obs->SatPVT[i].SatPos[1] = pos[1];
		obs->SatPVT[i].SatPos[2] = pos[2];
		obs->SatPVT[i].SatVel[0] = vel[0];
		obs->SatPVT[i].SatVel[1] = vel[1];
		obs->SatPVT[i].SatVel[2] = vel[2];
		delete[] pos; delete[] vel;
	}
	return 0;
}

/* compute GPS Sat Clk/dClk with signal transmission time -------------------------------------- */
int gpsSatClk(EPOCHOBS *obs, GPSEPHREC *ephgps, const int idx, const double dt)
{
	int n, ord_idx;
	double M, n0, tk, E, Ek, sinE, cosE, mu = GM_WGS84;
	double ClkParams[3], dts, ddts, E_dot;

	tk = TimeDiff(&(obs->Time), &(ephgps[idx].TOE)) - dt;
	n0 = sqrt(mu) / pow(ephgps[idx].SqrtA, 3);
	M = ephgps[idx].M0 + (n0 + ephgps[idx].DeltaN)*tk;
	for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER&&n < MAX_ITER_KEPLER; n++) {
		/* Newton's iteration method */
		Ek = E; E -= (E - ephgps[idx].e*sin(E) - M) / (1.0 - ephgps[idx].e*cos(E));
	}
	if (n >= MAX_ITER_KEPLER) {
		printf("eph2pos: kepler iteration overflow GPS_prn=%d\n", ephgps[idx].PRN);
		return -1; //iteration failed
	}
	sinE = sin(E); cosE = cos(E);

	ClkParams[0] = ephgps[idx].ClkBias;
	ClkParams[1] = ephgps[idx].ClkDrift;
	ClkParams[2] = ephgps[idx].ClkDriftRate;
	/* sat clock oft and sft! */
	dts = svClkCorr(TimeDiff(&(obs->Time), &(ephgps[idx].TOC)) - dt, ClkParams);
	ddts = svClkdCorr(TimeDiff(&(obs->Time), &(ephgps[idx].TOC)) - dt, ClkParams);
	/* relativity correction */
	dts -= 2.0*sqrt(mu*SQR(ephgps[idx].SqrtA))*ephgps[idx].e*sinE / SQR(CLIGHT);
	E_dot = (n0 + ephgps[idx].DeltaN) / (1 - ephgps[idx].e*cosE);
	ddts -= 2.0*sqrt(mu)*ephgps[idx].SqrtA*ephgps[idx].e*cosE*E_dot / SQR(CLIGHT);

	ord_idx = findSatidx(obs, idx + 1, GPS);
	obs->SatPVT[ord_idx].SatClkOft = dts;
	obs->SatPVT[ord_idx].SatClkSft = ddts;
	return 0;
}

/* compute BDS Sat Clk/dClk with signal transmission time -------------------------------------- */
int bdsSatClk(EPOCHOBS *obs, BDSEPHREC *ephbds, const int idx, const double dt)
{
	int n, ord_idx;
	double M, n0, tk, E, Ek, sinE, cosE, mu = GM_CGCS;
	double ClkParams[3], dts, ddts, E_dot;

	tk = TimeDiff(&(obs->Time), &(ephbds[idx].TOE)) - dt;
	n0 = sqrt(mu) / pow(ephbds[idx].SqrtA, 3);
	M = ephbds[idx].M0 + (n0 + ephbds[idx].DeltaN)*tk;
	for (n = 0, E = M, Ek = 0.0; fabs(E - Ek) > RTOL_KEPLER&&n < MAX_ITER_KEPLER; n++) {
		/* Newton's iteration method */
		Ek = E; E -= (E - ephbds[idx].e*sin(E) - M) / (1.0 - ephbds[idx].e*cos(E));
	}
	if (n >= MAX_ITER_KEPLER) {
		printf("eph2pos: kepler iteration overflow GPS_prn=%d\n", ephbds[idx].PRN);
		return -1; //iteration failed
	}
	sinE = sin(E); cosE = cos(E);

	ClkParams[0] = ephbds[idx].ClkBias;
	ClkParams[1] = ephbds[idx].ClkDrift;
	ClkParams[2] = ephbds[idx].ClkDriftRate;
	/* sat clock oft and sft! */
	dts = svClkCorr(TimeDiff(&(obs->Time), &(ephbds[idx].TOC)) - dt, ClkParams);
	ddts = svClkdCorr(TimeDiff(&(obs->Time), &(ephbds[idx].TOC)) - dt, ClkParams);
	/* relativity correction */
	dts -= 2.0*sqrt(mu*SQR(ephbds[idx].SqrtA))*ephbds[idx].e*sinE / SQR(CLIGHT);
	E_dot = (n0 + ephbds[idx].DeltaN) / (1 - ephbds[idx].e*cosE);
	ddts -= 2.0*sqrt(mu)*ephbds[idx].SqrtA*ephbds[idx].e*cosE*E_dot / SQR(CLIGHT);

	ord_idx = findSatidx(obs, idx + 1, BDS);
	obs->SatPVT[ord_idx].SatClkOft = dts;
	obs->SatPVT[ord_idx].SatClkSft = ddts;
	return 0;
}

/* spp-spv main function: LS -------------------------------------------------------------------- */
int SPP_SPV(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *pos)
{
	pos->Time.Week = obs->Time.Week;
	pos->Time.SecOfWeek = obs->Time.SecOfWeek;
#if !USE_FILTER
	//spp-LS
	if (spp_LS_Epoch(obs, ephgps, ephbds, pos) < 0) { return -1; }
	//spv-LS
	if (spv_LS_Epoch(obs, pos) < 0) { return -1; }

#else
	//spp-LS
	if (pos->ISFIRST) {
		//spp-LS
		if (spp_LS_Epoch(obs, ephgps, ephbds, pos) < 0) { return -1; }
		//spv-LS
		if (spv_LS_Epoch(obs, pos) < 0) { return -1; }
		eye(pos->Dxx, 6);
		scalar_matMul(InitVar, pos->Dxx, pos->Dxx, 36);
		pos->ISFIRST = false;
	}
	else {
		//spp-EKF
		if (spp_KF_Epoch(obs, ephgps, ephbds, pos) < 0) { return -1; }
		//spv-LS
		if (spv_LS_Epoch(obs, pos) < 0) { return -1; }
	}
#endif

	pos->IsSuccess = true;
	return 0;
}

/* spp main function --------------------------------------------------------------------------- */
int spp_LS_Epoch(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *pos)
{
	bool HasBeenSubed = false;
	int *locs = new int[obs->SatNum];
	int necessNum = 5, num, numgps, numbds, n = 0;
	double Pdop = 999.9, BLHr[3], sigmaPos = 999.9, dPos_norm = 1, Dxx[25], Qxx[25], v[40];
	//Pos[3] of the previous epoch 'pos' has not been cleared, and the first epoch is initialized to '0.0'
	double Xr[3] = { pos->Pos[0],pos->Pos[1],pos->Pos[2] };
	double rcvClkoft[2] = { pos->RcvClkOft[0], pos->RcvClkOft[1] };
	for (int i = 0; i < 25; i++) { Dxx[i] = 99.9; Qxx[i] = 999.9; }
	memset(v, 0, sizeof(double) * 40);

LOOP:
	while (n++<MAXITERS && dPos_norm>MAX_dPos_NORM) {
		double prePos[3], dPos[3];
		matdpcpy(prePos, Xr, 3);

		if (sppfindValidSat_locs(obs, locs, num, numgps, numbds) < necessNum + 1) { delete[] locs; return -1; }
		if ((!(numgps > 0 && numbds > 0)) && (!HasBeenSubed)) { necessNum--; HasBeenSubed = true; }
		if (spp_LS(obs, locs, Xr, rcvClkoft, Qxx, Dxx, necessNum, num, sigmaPos, Pdop, v) < 0) { delete[] locs; return -1; }
		tropCorr_getEl(obs, Xr, BLHr, locs, num);
		satClk_EarthRotCorr(obs, ephgps, ephbds, Xr);

		matSub(prePos, Xr, dPos, 3, 1);
		dPos_norm = norm(dPos, 3);
	}

	if ((fabs(maxElement(v, num, 1)) > ABS_V_MAX) || (fabs(minElement(v, num, 1)) > ABS_V_MAX))
	{
		for (int i = 0; i < num; i++)
		{
			if (fabs(v[i]) > ABS_V_MAX) {
				obs->SatPVT[locs[i]].Valid = false;
			}
		}
		n = 0;
		dPos_norm = 1;
		memset(v, 0, sizeof(double) * 40);
		goto LOOP;
	}

	matdpcpy(pos->Dxx, Dxx, necessNum*necessNum);
	matdpcpy(pos->Qxx, Qxx, necessNum*necessNum);
	pos->NecessaryNum = necessNum;
	pos->AllSatNum = num;
	pos->GPSSatNum = numgps;
	pos->BDSSatNum = numbds;
	pos->PDOP = Pdop;
	pos->SigmaPos = sigmaPos;
	pos->RcvClkOft[0] = rcvClkoft[0];
	pos->RcvClkOft[1] = rcvClkoft[1];
	pos->Pos[0] = Xr[0];
	pos->Pos[1] = Xr[1];
	pos->Pos[2] = Xr[2];

	delete[] locs;
	return 0;
}

/* spv main function --------------------------------------------------------------------------- */
int spv_LS_Epoch(EPOCHOBS *obs, POSRES *pos)
{
	int *locs = new int[obs->SatNum];
	int necessNums = 4, num;
	double sigmaVel = 999.9, rcvClksft = 0.0, Vel[3];
	double Xr[3] = { pos->Pos[0],pos->Pos[1],pos->Pos[2] };

	if (spvfindValidSat_locs(obs, locs, num) < necessNums + 1) { delete[] locs; return -1; }
	if (spv_LS(obs, locs, num, necessNums, Vel, Xr, sigmaVel, rcvClksft) < 0) { delete[] locs; return -1; }

	pos->RcvClkSft = rcvClksft;
	pos->SigmaVel = sigmaVel;
	pos->Vel[0] = Vel[0];
	pos->Vel[1] = Vel[1];
	pos->Vel[2] = Vel[2];

	delete[] locs;
	return 0;
}

int spp_KF_Epoch(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *pos)
{
	/* PART 1: backup Epoch k-1: Xk_1(rcvpos,rcvclock), Dxk_1 */
	double Xk_1[6] = { pos->Pos[0], pos->Pos[1], pos->Pos[2], pos->RcvClkOft[0], pos->RcvClkOft[1], pos->RcvClkSft };
	double Dxk_1[36], v[40];
	matdpcpy(Dxk_1, pos->Dxx, 36);
	memset(v, 0, sizeof(double) * 40);

	/* PART 2: LS for accurate satellite pos, vel, satclock, PDOP */
	bool HasBeenSubed = false;
	int *locs = new int[obs->SatNum];
	int necessNum = 5, num, numgps, numbds, n = 0;
	double Pdop = 999.9, BLHr[3], sigmaPos = 999.9, dPos_norm = 1, Dxx[36], Qxx[36];
	//Pos[3] of the previous epoch 'pos' has not been cleared, and the first epoch is initialized to '0.0'
	double Xr[3] = { pos->Pos[0],pos->Pos[1],pos->Pos[2] };
	double rcvClkoft[2] = { pos->RcvClkOft[0], pos->RcvClkOft[1] };
	for (int i = 0; i < 25; i++) { Dxx[i] = 99.9; Qxx[i] = 999.9; }

	LOOP:
	while (n++<MAXITERS && dPos_norm>MAX_dPos_NORM) {
		double prePos[3], dPos[3];
		matdpcpy(prePos, Xr, 3);

		if (sppfindValidSat_locs(obs, locs, num, numgps, numbds) < necessNum + 1) { delete[] locs; return -1; }
		if ((!(numgps > 0 && numbds > 0)) && (!HasBeenSubed)) { necessNum--; HasBeenSubed = true; }
		if (spp_LS(obs, locs, Xr, rcvClkoft, Qxx, Dxx, necessNum, num, sigmaPos, Pdop, v) < 0) { delete[] locs; return -1; }
		tropCorr_getEl(obs, Xr, BLHr, locs, num);
		satClk_EarthRotCorr(obs, ephgps, ephbds, Xr);

		matSub(prePos, Xr, dPos, 3, 1);
		dPos_norm = norm(dPos, 3);
	}

	if ((fabs(maxElement(v,num,1))>ABS_V_MAX)||(fabs(minElement(v, num, 1))>ABS_V_MAX))
	{
		for (int i = 0; i < num; i++)
		{
			if (fabs(v[i]) > ABS_V_MAX) {
				obs->SatPVT[locs[i]].Valid = false;
			}
		}
		n = 0;
		dPos_norm = 1;
		memset(v, 0, sizeof(double) * 40);
		goto LOOP;
	}

	/* PART 3: EKF */
	double Xkk_1[6], Dxkk_1[36], Xk[6], Dxk[36];
	oneStep_Prediction(Xk_1, Dxk_1, Xkk_1, Dxkk_1);
	measurement_Update(Dxkk_1, Xkk_1, Xk, Dxk, obs, locs, num);

	/* PART 4: save results */
	pos->Pos[0] = Xk[0];
	pos->Pos[1] = Xk[1];
	pos->Pos[2] = Xk[2];
	pos->RcvClkOft[0] = Xk[3];
	pos->RcvClkOft[1] = Xk[4];
	pos->RcvClkSft = Xk[5];
	pos->AllSatNum = num;
	pos->BDSSatNum = numbds;
	pos->GPSSatNum = numgps;
	matdpcpy(pos->Dxx, Dxk, 36);
	pos->IsSuccess = true;
	pos->NecessaryNum = 6;
	pos->SigmaPos = sqrt(Dxk[0] + Dxk[7] + Dxk[14]);
	pos->PDOP = Pdop;

	delete[] locs;

	return 0;
}

/* find valid locs for spp --------------------------------------------------------------------- */
int sppfindValidSat_locs(EPOCHOBS *obs, int *locs, int &num, int &numgps, int &numbds)
{
	int i;
	num = 0, numgps = 0, numbds = 0;
	//mark the valid idx
	for (i = 0; i < obs->SatNum; i++) {
		if (obs->SatPVT[i].Valid) {
			locs[num] = i;	num++;
			if(obs->SatPVT[i].Sys==GPS) { numgps++; }
			if(obs->SatPVT[i].Sys==BDS) { numbds++; }
		}
	}
	return num;
}

/* find valid locs for spv --------------------------------------------------------------------- */
int spvfindValidSat_locs(EPOCHOBS *obs, int *locs, int &num)
{
	int i;
	num = 0;
	//mark the valid idx
	for (i = 0; i < obs->SatNum; i++) {
		if ((obs->SatPVT[i].Valid) && (obs->SatObs[i].D[0]) != 0) {//only use L1C,B1I
			locs[num] = i;	num++;
		}
	}
	return num;
}

/* spp, Method: LS ----------------------------------------------------------------------------- */
int spp_LS(EPOCHOBS *obs, const int *locs, double *Xr, double *rcvClkoft, double *Qxx, double *Dxx,
	int necessNum, int num, double &sigmaPos, double &Pdop, double *v)
{
	int iters = 0; double maxCorr = 1.0;

	double *B = new double[num * necessNum];
	double *w = new double[num * 1];
	double *x = new double[necessNum * 1];
	double *BTPB = new double[necessNum*necessNum];
	double *inv_BTPB = new double[necessNum*necessNum];
	double *BTPw = new double[necessNum * 1];
	double *BT = new double[necessNum * num];
	double *xc = new double[necessNum - 3];
	double *Bx = new double[num * 1];

	while (iters++<MAXITERS&&maxCorr>LS_THRES)
	{
		double dXsr[3], d0, rClkOft0 = 0.0, sClkOft;
		int i, idx;
		// v=Bx-w: get B, w
		for (i = 0; i < num; i++) {
			idx = locs[i];
			double Xs[3] = { obs->SatPVT[idx].SatPos[0],obs->SatPVT[idx].SatPos[1],obs->SatPVT[idx].SatPos[2] };
			matSub(Xs, Xr, dXsr, 3, 1);
			d0 = norm(dXsr, 3);
			B[i * necessNum + 0] = (Xr[0] - Xs[0]) / d0;
			B[i * necessNum + 1] = (Xr[1] - Xs[1]) / d0;
			B[i * necessNum + 2] = (Xr[2] - Xs[2]) / d0;
			if (obs->SatPVT[idx].Sys == GPS) {
				B[i * necessNum + 3] = 1;
				if (necessNum == 5) { B[i * necessNum + 4] = 0; }
				rClkOft0 = rcvClkoft[0];
				sClkOft = CLIGHT * obs->SatPVT[idx].SatClkOft;
				w[i * 1] = getPIF(obs, idx, GPS) - (d0 + rClkOft0 - sClkOft + obs->SatPVT[idx].TropCorr);
			} else if (obs->SatPVT[idx].Sys == BDS) {
				B[i * necessNum + 3] = 0;
				if (necessNum == 5) { B[i * necessNum + 4] = 1; rClkOft0 = rcvClkoft[1]; }
				if (necessNum == 4) { B[i * necessNum + 3] = 1; rClkOft0 = rcvClkoft[0]; }//override it
				sClkOft = CLIGHT * obs->SatPVT[idx].SatClkOft;
				w[i * 1] = getPIF(obs, idx, BDS) - (d0 + rClkOft0 - sClkOft + obs->SatPVT[idx].TropCorr);
			}
		}
		//x
		matTran(B, num, necessNum, BT);
		matMul(BT, necessNum, num, B, num, necessNum, BTPB);
		matInv(BTPB, necessNum, inv_BTPB);
		matMul(BT, necessNum, num, w, num, 1, BTPw);
		matMul(inv_BTPB, necessNum, necessNum, BTPw, necessNum, 1, x);
		double x3[3] = { x[0],x[1],x[2] };
		if (necessNum == 5) { xc[0] = x[3], xc[1] = x[4]; }
		if (necessNum == 4) { xc[0] = x[3]; }
		//v	
		matMul(B, num, necessNum, x, necessNum, 1, Bx);
		matSub(Bx, w, v, num, 1);
		//renew maxCorr
		maxCorr = fabs(maxElement(x, necessNum, 1));
		//renew: Xr=Xr0+x, rcvClkbias=rcvClkbias0+x
		matAdd(Xr, x3, Xr, 3, 1);
		matAdd(rcvClkoft, xc, rcvClkoft, necessNum - 3, 1);
	}
	//sigma
	sigmaPos = sqrt(dot(v, v, num) / (num - necessNum));
	//Qxx=inv_BTPB
	Pdop = sqrt(SQR(inv_BTPB[0]) + SQR(inv_BTPB[1 * necessNum + 1]) + SQR(inv_BTPB[2 * necessNum + 2]));
	matdpcpy(Qxx, inv_BTPB, necessNum*necessNum);
	//Dxx=SQR(sigma)*Qxx
	scalar_matMul(SQR(sigmaPos), inv_BTPB, Dxx, necessNum*necessNum);

	//delete
	delete[] B;    delete[] w;        delete[] x;	 delete[] BT;
	delete[] BTPB; delete[] inv_BTPB; delete[] BTPw; delete[] Bx; delete[] xc;
	//LS failed
	if (iters > 9 || maxCorr > LS_THRES) { return -1; }
	return 0;
}

/* spv, Method: LS ----------------------------------------------------------------------------- */
int spv_LS(EPOCHOBS *obs, const int *locs, const int num, const int necessNum, double *vel, double *Xr,
	double &sigmaVel, double &rcvClksft)
{
	double *B = new double[num * necessNum], *w = new double[num * 1];
	double *x = new double[necessNum * 1], *v = new double[num * 1];
	double *BTPB = new double[necessNum*necessNum];
	double *inv_BTPB = new double[necessNum*necessNum];
	double *BTPw = new double[necessNum * 1];
	double *BT = new double[necessNum * num];
	double *Bx = new double[num * 1];
	double dXsr[3], d0, d0_dot, Dopp, wavelength;
	int i, idx;
	for (i = 0; i < num; i++) {
		idx = locs[i];
		wavelength = (obs->SatPVT[idx].Sys == GPS) ? WL1_GPS : WL1_BDS;
		double Xs[3] = { obs->SatPVT[idx].SatPos[0],obs->SatPVT[idx].SatPos[1],obs->SatPVT[idx].SatPos[2] };
		matSub(Xs, Xr, dXsr, 3, 1);
		Dopp = obs->SatObs[idx].D[0] * wavelength;
		d0 = norm(dXsr, 3);
		d0_dot = (obs->SatPVT[idx].SatPos[0] - Xr[0])*obs->SatPVT[idx].SatVel[0];
		d0_dot += (obs->SatPVT[idx].SatPos[1] - Xr[1])*obs->SatPVT[idx].SatVel[1];
		d0_dot += (obs->SatPVT[idx].SatPos[2] - Xr[2])*obs->SatPVT[idx].SatVel[2];
		d0_dot /= d0;
		B[i * necessNum + 0] = (-obs->SatPVT[idx].SatPos[0] + Xr[0]) / d0;
		B[i * necessNum + 1] = (-obs->SatPVT[idx].SatPos[1] + Xr[1]) / d0;
		B[i * necessNum + 2] = (-obs->SatPVT[idx].SatPos[2] + Xr[2]) / d0;
		B[i * necessNum + 3] = 1.0;
		w[i * 1] = Dopp - (d0_dot - CLIGHT * obs->SatPVT[idx].SatClkSft);
	}
	//x	
	matTran(B, num, necessNum, BT);
	matMul(BT, necessNum, num, B, num, necessNum, BTPB);
	matInv(BTPB, necessNum, inv_BTPB);
	matMul(BT, necessNum, num, w, num, 1, BTPw);
	matMul(inv_BTPB, necessNum, necessNum, BTPw, necessNum, 1, x);
	//v
	matMul(B, num, necessNum, x, necessNum, 1, Bx);
	matSub(Bx, w, v, num, 1);
	//sigma
	sigmaVel = sqrt(dot(v, v, num) / (num - necessNum));
	//renew Vel[], rcvClksft
	vel[0] = x[0];
	vel[1] = x[1];
	vel[2] = x[2];
	rcvClksft = x[3];
	//delete
	delete[] B;  delete[] w; delete[] v; delete[] x; delete[] Bx;
	delete[] BT; delete[] BTPB; delete[] inv_BTPB; delete[] BTPw;

	return 0;
}
