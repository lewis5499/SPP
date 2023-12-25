/*------------------------------------------------------------------------------
* strclib.h : SPP software structs' library
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references :  [1]"RTK_Structs.h"
*
* version : $Revision: 1.1 $ $Date: 2023/10/23 11:53:00 $
*
* history : 2023/10/10 1.0 new
*
*-----------------------------------------------------------------------------*/
#ifndef _STRUCT_H_
#define _STRUCT_H_

#pragma warning(disable:4996)

#include "consts.h"
#include <stdio.h>
#include <winsock.h>
#include <io.h>

/* Navigation satellite system definition */
enum GNSSSys { UNKS = 0, GPS, BDS, GLONASS, GALILEO, QZSS };

/* MJDT definition */
struct MJDTIME {             
	int Days;
	double FracDay;

	MJDTIME()
	{
		Days = 0;
		FracDay = 0.0;
	}
	MJDTIME(int d, double f)
	{
		Days = d;
		FracDay = f;
	}
};

/* Universal time definition */
struct COMMONTIME {   
	short Year;
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Minute;
	double         Second;

	COMMONTIME()
	{
		Year = 0;
		Month = 0;
		Day = 0;
		Hour = 0;
		Minute = 0;
		Second = 0.0;
	}
	COMMONTIME(short y,unsigned short mon, unsigned short d, unsigned short h, unsigned short min, double s)
	{
		Year = y;
		Month = mon;
		Day = d;
		Hour = h;
		Minute = min;
		Second = s;
	}
};

/* GPST definition */
struct GPSTIME {              
	unsigned short Week;
	double         SecOfWeek;

	GPSTIME()
	{
		Week = 0;
		SecOfWeek = 0.0;
	}
	GPSTIME(unsigned short w, double s)
	{
		Week = w;
		SecOfWeek = s;
	}
};

/* GPS broadcast ephemeris */
struct GPSEPHREC {			
	unsigned short PRN;
	GNSSSys     Sys;
	GPSTIME  	TOC, TOE;
	short		SVHealth;
	double		ClkBias, ClkDrift, ClkDriftRate;
	double		IODE1, IODC, IODE2;
	double      TGD;
	double		SqrtA, e, M0, OMEGA, i0, omega, OMEGADot, iDot, DeltaN;
	double		Crs, Cuc, Cus, Cic, Cis, Crc;
	double		SVAccuracy;

	GPSEPHREC() 
	{
		PRN = SVHealth = 0;
		Sys = UNKS;
		ClkBias = ClkDrift = ClkDriftRate = IODE1= IODE2 = IODC = TGD = 0.0;
		SqrtA = e = M0 = OMEGA = i0 = omega = OMEGADot = iDot = DeltaN = 0.0;
		Crs = Cuc = Cus = Cic = Cis = Crc = SVAccuracy = 0.0;
	}
};

/* Beidou broadcast ephemeris */
struct BDSEPHREC {
	unsigned short PRN;
	GNSSSys     Sys;
	GPSTIME  	TOC, TOE;
	short		SVHealth;
	double		ClkBias, ClkDrift, ClkDriftRate;
	double		IODE, IODC;
	double      TGD1, TGD2;
	double		SqrtA, e, M0, OMEGA, i0, omega, OMEGADot, iDot, DeltaN;
	double		Crs, Cuc, Cus, Cic, Cis, Crc;
	double		SVAccuracy;

	BDSEPHREC() 
	{
		PRN = SVHealth = 0;
		Sys = UNKS;
		ClkBias = ClkDrift = ClkDriftRate = IODE = IODC = TGD1 = TGD2 = 0.0;
		SqrtA = e = M0 = OMEGA = i0 = omega = OMEGADot = iDot = DeltaN = 0.0;
		Crs = Cuc = Cus = Cic = Cis = Crc = SVAccuracy = 0.0;
	}
};

/*  Observation data definition for each satellite  */
struct SATOBS {
	short    Prn;
	GNSSSys  System;
	double   P[2], L[2], D[2];   // m
	double   cn0[2], LockTime[2];
	unsigned char half[2];
	int		 parity[2], track[2];
	bool     Valid, validP[2], validL[2];

	SATOBS()
	{
		Prn = 0; 
		System = UNKS;
		parity[0] = parity[1] = 0; 
		track[0] = track[1] = 0;
		half[0] = half[1] = 0;
		LockTime[0] = LockTime[1] = 0.0;
		cn0[0] = cn0[1] = 0.0;
		P[0] = P[1] = L[0] = L[1] = D[0] = D[1] = 0.0;
		Valid = false;
		validP[0] = validP[1] = false;
		validL[0] = validL[1] = false;
	}
};

/* MW+GF observations */
struct MWGF {
	short Prn;//satellite number
	GNSSSys Sys;
	double LMW;
	double PGF, LGF;
	double PIF, LIF;
	int n;

	MWGF()
	{
		Prn = n = 0;
		Sys = UNKS;
		LMW = PGF = LGF = PIF = LIF = 0.0;
	}
	void reset() {
		Prn = n = 0;
		Sys = UNKS;
		LMW = PGF = LGF = PIF = LIF = 0.0;
	}
};

/* Intermediate calculations of each satellite's position, velocity, clock offset, etc. */
struct SATRES {
	short prn;
	GNSSSys Sys;
	double SatPos[3], SatVel[3];
	double SatClkOft, SatClkSft;
	double Elevation, Azimuth;
	double TropCorr;
	double Tgd1, Tgd2;
	double PIF, LIF;
	bool Valid;  //false=no ephemeris or ephemeris expired, or cycle slip
				 //true=calculation successful

	SATRES()
	{
		prn = 0; Sys = UNKS;
		SatPos[0] = SatPos[1] = SatPos[2] = 0.0;
		SatVel[0] = SatVel[1] = SatVel[2] = 0.0;
		Elevation = PI / 2.0;
		Azimuth = 0.0;
		SatClkOft = SatClkSft = 0.0;
		Tgd1 = Tgd2 = TropCorr = 0.0;
		PIF = LIF = 0.0;
		Valid = false;
	}
};

/* spp best pos */
struct BESTPOS {
	GPSTIME Time;
	double XYZ[3];
	double ENU[3];

	BESTPOS()
	{
		XYZ[0] = rcvXref;
		XYZ[1] = rcvYref;
		XYZ[2] = rcvZref;
		memset(ENU, 0, sizeof(double) * 3);
	}
};

/*  Observation data definition for each epoch  */
struct EPOCHOBS {
	GPSTIME Time;				//Current obs time: GPST.
	short   allobsNum;			//The total number of obs values from all obs satellites.
	short   obsNum;				//The total number of observations from observational satellites, if dual-frequency data
	                            //is complete, should be equal to 2 times the number of satellites (2*SatNum)
	short   SatNum;				//obs sats nums(GPS+BDS)
	SATOBS  SatObs[MAXCHANNUM]; // GPS, BDS in order
	SATRES  SatPVT[MAXCHANNUM]; // Calculation results such as satellite positions, the array index is the same as 'SatObs'
	MWGF    ComObs[MAXCHANNUM]; // Combined observations for the current epoch, the array index is the same as 'SatObs'
	BESTPOS BestPos;

	EPOCHOBS()
	{
		SatNum = obsNum = allobsNum = 0;
	}
	void reset() {
		SatNum = obsNum = allobsNum = 0;
		memset(&Time, 0, sizeof(GPSTIME));
		memset(SatObs, 0, sizeof(SATOBS)*MAXCHANNUM);
		memset(SatPVT, 0, sizeof(SATRES)*MAXCHANNUM);
	}
};

/* Definition of the positioning result for each epoch */
struct POSRES {
	GPSTIME Time;
	double Pos[3];
	double Vel[3];
	double RcvClkOft[2];              /* 0=GPS clock offset; 1=BDS clock offset */
	double RcvClkSft;
	double PDOP, SigmaPos, SigmaVel;  /* Accuracy index */
	double Qxx[6*6], Dxx[6*6];
	short  NecessaryNum;
	short  GPSSatNum, BDSSatNum;      /* Number of GPS/BDS satellites used for SPP */
	short  AllSatNum;                 /* Number of all satellites in the observation epoch */
	bool   IsSuccess;                 /* Record whether SPP is successful or not */
	bool   ISFIRST;

	POSRES()
	{
		Pos[0] = Pos[1] = Pos[2] = 0.0;
		Vel[0] = Vel[1] = Vel[2] = 0.0;
		for (int i = 0; i < 36; i++) { Dxx[i] = 99.9; Qxx[i] = 999.9; }
		RcvClkOft[0] = RcvClkOft[1] = RcvClkSft = 0.0;
		PDOP = SigmaPos = SigmaVel = 999.9;
		GPSSatNum = BDSSatNum = AllSatNum = 0;
		NecessaryNum = 5;
		IsSuccess = false;
		ISFIRST = true;
		if(USE_FILTER){ NecessaryNum = 6;}
	}
	void reset() {
		//Pos[0] = Pos[1] = Pos[2] = 0.0;
		//Vel[0] = Vel[1] = Vel[2] = 0.0;
		//RcvClkOft[0] = RcvClkOft[1] = RcvClkSft = 0.0;
		PDOP = SigmaPos = SigmaVel = 999.9;
		GPSSatNum = BDSSatNum = AllSatNum = 0;
		NecessaryNum = 5;
		IsSuccess = false;
		if(USE_FILTER){ NecessaryNum = 6;}
	}
};

/* Timer Callback Params */
struct TimerCallbackParams {
	SOCKET		NetGps;
	int8_t		*buff;
	int8_t		*Buff;
	EPOCHOBS	*obs;
	POSRES		*sppPos;
	GPSEPHREC	*ephgps;
	BDSEPHREC	*ephbds;
	FILE		*rawb;
	bool rawbIsOpen;
	int lenR;
	int lenD;
	int idx;
	const char* filewlog;
	const char* filewpos;
	const char* rawbinary;

	TimerCallbackParams() 
	{
		rawbIsOpen = false;
		lenR = 0, lenD = 0, idx = 0;
		obs 	= new EPOCHOBS;
		sppPos 	= new POSRES;
		buff 	= new int8_t[MAXRAWLEN];
		Buff 	= new int8_t[MAXRAWLEN * 3];
		ephgps 	= new GPSEPHREC[MAXGPSNUM];
		ephbds 	= new BDSEPHREC[MAXBDSNUM];

		if (access(CMAKELIST_PATH, 0) == -1) {
			filewlog = { LOGFILE0_REALTIME };
			filewpos = { POSFILE0_REALTIME };	
			rawbinary = { RAWFILE_BINARY0 };
		} else { 
			filewlog = { LOGFILE_RealTime };
			filewpos = { POSFILE_RealTime };	
			rawbinary = { RAWFILE_BINARY };
		}
	#if SAVE_REALTIME_RAW_BINARY
		if ((rawb = fopen(rawbinary, "wb")) == NULL) {
			printf("Failed to open binary file %s\n", rawbinary);
		} else {
			rawbIsOpen = true;
		}
	#endif
	}

	~TimerCallbackParams() 
	{
		if (rawbIsOpen) { fclose(rawb); }
		delete obs;	   		delete sppPos;
		delete[] buff; 		delete[] Buff;
		delete[] ephgps; 	delete[] ephbds;
	}
};
#endif