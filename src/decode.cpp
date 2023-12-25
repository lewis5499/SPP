#include "decode.h"
#include "dataIO.h"
#include "solution.h"
#include <cmath>

/* sync header ----------------------------------------------------------------------------------*/
int sync_oem7(uint8_t *buff, uint8_t data) {
	buff[0] = buff[1]; buff[1] = buff[2]; buff[2] = data;
	return buff[0] == OEM4SYNC1 && buff[1] == OEM4SYNC2 && buff[2] == OEM4SYNC3;
}

/* the whole decoding and algorithms for oem4/5/6/7 files----------------------------------------*/
int procbNovOem7f(FILE* fp, uint8_t *buff, const char *filewlog, const char *filewfinal) {
	EPOCHOBS  obs;
	POSRES	  sppPos;
	GPSEPHREC ephgps[MAXGPSNUM];
	BDSEPHREC ephbds[MAXBDSNUM];
	int data, idx;
	while (1) {
		memset(buff, 0, sizeof(uint8_t)*MAXRAWLEN);
		/* synchronize frame */
		for (int i = 0;; i++) {
			if ((data = fgetc(fp)) == EOF) { 
				printf("End of binary file reached.\n"); 
				fclose(fp);
				return 0; 
			}
			if (sync_oem7(buff, (uint8_t)data)) break;
		}
		idx = decodebNovOem7f(fp, buff, filewlog, filewfinal, &obs, ephgps, ephbds, &sppPos);
		if (!idx) {
			calculate_oem7(&obs, ephgps, ephbds, &sppPos);
			output_oem7f(filewlog, filewfinal, &sppPos, &obs);// Output SATpos/SPPpos to log/final(pos) file.
			reset_oem7(&obs, &sppPos);
		}
	}
}

/* decode one "AA4412" in buff ------------------------------------------------------------------
* it should be confirmed that buff[0]='AA',buff[1]='44',buff[2]='12'
* ---------------------------------------------------------------------------------------------- */
int decodebNovOem7f(FILE* fp, uint8_t *buff, const char *filewlog, const char *filewfinal, EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos) {
	double tow; 
	int msgtype; 
	uint16_t week, msglen, msgID;
	
	fread(buff + 3, 25, 1, fp);// has read 'header', buff[0]=AA, buff have 28 elements

	tow = U4(buff + 16)*0.001;
	week = U2(buff + 14);
	msglen = U2(buff + 8);
	msgID = U2(buff + 4);
	msgtype = (U1(buff + 6) >> 4) & 0x3; /* message type: 0=binary,1=ascii */

	fread(buff + OEM4HLEN, msglen + CRC32LEN, 1, fp);// read msg+CRC32, then 'buff' has OEM4HLEN+msglen+CRC32LEN elements

	return switchMsg(buff, obs, ephgps, ephbds, sppPos, msgID, msglen, msgtype, week, tow);
}

/* SOCKET decode one "AA4412" in buff----------------------------------------------------------- */
int decodebNovOem7s(int8_t *Buff, uint8_t *onebuff, EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos, int i, int &lenD) {
	double tow;
	int len_used, msgtype;	
	uint16_t week, msglen, msgID;

	bread(Buff + i + 1, onebuff + 3, 25);// has read 'header', onebuff[0]='AA', buff have 28 elements

	tow = U4(onebuff + 16)*0.001;
	week = U2(onebuff + 14);
	msglen = U2(onebuff + 8);
	msgID = U2(onebuff + 4);
	msgtype = (U1(onebuff + 6) >> 4) & 0x3; /* message type: 0=binary,1=ascii */

	len_used = OEM4HLEN + msglen + CRC32LEN + i - 2;
	if (len_used > lenD) { return -1; }
	
	bread(Buff + i + 26, onebuff + 28, msglen + CRC32LEN);// read msg+CRC32, then 'buff' has OEM4HLEN+msglen+CRC32LEN elements

	lenD -= len_used;
	memmove(Buff, Buff + len_used, lenD);

	return switchMsg(onebuff, obs, ephgps, ephbds, sppPos, msgID, msglen, msgtype, week, tow);
}

/* SOCKET decode one Message in buff------------------------------------------------------------ */
int switchMsg(uint8_t *buff, EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos, 
			  const uint16_t msgID, const uint16_t msglen, const int msgtype, const uint16_t week, const double tow) {
	/* check crc32 */
	if (rtk_crc32(buff, OEM4HLEN + msglen) != U4(buff + OEM4HLEN + msglen)) return 2;//'return 2'->'continue' in the outer loop

	// if crc is right,then:
	if (msgID == ID_RANGE && msgtype == 0) {
		int idx = decode_rangeb(buff, obs, msglen, week, tow);
		//printf("decode_rangeb, with [return %d]\n", idx);
		return 0;
	}else if (msgID == ID_GPSEPHEM && msgtype == 0) {
		decode_gpsephem(buff, ephgps, msglen);
		//printf("decode_gpsephem\n");
		return 1;
	}else if (msgID == ID_BDSEPHEM && msgtype == 0) {
		decode_bdsephem(buff, ephbds, msglen);
		//printf("decode_bdsephem\n");
		return 1;
	}else if (msgID == ID_BESTPOS && msgtype == 0) {
		decode_bestpos(buff, obs);
		//printf("decode_bestpos\n");
		return 1;
	}
	return 2;
}

/* SOCKET compute one OBSEPOCH in buff---------------------------------------------------------- */
int calculate_oem7(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos) {
	// perform satPos/Spp_pvt solution and save the result in a POSRES.
	compSatPosVel(ephgps, ephbds, obs, false);
	inOrder(obs);
	markOutlier(obs);
	SPP_SPV(obs, ephgps, ephbds, sppPos);	
	return 1;
}

/* check crc32 --------------------------------------------------------------------------------- */
uint32_t rtk_crc32(const uint8_t *buff, const int len){
	uint32_t crc = 0;
	int i, j;
	for (i = 0; i < len; i++) {
		crc ^= buff[i];
		for (j = 0; j < 8; j++) {
			if (crc & 1) crc = (crc >> 1) ^ POLYCRC32; else crc >>= 1;
		}
	}
	return crc;
}

/* signal type to obs code --------------------------------------------------------------------- */
int sig2code(int sys, int sigtype){
	if (sys == SYS_GPS) {
		switch (sigtype) {
		case  0: return CODE_L1C; /* L1C/A */
		case  5: return CODE_L2P; /* L2P    (OEM7) */
		case  9: return CODE_L2W; /* L2P(Y),semi-codeless */
		case 14: return CODE_L5Q; /* L5Q    (OEM6) */
		case 16: return CODE_L1L; /* L1C(P) (OEM7) */
		case 17: return CODE_L2S; /* L2C(M) (OEM7) */
		}
	}else if (sys == SYS_GLO) {
		switch (sigtype) {
		case  0: return CODE_L1C; /* L1C/A */
		case  1: return CODE_L2C; /* L2C/A (OEM6) */
		case  5: return CODE_L2P; /* L2P */
		case  6: return CODE_L3Q; /* L3Q   (OEM7) */
		}
	}else if (sys == SYS_GAL) {
		switch (sigtype) {
		case  2: return CODE_L1C; /* E1C  (OEM6) */
		case  6: return CODE_L6B; /* E6B  (OEM7) */
		case  7: return CODE_L6C; /* E6C  (OEM7) */
		case 12: return CODE_L5Q; /* E5aQ (OEM6) */
		case 17: return CODE_L7Q; /* E5bQ (OEM6) */
		case 20: return CODE_L8Q; /* AltBOCQ (OEM6) */
		}
	}else if (sys == SYS_QZS) {
		switch (sigtype) {
		case  0: return CODE_L1C; /* L1C/A */
		case 14: return CODE_L5Q; /* L5Q    (OEM6) */
		case 16: return CODE_L1L; /* L1C(P) (OEM7) */
		case 17: return CODE_L2S; /* L2C(M) (OEM7) */
		case 27: return CODE_L6L; /* L6P    (OEM7) */
		}
	}else if (sys == SYS_BDS) {
		switch (sigtype) {
		case  0: return CODE_L2I; /* B1I with D1 (OEM6) */
		case  1: return CODE_L7I; /* B2I with D1 (OEM6) */
		case  2: return CODE_L6I; /* B3I with D1 (OEM7) */
		case  4: return CODE_L2I; /* B1I with D2 (OEM6) */
		case  5: return CODE_L7I; /* B2I with D2 (OEM6) */
		case  6: return CODE_L6I; /* B3I with D2 (OEM7) */
		case  7: return CODE_L1P; /* B1C(P) (OEM7) */
		case  9: return CODE_L5P; /* B2a(P) (OEM7) */
		case 11: return CODE_L7D; /* B2b(I) (OEM7,F/W 7.08) */
		}
	}else if (sys == SYS_IRN) {
		switch (sigtype) {
		case  0: return CODE_L5A; /* L5 (OEM7) */
		}
	}else if (sys == SYS_SBS) {
		switch (sigtype) {
		case  0: return CODE_L1C; /* L1C/A */
		case  6: return CODE_L5I; /* L5I (OEM6) */
		}
	}
	return 0;
}

/* decode receiver tracking status --------------------------------------------------------------
* decode receiver tracking status
* args   : uint32_t stat I  tracking status field
*          int    *sys   O      system (SYS_???)
*          int    *code  O      signal code (CODE_L??)
*          int    *track O      tracking state
*                         (OEM4/5)
*                         0=L1 idle                   8=L2 idle
*                         1=L1 sky search             9=L2 p-code align
*                         2=L1 wide freq pull-in     10=L2 search
*                         3=L1 narrow freq pull-in   11=L2 pll
*                         4=L1 pll                   12=L2 steering
*                         5=L1 reacq
*                         6=L1 steering
*                         7=L1 fll
*                         (OEM6/7)
*                         0=idle                      7=freq-lock loop
*                         2=wide freq band pull-in    9=channel alignment
*                         3=narrow freq band pull-in 10=code search
*                         4=phase lock loop          11=aided phase lock loop
*          int    *plock O      phase-lock flag   (0=not locked, 1=locked)
*          int    *clock O      code-lock flag    (0=not locked, 1=locked)
*          int    *parity O     parity known flag (0=not known,  1=known)
*          int    *halfc O      phase measurement (0=half-cycle not added,
*                                                  1=added)
* return : freq-index (-1: when exclude[GPS:L1C/L2W BDS:BII/B3I] )
* notes  : null
*----------------------------------------------------------------------------------------------- */
int decode_track_stat(uint32_t stat, int *sys, int *code, int *track,
	int *plock, int *clock, int *parity, int *halfc, int *sigtype)
{
	int satsys, idx_code = -1;

	*code    = CODE_NONE;
	*track   = stat & 0x1F;
	*plock   = (stat >> 10) & 1;
	*parity  = (stat >> 11) & 1;
	*clock   = (stat >> 12) & 1;
	satsys   = (stat >> 16) & 7;	//0:GPS, 4:BDS
	*halfc   = (stat >> 28) & 1;
	*sigtype = (stat >> 21) & 0x1F; //sigtype = GPS(0:L1C,9:L2W), BDS(0,4:L2I,2,6:L6I)

	switch (satsys) {
	case 0: *sys = SYS_GPS; break;
	case 1: *sys = SYS_GLO; return -1; 
	case 2: *sys = SYS_SBS; return -1; 
	case 3: *sys = SYS_GAL; return -1;  /* OEM6 */
	case 4: *sys = SYS_BDS;	break;		/* OEM6 F/W 6.400 */
	case 5: *sys = SYS_QZS; return -1;  /* OEM6 */
	case 6: *sys = SYS_IRN; return -1;  /* OEM7 */
	default:
		printf("oem7 unknown system: sys=%d\n", satsys);
	}
	*code = sig2code(*sys, *sigtype);
	if (!(*code)) {
		printf("oem7 signal type error: sys=%d sigtype=%d\n", *sys, *sigtype);
		return -1;
	}
	//satsys = 0:GPS, 4:BDS  //sigtype = GPS(0:L1C,9:L2W), BDS(0,4:L2I,2,6:L6I)
	if ((satsys==0&&(*sigtype==0||*sigtype==9))||(satsys==4&&(*sigtype==0||*sigtype==4||*sigtype==2||*sigtype==6))) {
		idx_code = 0;
	}else{
		return -1;
	}
	return idx_code;
}

/* decode binary range message------------------------------------------------------------------ */
int decode_rangeb(uint8_t *buff, EPOCHOBS *obs, const uint16_t msglen, const uint16_t week, const double tow) {
	uint8_t *p = buff + OEM4HLEN;
	double psr, adr, dopp, snr, lockt;
	int i, nobs, prn, sys, code, idx, track, plock, clock, parity, halfc, sigtype, gps_bds_obscount=0;

	nobs = U4(p);
	if (msglen != 4 + nobs * 44) {
		printf("oem4 rangeb length error: len=%d nobs=%d\n", msglen, nobs);
		return -1;
	}
	for (i = 0, p += 4; i < nobs; i++, p += 44) {
		if ((idx = decode_track_stat(U4(p + 40), &sys, &code, &track, &plock, &clock,
			&parity, &halfc, &sigtype)) < 0) {
			continue;
		}
		gps_bds_obscount++;

		//sigtype = GPS(0:L1C,9:L2W), BDS(0,4:L2I,2,6:L6I)
		//P[loc] L[loc] D[loc] cn0[loc] LockTime[loc] half[loc]
		//The value of 'loc' can also be judged by 'code'
		//sys:SYS_GPS(0x01)(1) or SYS_BDS(0x20)(32)
		int loc = 0, sys_enum = 0, prnidx;
		double wavelength = 0.0;

		if (sys == SYS_GPS) {
			sys_enum = GPS;
			loc = sigtype == 0 ? 0 : 1;
			wavelength = (loc == 0) ? WL1_GPS : WL2_GPS;
		}else if (sys == SYS_BDS) {
			sys_enum = BDS;
			loc = (sigtype == 0 || sigtype == 4) ? 0 : 1;
			wavelength = (loc == 0) ? WL1_BDS : WL3_BDS;
		}
		prn = U2(p);
		psr = R8(p + 4);
		adr = -R8(p + 16)*wavelength;//adr = -R8(p + 16);
		dopp = -R4(p + 28);
		snr = R4(p + 32);//(uint16_t)(snr / SNR_UNIT + 0.5)
		lockt = R4(p + 36);

		if (!clock) { psr = 0.0; printf("psr-code lock: %d!\n", clock); }          /*  code unlock */
		if (!plock) { adr = dopp = 0.0;  printf("adr-phase lock: %d!\n", plock); } /* phase unlock */

		prnidx = sys == SYS_GPS ? prn : prn + MAXGPSNUM;

		((obs->SatObs) + prnidx - 1)->System = (GNSSSys)sys_enum;
		((obs->SatObs) + prnidx - 1)->Prn = prn;
		((obs->SatObs) + prnidx - 1)->P[loc] = psr;
		((obs->SatObs) + prnidx - 1)->L[loc] = adr;
		((obs->SatObs) + prnidx - 1)->D[loc] = dopp;
		((obs->SatObs) + prnidx - 1)->cn0[loc] = snr;
		((obs->SatObs) + prnidx - 1)->half[loc] = halfc;
		((obs->SatObs) + prnidx - 1)->LockTime[loc] = lockt;
		((obs->SatObs) + prnidx - 1)->parity[loc] = parity;
		((obs->SatObs) + prnidx - 1)->track[loc] = track;
		((obs->SatObs) + prnidx - 1)->Valid = true; //
	}
	int SatNumCount = 0;
	for (int i = 0; i < MAXCHANNUM; i++) {
		if (((obs->SatObs) + i)->System) {
			SatNumCount++;
		}
	}
	obs->Time.Week = week;
	obs->Time.SecOfWeek = tow;
	obs->allobsNum = nobs;
	obs->obsNum = gps_bds_obscount;
	obs->SatNum = SatNumCount;
	return 1;
}

/* decode binary gps ephemeris message---------------------------------------------------------- */
int decode_gpsephem(uint8_t *buff, GPSEPHREC *eph, const uint16_t msglen)
{
	uint8_t *p = buff + OEM4HLEN;
	if (msglen != GPSEPHEMLEN) {
		printf("oem7 gpsephemrisb length error: len=%d\n", msglen);
		return -1;
	}
	uint32_t prn = U4(p);		p += 4;
	double tow = R8(p);			p += 8;
	uint32_t health = U4(p);	p += 4;
	uint32_t IODE1 = U4(p);     p += 4;
	uint32_t IODE2 = U4(p);	    p += 4;
	uint32_t week = U4(p);		p += 4;
	uint32_t zweek = U4(p);		p += 4;
	double toe = R8(p);			p += 8;
	double sqrtA = sqrt(R8(p)); p += 8;
	double deln = R8(p);		p += 8;
	double M0 = R8(p);			p += 8;
	double e = R8(p);			p += 8;
	double omg = R8(p);			p += 8;
	double cuc = R8(p);			p += 8;
	double cus = R8(p);			p += 8;
	double crc = R8(p);			p += 8;
	double crs = R8(p);			p += 8;
	double cic = R8(p);			p += 8;
	double cis = R8(p);			p += 8;
	double i0 = R8(p);			p += 8;
	double idot = R8(p);		p += 8;
	double OMG0 = R8(p);		p += 8;
	double omgdot = R8(p);		p += 8;
	uint32_t iodc = U4(p);		p += 4;
	double toc = R8(p);			p += 8;
	double tgd = R8(p);			p += 8;
	double af0 = R8(p);			p += 8;
	double af1 = R8(p);			p += 8;
	double af2 = R8(p);			p += 8;
	/* omit AS */				p += 4;
	double N = R8(p);			p += 8;
	double ura = R8(p);			p += 8;

	eph[prn - 1].Cic = cic;
	eph[prn - 1].Cis = cis;
	eph[prn - 1].Crc = crc;
	eph[prn - 1].Crs = crs;
	eph[prn - 1].Cuc = cuc;
	eph[prn - 1].Cus = cus;
	eph[prn - 1].ClkBias = af0;
	eph[prn - 1].ClkDrift = af1;
	eph[prn - 1].ClkDriftRate = af2;
	eph[prn - 1].DeltaN = deln;
	eph[prn - 1].e = e;
	eph[prn - 1].i0 = i0;
	eph[prn - 1].iDot = idot;
	eph[prn - 1].IODC = iodc;
	eph[prn - 1].IODE1 = IODE1;
	eph[prn - 1].IODE2 = IODE2;
	eph[prn - 1].M0 = M0;
	eph[prn - 1].OMEGA = OMG0;
	eph[prn - 1].omega = omg;
	eph[prn - 1].OMEGADot = omgdot;
	eph[prn - 1].PRN = prn;
	eph[prn - 1].SqrtA = sqrtA;
	eph[prn - 1].SVAccuracy = ura;
	eph[prn - 1].SVHealth = health;
	eph[prn - 1].Sys = GPS;
	eph[prn - 1].TGD = tgd;
	eph[prn - 1].TOC.Week = week;
	eph[prn - 1].TOC.SecOfWeek = toc;
	eph[prn - 1].TOE.Week = week;
	eph[prn - 1].TOE.SecOfWeek = toe;

	return 1;
}

/* decode binary bds ephemeris message---------------------------------------------------------- */
int decode_bdsephem(uint8_t *buff, BDSEPHREC *eph, const uint16_t msglen)
{
	uint8_t *p = buff + OEM4HLEN;

	if (msglen != BDSEPHEMLEN) {
		printf("oem7 bdsephemrisb length error: len=%d\n", msglen);
		return -1;
	}
	int prn = U4(p);			 p += 4;
	uint32_t week = U4(p);		 p += 4;
	double ura = R8(p);			 p += 8;
	uint32_t health = U4(p)&1;   p += 4;  /* 0:good, 1:not good: broadcasting satellite */
	double tgd1 = R8(p);		 p += 8;  /* TGD1 for B1 (s) */
	double tgd2 = R8(p);		 p += 8;  /* TGD2 for B2 (s) */
	uint32_t AODC = U4(p);		 p += 4;  /* AODC */
	uint32_t toc = U4(p);		 p += 4;
	double a0 = R8(p);			 p += 8;
	double a1 = R8(p);			 p += 8;
	double a2 = R8(p);			 p += 8;
	uint32_t AODE = U4(p);		 p += 4;  /* AODE */
	uint32_t toe = U4(p);	  	 p += 4;
	double sqrtA = R8(p);		 p += 8;
	double e = R8(p);			 p += 8;
	double omg = R8(p);			 p += 8;
	double deln = R8(p);	     p += 8;
	double M0 = R8(p);		 	 p += 8;
	double OMG0 = R8(p);		 p += 8;
	double OMGd = R8(p);		 p += 8;
	double i0 = R8(p);			 p += 8;
	double IDOT = R8(p);		 p += 8;
	double cuc = R8(p);			 p += 8;
	double cus = R8(p);			 p += 8;
	double crc = R8(p);			 p += 8;
	double crs = R8(p);			 p += 8;
	double cic = R8(p);		     p += 8;
	double cis = R8(p);			 p += 8;

	// GPST-1356weeks-14sec=BDST
	eph[prn - 1].TOE.Week = week + GPST_BDT_WEEKS;
	eph[prn - 1].TOE.SecOfWeek = toe + GPST_BDT;
	eph[prn - 1].TOC.Week = week + GPST_BDT_WEEKS;
	eph[prn - 1].TOC.SecOfWeek = toc + GPST_BDT;
	eph[prn - 1].TGD1 = tgd1;
	eph[prn - 1].TGD2 = tgd2;
	eph[prn - 1].Sys = BDS;
	eph[prn - 1].SVHealth = health;
	eph[prn - 1].SVAccuracy = ura;
	eph[prn - 1].SqrtA = sqrtA;
	eph[prn - 1].PRN = prn;
	eph[prn - 1].OMEGADot = OMGd;
	eph[prn - 1].omega = omg;
	eph[prn - 1].OMEGA = OMG0;
	eph[prn - 1].M0 = M0;
	eph[prn - 1].IODE = AODE;
	eph[prn - 1].IODC = AODC;
	eph[prn - 1].iDot = IDOT;
	eph[prn - 1].i0 = i0;
	eph[prn - 1].e = e;
	eph[prn - 1].DeltaN = deln;
	eph[prn - 1].Cus = cus;
	eph[prn - 1].Cuc = cuc;
	eph[prn - 1].Crs = crs;
	eph[prn - 1].Crc = crc;
	eph[prn - 1].Cic = cic;
	eph[prn - 1].Cis = cis;
	eph[prn - 1].ClkBias = a0;
	eph[prn - 1].ClkDrift = a1;
	eph[prn - 1].ClkDriftRate = a2;

	return 1;
}

int decode_bestpos(uint8_t *buff, EPOCHOBS *obs)
{
	uint8_t *p = buff + OEM4HLEN + 8;
	double BLH[3];
	float undulation;

	BLH[0] = R8(p)*Rad;	p += 8;
	BLH[1] = R8(p)*Rad;	p += 8;
	BLH[2] = R8(p);		p += 8;
	undulation = R4(p);

	BLH[2] += undulation;

	blh2xyz(BLH, obs->BestPos.XYZ);
	/*printf("BESTPOS: %.4f\t%.4f\t%.4f\t\t%.4f\t%.4f\t%.4f\n", obs->BestPos.XYZ[0], obs->BestPos.XYZ[1], obs->BestPos.XYZ[2],
		BLH[0]*Deg, BLH[1]*Deg, BLH[2]);*/
	return 1;
}