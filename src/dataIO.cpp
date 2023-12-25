#include "dataIO.h"
#include "consts.h"
#include "sockets.h"
#include "decode.h"
#include "solution.h"
#include <sys/stat.h>
#include <direct.h>
#include <io.h>
#include <string>

#pragma warning(disable:4996)

using namespace std;

/* start MENU ---------------------------------------------------------------------------- */
void Menu(int mode)
{
	switch (mode)
	{
	case MODE_FILE:	   sppFileThread();   break;
	case MODE_SOCKET:  sppSocketThread(); break;
	default:							  return;			   
	}
	system("pause");
}

/* start file program ---------------------------------------------------------------------*/
void sppFileThread()
{
	bool notInProj = (access(CMAKELIST_PATH, 0) == -1);
	string logpath = notInProj ? OUTPATH0 : OUTPATH;
	string filerb = notInProj ? DATAPATH0 : DATAPATH;
	string filewlog = notInProj ? LOGFILE_PostProc0 : LOGFILE_PostProc;
	string filewpos = notInProj ? POSFILE_PostProc0 : POSFILE_PostProc;

	if (access(logpath.c_str(), 0) == -1) {
		if (mkdir(logpath.c_str())) {
			printf("Unable to create directory.\n");
		}
	}
	printf("%s", INIT_INFO_FILE);
	printf("Start decoding %s...\n", filerb.c_str());

	if (input_oem7f(filerb.c_str(), filewlog.c_str(), filewpos.c_str()) < 0) {
		printf("fail to start decoding process: FILE ERROR!\n");
	}
	return;
}

/* init the files to write ----------------------------------------------------------------*/
int initNewFile(const char *filewlog, const char *filewfinal)
{
	FILE *fp_pos, *fp_log;
	fp_pos = fopen(filewfinal, "w");
	fp_log = fopen(filewlog, "w");
	if (fp_pos && fp_log) {
		fprintf(fp_pos, "%s\n", FINALPOSHEADER_PosFile);
	} else {
		printf("fail to create log/pos file!\n");
		return -1;
	}
	fclose(fp_pos); 
	fclose(fp_log);
	return 0;
}

/* open file and start decode -------------------------------------------------------------*/
int input_oem7f(const char *filer, const char *filewlog, const char *filewpos)
{
	FILE* fp;
	uint8_t buff[MAXRAWLEN];
	if ((fp = fopen(filer, "rb")) == NULL) {
		printf("Cannot open raw obsfile. \n");
		return -1;
	}
	if (initNewFile(filewlog, filewpos) < 0) { return -1; }
	printf("\nDecoding...\n\n");
	return procbNovOem7f(fp, buff, filewlog, filewpos);
}

/* input socket Buff data, start to decode ------------------------------------------------*/
int input_oem7(int8_t *Buff, int &lenD, EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos)
{
	uint8_t *onebuff = new uint8_t[MAXONERAWLEN];
	int data, i, idx;
	bool haveAA4412;
	while (1) {
		haveAA4412 = false;	idx = 2;
		/* synchronize frame */
		for (i=0;i<lenD; i++) {
			data = I1(Buff + i);
			if (sync_oem7(onebuff, (uint8_t)data)) { haveAA4412 = true; break; }//onebuff[0]~[2]:'AA4412'
		}																		//(Buff+i-2)begins with 'AA4412'
		if (!haveAA4412) { break; }
		idx=decodebNovOem7s(Buff, onebuff, obs, ephgps, ephbds, sppPos, i, lenD);
		memset(onebuff, 0, sizeof(uint8_t)*MAXONERAWLEN);

		if (idx < 0) { break; }//'Buff': insufficient length!

		if (!idx) { delete[] onebuff; return idx; } /* idx==0 means OBS is decoded */
	}
	delete[] onebuff;
	return 1;
}

/* print spp/spv and satellite's pvt results to file --------------------------------------*/
int output_oem7f(const char *filewlog, const char *filewpos, POSRES *sppPos, EPOCHOBS *obs)
{
	FILE *fp_log, *fp_pos;
	if ((fp_log = fopen(filewlog, "a+")) == NULL) {
		printf("Cannot open log result file. \n");
		system("pause");
		return -1;
	}
	if ((fp_pos = fopen(filewpos, "a+")) == NULL) {
		printf("Cannot open final 'pos' result file. \n");
		system("pause");
		return -1;
	}
	if (!sppPos->IsSuccess) { return -1; }
	//log
	PrintSPPlog(sppPos, fp_log);
	int idx = 0;
	while (1) {
		if (obs->SatPVT[idx++].Valid) {	PrintSatPVTlog(&(obs->SatPVT[idx-1]), fp_log);	}
		if (idx >= obs->SatNum) { break; }
	}
	//spp+spv
	Print_spp_spv(filewpos, fp_pos, obs, sppPos);
	fclose(fp_log); fclose(fp_pos);
	return 0;
}

/* output the SOCKET spp/spv results ------------------------------------------------------*/
int output_oem7(EPOCHOBS *obs, POSRES *sppPos)
{
	if (!sppPos->IsSuccess) { return -1; }
	double BLH[3], ENU[3];
	int i = 0; 
	string str = "";

	xyz2blh(sppPos->Pos, BLH);
	xyz2enu(sppPos->Pos, obs->BestPos.XYZ, ENU);
	while (i++ < MAXCHANNUM){
		if (obs->SatPVT[i - 1].Valid){
			if (obs->SatPVT[i - 1].Sys == GPS){
				str += "G";
				str += (obs->SatPVT[i - 1].prn < 10) ? "0" + to_string(obs->SatPVT[i - 1].prn) : to_string(obs->SatPVT[i - 1].prn);
			}else if(obs->SatPVT[i - 1].Sys == BDS){
				str += "C";
				str += (obs->SatPVT[i - 1].prn < 10) ? "0" + to_string(obs->SatPVT[i - 1].prn) : to_string(obs->SatPVT[i - 1].prn);
			}
		}
	}
	printf("%4d %.3f %13.3f %13.3f %13.3f %9.3f %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %7.3f %8.3f %8.3f %8.3f %3d %3d %3d\n",
			sppPos->Time.Week, sppPos->Time.SecOfWeek, sppPos->Pos[0], sppPos->Pos[1], sppPos->Pos[2],
			BLH[0] * Deg, BLH[1] * Deg, BLH[2], ENU[0], ENU[1], ENU[2], sppPos->Vel[0], sppPos->Vel[1], sppPos->Vel[2],
			sppPos->PDOP, sppPos->SigmaPos, sppPos->SigmaVel, sppPos->GPSSatNum, sppPos->BDSSatNum, sppPos->AllSatNum);
	return 0;
}

void reset_oem7(EPOCHOBS *obs, POSRES *sppPos)
{
	(*obs).reset();
	(*sppPos).reset();
}

/* output the FILE SatPVT results ---------------------------------------------------------*/
void PrintSatPVTlog(SATRES *obsPVT, FILE* fp)
{
	string PRN = "";
	if (obsPVT->Sys == GPS) { PRN += "G"; }
	if (obsPVT->Sys == BDS) { PRN += "C"; }
	if (obsPVT->prn < 10) { PRN += "0" + to_string(obsPVT->prn); }
	if (obsPVT->prn >= 10) { PRN += to_string(obsPVT->prn); }
	double X, Y, Z, Vx, Vy, Vz, Clk, Clkd, PIF, Trop, E;

	X = obsPVT->SatPos[0], Y = obsPVT->SatPos[1], Z = obsPVT->SatPos[2];
	Vx = obsPVT->SatVel[0], Vy = obsPVT->SatVel[1], Vz = obsPVT->SatVel[2];
	Clk = obsPVT->SatClkOft, Clkd = obsPVT->SatClkSft;
	PIF = obsPVT->PIF, Trop = obsPVT->TropCorr, E = (obsPVT->Elevation)*Deg;
	fprintf(fp, "%s X=%13.3f Y=%13.3f Z=%13.3f Clk=%13.6e Vx=%13.4f Vy=%13.4f Vz=%13.4f Clkd=%13.6e PIF=%13.4f Trop=%7.3f E=%7.3fdeg\n",
				PRN.c_str(), X, Y, Z, Clk, Vx, Vy, Vz, Clkd, PIF, Trop, E);
}

/* output the FILE SPP results ------------------------------------------------------------*/
void PrintSPPlog(POSRES *spp, FILE *fp)
{
	int week = spp->Time.Week;
	double toe = spp->Time.SecOfWeek;
	double gpsClk, bdsClk, PDOP, SigmaPos, BLHr[3];
	short GPSSats, BDSSats, Sats;

	xyz2blh(spp->Pos, BLHr);
	gpsClk = spp->RcvClkOft[0]; bdsClk = spp->RcvClkOft[1];
	PDOP = spp->PDOP; SigmaPos = spp->SigmaPos;
	GPSSats = spp->GPSSatNum; BDSSats = spp->BDSSatNum;
	Sats = spp->AllSatNum;
	fprintf(fp, "> %d %8.3f Sats:%2d GPSSats:%2d BDSSats:%2d \nSPP X:%13.3f Y:%13.3f Z:%13.3f B:%7.3f L:%7.3f H:%7.3f gpsClk:%7.3f bdsClk:%7.3f PDOP:%7.3f SigmaPos:%7.3f\n",
			week, toe, Sats, GPSSats, BDSSats, spp->Pos[0], spp->Pos[1], spp->Pos[2], BLHr[0]*Deg, BLHr[1]*Deg, BLHr[2], gpsClk, bdsClk, PDOP, SigmaPos);
}

/* output the FILE SPP/SPV pos results ----------------------------------------------------*/
void Print_spp_spv(const char *filewpos, FILE *fp, EPOCHOBS *obs, POSRES *sppPos)
{
	if (!sppPos->IsSuccess) { return; }
	int i = 0; 
	string str = "";
	double BLH[3], ENU[3];

	xyz2blh(sppPos->Pos, BLH);
	xyz2enu(sppPos->Pos, obs->BestPos.XYZ, ENU);
	while (i++<MAXCHANNUM){
		if (obs->SatPVT[i-1].Valid){
			if (obs->SatPVT[i-1].Sys==GPS){
				str += "G";
				str+=(obs->SatPVT[i-1].prn<10)?"0"+to_string(obs->SatPVT[i-1].prn):to_string(obs->SatPVT[i-1].prn);
			}else if (obs->SatPVT[i-1].Sys==BDS){
				str += "C";
				str+=(obs->SatPVT[i-1].prn<10)?"0"+to_string(obs->SatPVT[i-1].prn):to_string(obs->SatPVT[i-1].prn);
			}
		}
	}
	fprintf(fp, "%4d %.3f %14.4f %14.4f %14.4f %14.4f %14.4f %14.4f %7.3f %7.3f %7.3f %14.8f %14.8f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %3d %3d %3d %s\n",
			sppPos->Time.Week,sppPos->Time.SecOfWeek, sppPos->Pos[0], sppPos->Pos[1], sppPos->Pos[2], obs->BestPos.XYZ[0], obs->BestPos.XYZ[1], obs->BestPos.XYZ[2],
			ENU[0],ENU[1],ENU[2],BLH[0]*Deg,BLH[1]*Deg,BLH[2],sppPos->Vel[0], sppPos->Vel[1], sppPos->Vel[2], 
			sppPos->PDOP,sppPos->SigmaPos,sppPos->SigmaVel,sppPos->GPSSatNum,sppPos->BDSSatNum,sppPos->AllSatNum, str.c_str());
}
