/*------------------------------------------------------------------------------
* decode.h : SPP software decoding operations
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references :  [1]"RTK_Structs.h" [2] rtklib
*
* version : $Revision: 1.2 $ $Date: 2023/11/09 14:12:28 $
*
* history : 2023/11/07 1.2 new
*
*-----------------------------------------------------------------------------*/
#ifndef _DECODE_
#define _DECODE_

#include "strclib.h"
#include <stdio.h>
#include <cstring>

/* get fields (little-endian) ------------------------------------------------*/
#define U1(p) (*((uint8_t *)(p)))
#define I1(p) (*((int8_t  *)(p)))
static uint16_t U2(uint8_t *p) { uint16_t u; memcpy(&u, p, 2); return u; }
static uint32_t U4(uint8_t *p) { uint32_t u; memcpy(&u, p, 4); return u; }
static int32_t  I4(uint8_t *p) { int32_t  i; memcpy(&i, p, 4); return i; }
static float    R4(uint8_t *p) { float    r; memcpy(&r, p, 4); return r; }
static double   R8(uint8_t *p) { double   r; memcpy(&r, p, 8); return r; }

/* decode NovOem7 raw binary file/socket data --------------------------------*/
int sync_oem7(uint8_t *buff, uint8_t data);
int procbNovOem7f  (FILE* fp, uint8_t *buff, const char *filewlog, const char *filewfinal);
int decodebNovOem7f(FILE* fp, uint8_t *buff, const char *filewlog, const char *filewfinal, 
					EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos);
int decodebNovOem7s(int8_t *Buff, uint8_t *onebuff, EPOCHOBS *obs, GPSEPHREC *ephgps, 
					BDSEPHREC *ephbds, POSRES *sppPos, int i, int &lenD);
int switchMsg(uint8_t *buff, EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos, 
	const uint16_t msgID, const uint16_t msglen, const int msgtype, const uint16_t week, const double tow);

/* calculate NovOem7 raw binary file/socket data -----------------------------*/
int calculate_oem7(EPOCHOBS *obs, GPSEPHREC *ephgps, BDSEPHREC *ephbds, POSRES *sppPos);

/* check crc-32bit -----------------------------------------------------------*/
uint32_t rtk_crc32(const uint8_t *buff, const int len);

/* decode rangeb: based on rtklib --------------------------------------------*/
int decode_track_stat(uint32_t stat, int *sys, int *code, int *track, int *plock, int *clock,
						int *parity, int *halfc, int *sigtype);
int sig2code(int sys, int sigtype);
int decode_rangeb(uint8_t *buff, EPOCHOBS *obs, const uint16_t msglen, const uint16_t week, const double tow);

/* decode ephemeris ----------------------------------------------------------*/
int decode_gpsephem(uint8_t *buff, GPSEPHREC *eph, const uint16_t msglen);
int decode_bdsephem(uint8_t *buff, BDSEPHREC *eph, const uint16_t msglen);

/* decode bestpos ------------------------------------------------------------*/
int decode_bestpos(uint8_t *buff, EPOCHOBS *obs);

#endif