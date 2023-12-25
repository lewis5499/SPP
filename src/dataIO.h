/*------------------------------------------------------------------------------
* dataIO.h : SPP software data stream operations
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references :  [1]"RTK_Structs.h" [2] rtklib
*
* version : $Revision: 1.1 $ $Date: 2023/11/3 13:53:00 $
*
* history : 2023/11/03 1.0 new
*
*-----------------------------------------------------------------------------*/
#ifndef _FILEIO_H_
#define _FILEIO_H_

#include "strclib.h"
#include <stdio.h>

/* start the Menu ------------------------------------------------------------------ */
void Menu(int mode);

/* start the file thread ----------------------------------------------------------- */
void sppFileThread();
int initNewFile(const char *filewlog, const char *filewpos);

/* Open the file/socket and proceed with decoding ---------------------------------- */
int input_oem7f(const char* filer, const char *filewlog, const char *filewpos);
int input_oem7(int8_t *Buff, int &lenD, EPOCHOBS *obs, GPSEPHREC* ephgps, BDSEPHREC* ephbds, POSRES *sppPos);

/* Output satPos and SPP information to the file/console --------------------------- */
int output_oem7f(const char *filewlog, const char *filewpos, POSRES *sppPos, EPOCHOBS *obs);
int output_oem7(EPOCHOBS *obs, POSRES *sppPos);

/* Reset sppPos results and observations ------------------------------------------- */
void reset_oem7(EPOCHOBS *obs, POSRES *sppPos);

/* print the log/final results ----------------------------------------------------- */
void PrintSatPVTlog(SATRES *obsPVT, FILE* fp);
void PrintSPPlog(POSRES *spp, FILE* fp);
void Print_spp_spv(const char *filewpos, FILE* fp, EPOCHOBS *obs, POSRES *sppPos);

#endif