/*------------------------------------------------------------------------------
* matrix.h : SPP software socket operations
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references :  [1]"RTK_Structs.h"
*
* version : $Revision: 1.1 $ $Date: 2023/11/03 11:53:00 $
*
* history : 2023/11/03 1.0 new
*
*-----------------------------------------------------------------------------*/
#ifndef _SOCKET_H_
#define _SOCKET_H_

#include <WinSock.h>
#include <winnt.h>
#include <minwindef.h>

#pragma comment(lib,"WS2_32.lib")

/* socket ON/OFF----------------------------------------------------------- */
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port);
void CloseSocket(SOCKET& sock);

/* start the socket thread -------------------------------------------------- */
int sppSocketThread();

/* callback function */
VOID CALLBACK TimerCallback(PVOID lpParam, BOOLEAN TimerOrWaitFired);

#endif