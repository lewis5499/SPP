#include "sockets.h"
#include "strclib.h"
#include "dataIO.h"
#include "decode.h"
#include <sys/stat.h>
#include <direct.h>
#include <io.h>
#include <stdio.h>

/* Open Socket -------------------------------------------------------------- */
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
	WSADATA wsaData;
	SOCKADDR_IN addrSrv;

	if (!WSAStartup(MAKEWORD(1, 1), &wsaData))
	{
		if ((sock = socket(AF_INET, SOCK_STREAM, 0)) != INVALID_SOCKET)
		{
			addrSrv.sin_addr.S_un.S_addr = inet_addr(IP);
			addrSrv.sin_family = AF_INET;
			addrSrv.sin_port = htons(Port);
			connect(sock, (SOCKADDR*)&addrSrv, sizeof(SOCKADDR));
			return true;
		}
	}
	return false;
}

/* Close Socket ------------------------------------------------------------- */
void CloseSocket(SOCKET& sock)
{
	closesocket(sock);
	WSACleanup();
}

/* SPP mode: Socket --------------------------------------------------------- */
int sppSocketThread() 
{
	/* Create a Timer Queue */
	HANDLE g_TimerQueue = CreateTimerQueue();
	
	if (NULL == g_TimerQueue) {
		printf("Failed to create timer queue!\n");
		return -1;
	}
	
	/* Initialize params */
	TimerCallbackParams params;

	if (!OpenSocket(params.NetGps, IP1, PORT1)) {
		printf("This ip & port was not opened.\n");
		return -1;
	}

#if SAVE_REALTIME_TO_FILE
	// Running the .exe file outside the project will result in the output being 
	// directed to the 'log' folder located in the same directory as the current .exe file.
	if (access(CMAKELIST_PATH, 0) == -1) {
		if (access(OUTPATH0, 0) == -1) {
			if (mkdir(OUTPATH0)) {
				printf("Unable to create directory.\n");
			}
		}
	} else { 
	// When compiling and running the .exe file within the project, the output
	// will be directed to the 'log' folder located in the root directory of the project.
		if (access(OUTPATH, 0) == -1) {
			if (mkdir(OUTPATH)) {
				printf("Unable to create directory.\n");
			}
		}
	}
	if (initNewFile(params.filewlog, params.filewpos) < 0) { return -1;	}
#endif
	
	/* Create a Timer */
	HANDLE hTimer = NULL;
	if (!CreateTimerQueueTimer(
		&hTimer,             // Output: Timer handle
		g_TimerQueue,        // Timer queue handle
		TimerCallback,		 // Timer callback function
		(PVOID)&params,      // User data (can pass on required data)
		0,					 // Initial delay (milliseconds)
		1000,                // Interval (milliseconds)
		0                    // Flag (0 indicates that the timer is deleted 
							 //immediately after execution of the timer task)
	))  {
		printf("Failed to create timer!\n");
		return -1;
	} else {
		printf("\n[NOTICE!]:\n\n======================\nTimer started.\n\nPress Enter to exit.\n======================\n");
	}

	/* Initialize the console output */
	Sleep(1500);
	system("cls");
	printf("%s", INIT_INFO_SOKECT);
	printf("%s", FINALPOSHEADER_Console);

	/* Exit: Wait for the user to press Enter */
	getchar(); 

	/* Clean up system resources  */
	DeleteTimerQueueTimer(g_TimerQueue, hTimer, NULL);
	DeleteTimerQueue(g_TimerQueue);
	CloseSocket(params.NetGps);

	return 0;
}

/* Timer callback function -------------------------------------------------- */
VOID CALLBACK TimerCallback(PVOID lpParam, BOOLEAN TimerOrWaitFired)
{
	TimerCallbackParams* params = (TimerCallbackParams*)lpParam;

	if ((params->lenR = recv(params->NetGps, (char*)params->buff, MAXRAWLEN, 0)) > 0)
	{
		memcpy(params->Buff + params->lenD, params->buff, params->lenR);
		params->lenD += params->lenR;

	#if (SAVE_REALTIME_RAW_BINARY)
		if (params->rawbIsOpen) {
			fwrite(params->buff, sizeof(int8_t), params->lenR, params->rawb);
		}
	#endif

		do {
			params->idx = input_oem7(params->Buff, params->lenD, params->obs, params->ephgps, params->ephbds, params->sppPos);
			/* idx==0 means OBS is decoded */
			if (!params->idx) {
				calculate_oem7(params->obs, params->ephgps, params->ephbds, params->sppPos);
				output_oem7(params->obs, params->sppPos);

			#if SAVE_REALTIME_TO_FILE
				output_oem7f(params->filewlog, params->filewpos, params->sppPos, params->obs);
			#endif

				reset_oem7(params->obs, params->sppPos);
			}
		} while (!params->idx);

		memset(params->buff, 0, sizeof(int8_t) * MAXRAWLEN);
		//printf("%d\t%d\t\n", params->lenR, params->lenD);
	}
}