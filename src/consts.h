/*------------------------------------------------------------------------------
* consts.h : SPP software constants
*
*          Copyright (C) 2023 by H.Z. Liu, All rights reserved.
*
* options : none
*
* references :  [1]"RTK_Structs.h"
*
* version : $Revision: 1.1 $ $Date: 2023/10/01 13:57:00 $
*
* history : 2023/10/01 1.0 new
*                          
*-----------------------------------------------------------------------------*/

#ifndef _CONSTS_H_
#define _CONSTS_H_

#define USE_FILTER	  1

#define MODE		  1

#define MODE_FILE	  0
#define MODE_SOCKET   1

#define METHOD		  0

#define METHOD_SPP	  0
#define METHOD_RTK    1

#define SAVE_REALTIME_TO_FILE		1
#define SAVE_REALTIME_RAW_BINARY	1

#define OUTPATH			 "..\\log"
#define RAWFILE_BINARY	 "..\\log\\raw.oem719"
#define LOGFILE_RealTime "..\\log\\real-time.oem719.log"
#define POSFILE_RealTime "..\\log\\real-time.oem719.pos"
// #define DATAPATH		 "..\\log\\raw.oem719"
// #define LOGFILE_PostProc "..\\log\\raw.oem719.log"
// #define POSFILE_PostProc "..\\log\\raw.oem719.pos"
#define DATAPATH		 "..\\dataset\\202310301810.oem719"
#define LOGFILE_PostProc "..\\log\\202310301810.oem719.log"
#define POSFILE_PostProc "..\\log\\202310301810.oem719.pos"

#define OUTPATH0		  ".\\log"
#define RAWFILE_BINARY0	  ".\\log\\raw.oem719"
#define LOGFILE0_REALTIME ".\\log\\real-time.oem719.log"
#define POSFILE0_REALTIME ".\\log\\real-time.oem719.pos"
// #define DATAPATH0		  ".\\log\\raw.oem719"
// #define LOGFILE_PostProc0 ".\\log\\raw.oem719.log"
// #define POSFILE_PostProc0 ".\\log\\raw.oem719.pos"
#define DATAPATH0		  ".\\dataset\\202310301810.oem719"
#define LOGFILE_PostProc0 ".\\log\\202310301810.oem719.log"
#define POSFILE_PostProc0 ".\\log\\202310301810.oem719.pos"

#define CMAKELIST_PATH	  "..\\CMakeLists.txt"

#define ABS_V_MAX	2.0

/* kalman filter params */
//////////////////////////////////////////
/////      [Static configuration]:
////  GNSS receivers maintaining a relatively
////  static position with respect to the Earth's
////  geocentric reference system.
//init Dxk
#define InitVar		81.0
// Qk(De)
#define sigma_x		0.00010
#define sigma_y		0.00010
#define sigma_z		0.00010
#define sigma_bG	0.00030
#define sigma_bC	0.00030
#define sigma_d		0.00010
// Rk(Dr)
#define var_obs	    16.0
//////////////////////////////////////////

//////////////////////////////////////////
/////      [Kinematic configuration]:
////  Pedestrians moving at a constant speed
/*//init Dxk
#define InitVar		9.0
// Qk(De)
#define sigma_x		1.5
#define sigma_y		1.5
#define sigma_z		1.5
#define sigma_bG	1.5
#define sigma_bC	1.5
#define sigma_d		1.5
// Rk(Dr)
#define var_obs	    81.0*/
//////////////////////////////////////////

typedef signed char        int8_t;
typedef short              int16_t;
typedef int                int32_t;
typedef long long          int64_t;
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64_t;

/* OEM basic info */
#define OEM4SYNC1       0xAA	     /* oem7/6/4 message start sync code 1 */
#define OEM4SYNC2       0x44		 /* oem7/6/4 message start sync code 2 */
#define OEM4SYNC3       0x12		 /* oem7/6/4 message start sync code 3 */
#define OEM4SYNCLEN     3
#define OEM4HLEN        28			 /* oem7/6/4 message header length (bytes) */
#define CRC32LEN		4
#define POLYCRC32       0xEDB88320u  /* CRC32 polynomial */

/* message IDs */
#define ID_RANGECMP			 140     /* oem7/6/4 range compressed */
#define ID_RANGE			  43     /* oem7/6/4 range measurement */
#define ID_RAWEPHEM			  41     /* oem7/6/4 raw ephemeris */
#define ID_IONUTC			   8     /* oem7/6/4 iono and utc data */
#define ID_RAWWAASFRAME		 287     /* oem7/6/4 raw waas frame */
#define ID_RAWSBASFRAME		 973     /* oem7/6 raw sbas frame */
#define ID_GLOEPHEMERIS		 723     /* oem7/6/4 glonass ephemeris */
#define ID_GALEPHEMERIS		1122     /* oem7/6 decoded galileo ephemeris */
#define ID_GALIONO			1127     /* oem7/6 decoded galileo iono corrections */
#define ID_GALCLOCK			1121     /* oem7/6 galileo clock information */
#define ID_QZSSRAWEPHEM		1331     /* oem7/6 qzss raw ephemeris */
#define ID_QZSSRAWSUBFRAME  1330     /* oem7/6 qzss raw subframe */
#define ID_QZSSIONUTC		1347     /* oem7/6 qzss ion/utc parameters */
#define ID_BDSEPHEM			1696     /* oem7/6 decoded bds ephemeris */
#define ID_NAVICEPHEMERIS   2123     /* oem7 decoded navic ephemeris */
#define ID_GPSEPHEM			7		 /* oem7 decoded GPS L1 C/A ephemerides */
#define ID_BESTPOS			42		 /* oem7 decoded best position */

/* code IDs */
#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless,B1codeless (GPS,BDS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A,B1A    (GAL,BDS) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1S (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2bI  (GAL,BDS) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2bQ  (GAL,BDS) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
#define CODE_L6A    30                  /* obs code: E6A,B3A    (GAL,BDS) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
#define CODE_L6S    35                  /* obs code: L6S        (QZS) */
#define CODE_L6L    36                  /* obs code: L6L        (QZS) */
#define CODE_L8I    37                  /* obs code: E5abI      (GAL) */
#define CODE_L8Q    38                  /* obs code: E5abQ      (GAL) */
#define CODE_L8X    39                  /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
#define CODE_L2I    40                  /* obs code: B1_2I      (BDS) */
#define CODE_L2Q    41                  /* obs code: B1_2Q      (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) (obsolute) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) (obsolute) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define CODE_L1D    56                  /* obs code: B1D        (BDS) */
#define CODE_L5D    57                  /* obs code: L5D(L5S),B2aD (QZS,BDS) */
#define CODE_L5P    58                  /* obs code: L5P(L5S),B2aP (QZS,BDS) */
#define CODE_L5Z    59                  /* obs code: L5D+P(L5S) (QZS) */
#define CODE_L6E    60                  /* obs code: L6E        (QZS) */
#define CODE_L7D    61                  /* obs code: B2bD       (BDS) */
#define CODE_L7P    62                  /* obs code: B2bP       (BDS) */
#define CODE_L7Z    63                  /* obs code: B2bD+P     (BDS) */
#define CODE_L8D    64                  /* obs code: B2abD      (BDS) */
#define CODE_L8P    65                  /* obs code: B2abP      (BDS) */
#define CODE_L4A    66                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4B    67                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4X    68                  /* obs code: G1al1OCd+p (GLO) */
#define MAXCODE     68                  /* max number of obs code */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_BDS     0x20                /* navigation system: BeiDou */
#define SYS_IRN     0x40                /* navigation system: IRNS */
#define SYS_LEO     0x80                /* navigation system: LEO */
#define SYS_ALL     0xFF                /* navigation system: all */

/* mathmatic functions */
#define SQR(x)		((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

/* Hopefield Model */
#define H0			0.0					    /* m */
#define T0			(15.0+273.16)			/* K */
#define p0			1013.25				    /* mbar */
#define RH0			0.5					

/* COMMON UINTS */
#define EPSILON		1e-10			        /* A very small threshold close to 0 */					
#define PI			3.1415926535898		    /* PI */
#define Rad			(PI/180.0)              /* Radians per degree */
#define Deg			(180.0/PI)              /* Degrees per radian */
#define CLIGHT		299792458.0             /* Speed of light  [m/s]; IAU 1976  */

/* SPP params */
#define RTOL_KEPLER		 1E-13				    /* relative tolerance for Kepler equation */
#define MAX_ITER_KEPLER  30				        /* max number of iteration of Kelpler */
#define SIN_5			 -0.0871557427476582	/* sin(-5.0 deg) */
#define COS_5			 0.9961946980917456	    /* cos(-5.0 deg) */
#define THRES_tk		 7500
#define THRES_dLGF		 0.05
#define THRES_dLMW	     3.0
#define THRES_El		 10.0
#define LS_THRES		 1e-6
#define MAXITERS		 10
#define MAX_dPos_NORM    1e-4
#define GPST_BDT		 14.0			        /* Difference between GPS time and Beidou time[s] */
#define GPST_BDT_WEEKS   1356
#define MAXCHANNUM		 95
#define MAXGPSNUM		 32
#define MAXBDSNUM		 63
#define MAXRAWLEN		 40960			        /* max length of receiver raw message */
#define MAXONERAWLEN	 16384
#define BDSEPHEMLEN		 196
#define GPSEPHEMLEN		 224
#define MAXSTACKLEN		 16384

/* BASE ref POS */
// 2020 ref
#define rcvXref1   -2267794.937
#define rcvYref1	5009345.236
#define rcvZref1	3220980.312
// 2023 ref
#define rcvXref    -2267335.6694
#define rcvYref		5008649.1555
#define rcvZref		3222374.9736

/* ECEF params */
#define R_WGS84		6378137.0               /* Radius Earth [m]; WGS-84  */
#define F_WGS84		1.0/298.257223563       /* Flattening; WGS-84   */
#define Omega_WGS	7.2921151467e-5         /*[rad/s], the earth rotation rate */
#define GM_WGS84	398600.5e+9             /* [m^3/s^2]; WGS-84 */

#define R_CGCS2000  6378137.0               /* Radius Earth [m]; CGCS2000  */
#define F_CGCS2000  1.0/298.257222101       /* Flattening; CGCS2000   */
#define Omega_CGCS  7.2921150e-5            /*[rad/s], the earth rotation rate */
#define GM_CGCS     398600.4418e+9          /* [m^3/s^2]; CGCS2000  */

/* some constants about GPS satellite signal */
#define  FG1_GPS    1575.42E6               /* L1 signal frequency */
#define  FG2_GPS    1227.60E6               /* L2 signal frequency */
#define  FG12R      (77/60.0)               /* FG1_Freq/FG2_Freq */
#define  WL1_GPS	(CLIGHT/FG1_GPS)		/* L1 signal wavelength */
#define  WL2_GPS    (CLIGHT/FG2_GPS)		/* L2 signal wavelength */

/* some constants about Compass satellite signal */
#define  FG1_BDS	1561.098E6              /* B1 signal frequency */
#define  FG2_BDS    1207.140E6              /* B2 signal frequency */
#define  FG3_BDS    1268.520E6              /* B2 signal frequency */
#define  FC12R      (763/590.0)             /* B1 signal frequency/B2 signal frequency */
#define  FC13R      (763/620.0)				/* B1 signal frequency/B3 signal frequency */
#define  WL1_BDS    (CLIGHT/FG1_BDS)		/* B1 signal wavelength */
#define  WL2_BDS    (CLIGHT/FG2_BDS)		/* B2 signal wavelength */
#define  WL3_BDS    (CLIGHT/FG3_BDS)		/* B3 signal wavelength */

/* SOCKECT params */
#define IP1		    "47.114.134.129"
#define PORT1	    7190
#define IP2		    "8.140.46.126"
#define PORT2	    5002

/* HEADER to be printed */
#define FINALPOSHEADER_PosFile " Wk     SOW         ECEF-X/m       ECEF-Y/m       ECEF-Z/m    \
REF-ECEF-X/m    REF-ECEF-Y/m   REF-ECEF-Z/m   EAST/m  NORTH/m  UP/m       B/deg          \
L/deg        H/m      VX/m     VY/m     VZ/m    PDOP    SigmaP   SigmaV  GS  BS  n"

#define FINALPOSHEADER_Console " Wk     SOW        ECEF-X/m       ECEF-Y/m      ECEF-Z/m      \
B/deg    L/deg     H/m       E/m      N/m      U/m      VX/m     VY/m    VZ/m    PDOP    SigmaP   SigmaV  GS BS nSats\n"

#define INIT_INFO_SOKECT "/*------------------------------------------------------------------------------\n\
*          [REAL-TIME SPP software] main program\n\
*\n\
*          Copyright(C) 2023 by H.Z.Liu, All rights reserved.\n\
*\n\
* options : SOCKET\n\
*\n\
* references : [1] RTK_Structs.h	[2] rtklib \n\
*\n\
* version : $Revision: 1.6 $ $Date: 2023 / 12 / 15   18:29:52 $\n\
*\n\
* history : 2023 / 12 / 15   1.6 new	\n\
*           2023 / 12 / 14   1.5    	\n\
*           2023 / 12 / 11   1.4		\n\
*           2023 / 11 / 24   1.3		\n\
*           2023 / 11 / 06   1.2		\n\
*           2023 / 11 / 05   1.1		\n\
*           2023 / 11 / 03   1.0		\n\
*\n\
* Exit Program:  Press Enter to exit (KILL TIMER).\n\
*\n\
*------------------------------------------------------------------------------ */ \n\n"

#define INIT_INFO_FILE "/*------------------------------------------------------------------------------\n\
* SPP software main program\n\
*\n\
*          Copyright(C) 2023 by H.Z.Liu, All rights reserved.\n\
*\n\
* options : FILE\n\
*\n\
* references : [1] RTK_Structs.h	[2] rtklib \n\
*\n\
* version : $Revision: 1.6 $ $Date: 2023 / 12 / 15   18:29:52 $\n\
*\n\
* history : 2023 / 12 / 15   1.6 new	\n\
*           2023 / 12 / 14   1.5    	\n\
*           2023 / 12 / 11   1.4		\n\
*           2023 / 11 / 24   1.3		\n\
*           2023 / 11 / 06   1.2		\n\
*           2023 / 11 / 05   1.1		\n\
*           2023 / 11 / 03   1.0		\n\
*\n\
*------------------------------------------------------------------------------ */ \n\n"
#endif