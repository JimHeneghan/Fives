#ifndef _P_MACRO_H
#define _P_MACRO_H
#include <complex.h>

#define SizeX0 45
#define SizeY0 79
/****************************************************/
/*********************IMPORTANT!!********************/
/*SizeZ0 must include PML layers AND Simulationspace*/
#define SizeZ0 420
#define PMLs 150
/****************************************************/
#define MT 3e5
#define NFREQS 500
#define FRAMES 5000

#define DDx0 50e-9
#define DDy0 50e-9
#define DDz0 5e-9

#define pPMLHx_1(PP)       pPMLhx_1[PP]
#define pPMLHx_2(PP)       pPMLhx_2[PP]

#define pPMLHz_2(PP)       pPMLhz_2[PP]
#define pPMLHz_3(PP)       pPMLhz_3[PP]

#define pPMLDx_1(PP)       pPMLdx_1[PP]
#define pPMLDx_2(PP)       pPMLdx_2[PP]

#define pPMLDz_2(PP)       pPMLdz_2[PP]
#define pPMLDz_3(PP)       pPMLdz_3[PP]

#define pKRef(MM)		   pkref[MM]

#define pIXSensor(MM)      		   pixsensor[MM]

#define pPz(MM, NN, PP)     	   ppz[((MM)*SizeY0 + NN)*NFREQS + PP]
#define pPzT(MM, NN, PP)     	   ppzT[((MM)*SizeY0 + NN)*NFREQS + PP]


#define pRXSensor(MM, NN, PP)      prxsensor[((MM)*SizeY0 + NN)*NFREQS + PP]
#define pRXHSensor(MM, NN, PP)     prxhsensor[((MM)*SizeY0 + NN)*NFREQS + PP]

#define pRYSensor(MM, NN, PP)      prysensor[((MM)*SizeY0 + NN)*NFREQS + PP]
#define pRYHSensor(MM, NN, PP)     pryhsensor[((MM)*SizeY0 + NN)*NFREQS + PP]

#define pRZSensor(MM, NN, PP)      przsensor[((MM)*SizeY0 + NN)*NFREQS + PP]

#define pTXSensor(MM, NN, PP)      ptxsensor[((MM)*SizeY0 + NN)*NFREQS + PP]
#define pTXHSensor(MM, NN, PP)     ptxhsensor[((MM)*SizeY0 + NN)*NFREQS + PP]

#define pTYSensor(MM, NN, PP)      ptysensor[((MM)*SizeY0 + NN)*NFREQS + PP]
#define pTYHSensor(MM, NN, PP)     ptyhsensor[((MM)*SizeY0 + NN)*NFREQS + PP]

#define pTZSensor(MM, NN, PP)      ptzsensor[((MM)*SizeY0 + NN)*NFREQS + PP]


#define pIYSensor(MM)      piysensor[MM]

#define pIZSensor(MM)      pizsensor[MM]


#define pRXTime(MM)        prxtime[MM]
#define pTXTime(MM)        ptxtime[MM]


#define pHx(MM,NN,PP)      phx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pHy(MM,NN,PP)      phy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pHz(MM,NN,PP)      phz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pDx(MM,NN,PP)      pdx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pDy(MM,NN,PP)      pdy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pDz(MM,NN,PP)      pdz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pEx(MM,NN,PP)      pex[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pEy(MM,NN,PP)      pey[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pEz(MM,NN,PP)      pez[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pCurlHx(MM,NN,PP)      pCurlhx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlHy(MM,NN,PP)      pCurlhy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlHz(MM,NN,PP)      pCurlhz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pCurlEx(MM,NN,PP)  pCurlex[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlEy(MM,NN,PP)  pCurley[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlEz(MM,NN,PP)  pCurlez[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pICEz(MM,NN,PP)      pICez[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pICHz(MM,NN,PP)      pIChz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pSxx(MM,NN,PP)      psxx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pSyy(MM,NN,PP)      psyy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pSz(MM,NN,PP)       psz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pS1xx(MM,NN,PP)     ps1xx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS1yy(MM,NN,PP)     ps1yy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS1z(MM,NN,PP)      ps1z[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pS2xx(MM,NN,PP)     ps2xx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS2yy(MM,NN,PP)     ps2yy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS2z(MM,NN,PP)      ps2z[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pE1x(MM,NN,PP)      pe1x[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pE1y(MM,NN,PP)      pe1y[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pE1z(MM,NN,PP)      pe1z[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pMedia(MM,NN,PP)    pmedia[((MM)*SizeY0 + NN)*(SizeZ0) + PP]
#define pEpsR(MM,NN,PP) 	pepsR[((MM)*SizeY0 + NN)*(SizeZ0) + PP]

#define pExMovie(TT,MM,PP)        pexMovie[((TT)*SizeX0 + MM)*SizeZ0 + PP]
#define pExMovieXY(TT,MM,PP)      pexMoviexy[((TT)*SizeX0 + MM)*SizeY0 + PP]


#endif
