#ifndef _FDTD_MACRO_H
#define _FDTD_MACRO_H
#include <complex.h>
#include "fdtd-grid.h"
#include "p-macro.h"

/* macros that permit the "Grid" to be specified */
/* one-dimensional grid */
#define Hy1G(G,MM)     G->hy[MM]
#define Chyh1G(G,MM)   G->chyh[MM]
#define Chye1G(G,MM)   G->chye[MM]

#define Ez1G(G,MM)     G->ez[MM]
#define Ceze1G(G,MM)   G->ceze[MM]
#define Cezh1G(G,MM)   G->cezh[MM]

/* TMz grid */
#define Hx2G(G,MM,NN)   G->hx[(MM)*(SizeYG(G)) + NN]
#define Chxh2G(G,MM,NN) G->chxh[(MM)*(SizeYG(G)) + NN]
#define Chxe2G(G,MM,NN) G->chxe[(MM)*(SizeYG(G)) + NN]

#define Hy2G(G,MM,NN)   G->hy[(MM)*SizeYG(G) + NN]
#define Chyh2G(G,MM,NN) G->chyh[(MM)*SizeYG(G) + NN]
#define Chye2G(G,MM,NN) G->chye[(MM)*SizeYG(G) + NN]

#define Ez2G(G,MM,NN)   G->ez[(MM)*SizeYG(G) + NN]
#define Ceze2G(G,MM,NN) G->ceze[(MM)*SizeYG(G) + NN]
#define Cezh2G(G,MM,NN) G->cezh[(MM)*SizeYG(G) + NN]

/* TEz grid */
#define Ex2G(G,MM,NN)   G->ex[(MM)*SizeYG(G) + NN]
#define Cexe2G(G,MM,NN) G->cexe[(MM)*SizeYG(G) + NN]
#define Cexh2G(G,MM,NN) G->cexh[(MM)*SizeYG(G) + NN]

#define Ey2G(G,MM,NN)   G->ey[(MM)*(SizeYG(G)) + NN]
#define Ceye2G(G,MM,NN) G->ceye[(MM)*(SizeYG(G)) + NN]
#define Ceyh2G(G,MM,NN) G->ceyh[(MM)*(SizeYG(G)) + NN]

#define Hz2G(G,MM,NN)   G->hz[(MM)*(SizeYG(G)) + NN]
#define Chzh2G(G,MM,NN) G->chzh[(MM)*(SizeYG(G)) + NN]
#define Chze2G(G,MM,NN) G->chze[(MM)*(SizeYG(G)) + NN]

/* 3D grid */


#define HxG(G,MM,NN,PP)   G->hx[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define ChxhG(G,MM,NN,PP) G->chxh[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define ChxeG(G,MM,NN,PP) G->chxe[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

#define HyG(G,MM,NN,PP)   G->hy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define ChyhG(G,MM,NN,PP) G->chyh[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define ChyeG(G,MM,NN,PP) G->chye[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define HzG(G,MM,NN,PP)   G->hz[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]
#define ChzhG(G,MM,NN,PP) G->chzh[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]
#define ChzeG(G,MM,NN,PP) G->chze[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]

#define ExG(G,MM,NN,PP)   G->ex[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define CexeG(G,MM,NN,PP) G->cexe[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define CexhG(G,MM,NN,PP) G->cexh[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]

#define EyG(G,MM,NN,PP)   G->ey[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]
#define CeyeG(G,MM,NN,PP) G->ceye[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]
#define CeyhG(G,MM,NN,PP) G->ceyh[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]

#define EzG(G,MM,NN,PP)   G->ez[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define CezeG(G,MM,NN,PP) G->ceze[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define CezhG(G,MM,NN,PP) G->cezh[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define E1xG(G,MM,NN,PP)   G->e1x[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define E1yG(G,MM,NN,PP)   G->e1y[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]
#define E1zG(G,MM,NN,PP)   G->e1z[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define DxG(G,MM,NN,PP)   G->dx[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define DyG(G,MM,NN,PP)   G->dy[((MM)*(SizeYG(G)) + NN)*SizeZG(G) + PP]
#define DzG(G,MM,NN,PP)   G->dz[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define IxxG(G,MM,NN,PP)  G->ixx[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define IxyG(G,MM,NN,PP)  G->ixy[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]

#define SxxG(G,MM,NN,PP)  G->sxx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S1xxG(G,MM,NN,PP) G->s1xx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S2xxG(G,MM,NN,PP) G->s2xx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define SxyG(G,MM,NN,PP)  G->sxy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S1xyG(G,MM,NN,PP) G->s1xy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S2xyG(G,MM,NN,PP) G->s2xy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define QxxG(G,MM,NN,PP)  G->qxx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q1xxG(G,MM,NN,PP) G->q1xx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q2xxG(G,MM,NN,PP) G->q2xx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define QxyG(G,MM,NN,PP)  G->qxy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q1xyG(G,MM,NN,PP) G->q1xy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q2xyG(G,MM,NN,PP) G->q2xy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]


#define IyyG(G,MM,NN,PP)  G->iyy[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define IyxG(G,MM,NN,PP)  G->iyx[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]

#define SyyG(G,MM,NN,PP)  G->syy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S1yyG(G,MM,NN,PP) G->s1yy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S2yyG(G,MM,NN,PP) G->s2yy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define SyxG(G,MM,NN,PP)  G->syx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S1yxG(G,MM,NN,PP) G->s1yx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S2yxG(G,MM,NN,PP) G->s2yx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define QyyG(G,MM,NN,PP)  G->qyy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q1yyG(G,MM,NN,PP) G->q1yy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q2yyG(G,MM,NN,PP) G->q2yy[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define QyxG(G,MM,NN,PP)  G->qyx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q1yxG(G,MM,NN,PP) G->q1yx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Q2yxG(G,MM,NN,PP) G->q2yx[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define IzG(G,MM,NN,PP)   G->iz[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define SzG(G,MM,NN,PP)   G->sz[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S1zG(G,MM,NN,PP) G->s1z[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define S2zG(G,MM,NN,PP) G->s2z[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define Omega_pG(G,MM,NN,PP) G->omega_p[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Omega_bG(G,MM,NN,PP) G->omega_b[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define Nu_cG(G,MM,NN,PP) G->nu_c[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
#define SigmaG(G,MM,NN,PP) G->sigma[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define EpsRG(G,MM,NN,PP) G->epsR[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

#define MediaG(G,MM,NN,PP) G->media[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

/*PML arrays*/
/*only arrays in the 1D Z direction PML are uncommented*/


#define SigmaZDG(G,PP)	 G->sigmazD[PP]
#define SigmaZHG(G,PP)	 G->sigmazH[PP]

#define PMLHx_0G(G,PP)   G->PMLhx_0[PP]
#define PMLHx_1G(G,PP)   G->PMLhx_1[PP]
#define PMLHx_2G(G,PP)   G->PMLhx_2[PP]

#define PMLHz_2G(G,PP)   G->PMLhz_2[PP]
#define PMLHz_3G(G,PP)   G->PMLhz_3[PP]

#define PMLDx_0G(G,PP)   G->PMLdx_0[PP]
#define PMLDx_1G(G,PP)   G->PMLdx_1[PP]
#define PMLDx_2G(G,PP)   G->PMLdx_2[PP]

#define PMLDz_2G(G,PP)   G->PMLdz_2[PP]
#define PMLDz_3G(G,PP)   G->PMLdz_3[PP]

// #define SigmaXG(G,MM,NN,PP) G->sigmax[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]
// #define SigmaYG(G,MM,NN,PP) G->sigmay[((MM)*SizeYG(G) + NN)*(SizeZG(G)) + PP]

// #define PMLHy_0G(G,MM,NN,PP)   G->PMLhy_0[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLHz_0G(G,MM,NN,PP)   G->PMLhz_0[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]


// #define PMLHy_1G(G,MM,NN,PP)   G->PMLhy_1[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLHz_1G(G,MM,NN,PP)   G->PMLhz_1[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLHy_2G(G,MM,NN,PP)   G->PMLhy_2[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLHx_3G(G,MM,NN,PP)   G->PMLhx_3[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLHy_3G(G,MM,NN,PP)   G->PMLhy_3[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLHx_4G(G,MM,NN,PP)   G->PMLhx_4[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLHy_4G(G,MM,NN,PP)   G->PMLhy_4[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLHz_4G(G,MM,NN,PP)   G->PMLhz_4[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLDy_0G(G,MM,NN,PP)   G->PMLdy_0[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLDz_0G(G,MM,NN,PP)   G->PMLdz_0[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLDy_1G(G,MM,NN,PP)   G->PMLdy_1[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLDz_1G(G,MM,NN,PP)   G->PMLdz_1[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLDy_2G(G,MM,NN,PP)   G->PMLdy_2[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLDx_3G(G,MM,NN,PP)   G->PMLdx_3[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLDy_3G(G,MM,NN,PP)   G->PMLdy_3[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

// #define PMLDx_4G(G,MM,NN,PP)   G->PMLdx_4[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLDy_4G(G,MM,NN,PP)   G->PMLdy_4[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
// #define PMLDz_4G(G,MM,NN,PP)   G->PMLdz_4[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

/*PML integration terms*/
#define ICHxG(G,MM,NN,PP)   G->IChx[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define ICHyG(G,MM,NN,PP)   G->IChy[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define ICHzG(G,MM,NN,PP)   G->IChz[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

#define ICExG(G,MM,NN,PP)   G->ICex[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define ICEyG(G,MM,NN,PP)   G->ICey[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define ICEzG(G,MM,NN,PP)   G->ICez[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

/*Curl arrays*/
#define CurlHxG(G,MM,NN,PP)   G->Curlhx[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define CurlHyG(G,MM,NN,PP)   G->Curlhy[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define CurlHzG(G,MM,NN,PP)   G->Curlhz[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

#define CurlExG(G,MM,NN,PP)   G->Curlex[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define CurlEyG(G,MM,NN,PP)   G->Curley[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]
#define CurlEzG(G,MM,NN,PP)   G->Curlez[((MM)*(SizeYG(G)) + NN)*(SizeZG(G)) + PP]

/* Frequency Sensors */
#define IXSensorG(G,MM)    		   G->ixsensor[MM]

#define RXSensorG(G,MM, NN, PP)    G->rxsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]
#define RXHSensorG(G,MM, NN, PP)   G->rxhsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]

#define RYSensorG(G,MM, NN, PP)    G->rysensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]
#define RYHSensorG(G,MM, NN, PP)   G->ryhsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]

#define RZSensorG(G,MM, NN, PP)    G->rzsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]


#define TXSensorG(G,MM, NN, PP)    G->txsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]
#define TXHSensorG(G,MM, NN, PP)   G->txhsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]

#define TYSensorG(G,MM, NN, PP)    G->tysensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]
#define TYHSensorG(G,MM, NN, PP)   G->tyhsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]

#define TZSensorG(G,MM, NN, PP)    G->tzsensor[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]

#define IYSensorG(G,MM)    		   G->iysensor[MM]

#define IZSensorG(G,MM)    		   G->izsensor[MM]

#define PzG(G,MM, NN, PP)    	   G->pz[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]
#define PzTG(G,MM, NN, PP)    	   G->pzT[((MM)*(SizeYG(G)) + NN)*(NFREQS) + PP]


#define KRefG(G,MM)	   G-> kref[MM]

/*Time dependent sensors*/
#define RXTimeG(G,MM)  G->rxtime[MM]
#define TXTimeG(G,MM)  G->txtime[MM]

/* Movie Sensor */
#define ExMovieG(G,TT,MM,PP)     G->exMovie[((TT)*SizeXG(G) + MM)*SizeZG(G) + PP]
#define ExMovieXYG(G,TT,MM,PP)   G->exMoviexy[((TT)*SizeXG(G) + MM)*SizeYG(G) + PP]


#define SizeXG(G)      G->sizeX
#define SizeYG(G)      G->sizeY
#define SizeZG(G)      G->sizeZ
#define TimeG(G)       G->time
#define MaxTimeG(G)    G->maxTime
#define CdtdsG(G)      G->cdtds
#define TypeG(G)       G->type
#define PhixG(G)       G->phix
#define PhiyG(G)       G->phiy
#define PhizG(G)       G->phiz


/* macros that assume the "Grid" is "g" */
/* one-dimensional grid */
#define Hy1(MM)     Hy1G(g,MM)
#define Chyh1(MM)   Chyh1G(g,MM)   
#define Chye1(MM)   Chye1G(g,MM)   

#define Ez1(MM)     Ez1G(g,MM)     
#define Ceze1(MM)   Ceze1G(g,MM)   
#define Cezh1(MM)   Cezh1G(g,MM)   

/* TMz grid */
#define Hx2(MM,NN)   Hx2G(g,MM,NN)   
#define Chxh2(MM,NN) Chxh2G(g,MM,NN) 
#define Chxe2(MM,NN) Chxe2G(g,MM,NN) 

#define Hy2(MM,NN)   Hy2G(g,MM,NN)   
#define Chyh2(MM,NN) Chyh2G(g,MM,NN) 
#define Chye2(MM,NN) Chye2G(g,MM,NN) 

#define Ez2(MM,NN)   Ez2G(g,MM,NN)   
#define Ceze2(MM,NN) Ceze2G(g,MM,NN) 
#define Cezh2(MM,NN) Cezh2G(g,MM,NN) 

/* TEz grid */
#define Hz2(MM,NN)   Hz2G(g,MM,NN)   
#define Chzh2(MM,NN) Chzh2G(g,MM,NN) 
#define Chze2(MM,NN) Chze2G(g,MM,NN) 

#define Ex2(MM,NN)   Ex2G(g,MM,NN)   
#define Cexe2(MM,NN) Cexe2G(g,MM,NN) 
#define Cexh2(MM,NN) Cexh2G(g,MM,NN) 

#define Ey2(MM,NN)   Ey2G(g,MM,NN)   
#define Ceye2(MM,NN) Ceye2G(g,MM,NN) 
#define Ceyh2(MM,NN) Ceyh2G(g,MM,NN) 

/* 3D grid */



#define Hx(MM,NN,PP)   HxG(g,MM,NN,PP)   
#define Chxh(MM,NN,PP) ChxhG(g,MM,NN,PP) 
#define Chxe(MM,NN,PP) ChxeG(g,MM,NN,PP) 

#define Hy(MM,NN,PP)   HyG(g,MM,NN,PP)   
#define Chyh(MM,NN,PP) ChyhG(g,MM,NN,PP) 
#define Chye(MM,NN,PP) ChyeG(g,MM,NN,PP) 

#define Hz(MM,NN,PP)   HzG(g,MM,NN,PP)   
#define Chzh(MM,NN,PP) ChzhG(g,MM,NN,PP) 
#define Chze(MM,NN,PP) ChzeG(g,MM,NN,PP) 

#define Ex(MM,NN,PP)   ExG(g,MM,NN,PP)   
#define Cexe(MM,NN,PP) CexeG(g,MM,NN,PP) 
#define Cexh(MM,NN,PP) CexhG(g,MM,NN,PP) 

#define Ey(MM,NN,PP)   EyG(g,MM,NN,PP)   
#define Ceye(MM,NN,PP) CeyeG(g,MM,NN,PP) 
#define Ceyh(MM,NN,PP) CeyhG(g,MM,NN,PP) 

#define Ez(MM,NN,PP)   EzG(g,MM,NN,PP)
#define Ceze(MM,NN,PP) CezeG(g,MM,NN,PP)
#define Cezh(MM,NN,PP) CezhG(g,MM,NN,PP) 

#define E1x(MM,NN,PP)   E1xG(g,MM,NN,PP)
#define E1y(MM,NN,PP)   E1yG(g,MM,NN,PP) 
#define E1z(MM,NN,PP)   E1zG(g,MM,NN,PP)

#define Dx(MM,NN,PP)   DxG(g,MM,NN,PP) 
#define Dy(MM,NN,PP)   DyG(g,MM,NN,PP)
#define Dz(MM,NN,PP)   DzG(g,MM,NN,PP)

#define Ixx(MM,NN,PP)   IxxG(g,MM,NN,PP)
#define Ixy(MM,NN,PP)   IxyG(g,MM,NN,PP)

#define Sxx(MM,NN,PP)   SxxG(g,MM,NN,PP)
#define S1xx(MM,NN,PP)  S1xxG(g,MM,NN,PP)
#define S2xx(MM,NN,PP)  S2xxG(g,MM,NN,PP)

#define Sxy(MM,NN,PP)   SxyG(g,MM,NN,PP)
#define S1xy(MM,NN,PP)  S1xyG(g,MM,NN,PP)
#define S2xy(MM,NN,PP)  S2xyG(g,MM,NN,PP)

#define Qxx(MM,NN,PP)   QxxG(g,MM,NN,PP)
#define Q1xx(MM,NN,PP)  Q1xxG(g,MM,NN,PP)
#define Q2xx(MM,NN,PP)  Q2xxG(g,MM,NN,PP)

#define Qxy(MM,NN,PP)   QxyG(g,MM,NN,PP)
#define Q1xy(MM,NN,PP)  Q1xyG(g,MM,NN,PP)
#define Q2xy(MM,NN,PP)  Q2xyG(g,MM,NN,PP)

#define Iyy(MM,NN,PP)   IyyG(g,MM,NN,PP)
#define Iyx(MM,NN,PP)   IyxG(g,MM,NN,PP)

#define Syy(MM,NN,PP)   SyyG(g,MM,NN,PP)
#define S1yy(MM,NN,PP)  S1yyG(g,MM,NN,PP)
#define S2yy(MM,NN,PP)  S2yyG(g,MM,NN,PP)

#define Syx(MM,NN,PP)   SyxG(g,MM,NN,PP)
#define S1yx(MM,NN,PP)  S1yxG(g,MM,NN,PP)
#define S2yx(MM,NN,PP)  S2yxG(g,MM,NN,PP)

#define Qyy(MM,NN,PP)   QyyG(g,MM,NN,PP)
#define Q1yy(MM,NN,PP)  Q1yyG(g,MM,NN,PP)
#define Q2yy(MM,NN,PP)  Q2yyG(g,MM,NN,PP)

#define Qyx(MM,NN,PP)   QyxG(g,MM,NN,PP)
#define Q1yx(MM,NN,PP)  Q1yxG(g,MM,NN,PP)
#define Q2yx(MM,NN,PP)  Q2yxG(g,MM,NN,PP)

#define Iz(MM,NN,PP)   IzG(g,MM,NN,PP)
#define Sz(MM,NN,PP)   SzG(g,MM,NN,PP)
#define S1z(MM,NN,PP)  S1zG(g,MM,NN,PP)
#define S2z(MM,NN,PP)  S2zG(g,MM,NN,PP)

#define EpsR(MM,NN,PP) EpsRG(g,MM,NN,PP)
#define Omega_p(MM,NN,PP) Omega_pG(g,MM,NN,PP)
#define Nu_c(MM,NN,PP) Nu_cG(g,MM,NN,PP)
#define Omega_b(MM,NN,PP) Omega_bG(g,MM,NN,PP)
#define Sigma(MM,NN,PP) SigmaG(g,MM,NN,PP)
#define Media(MM,NN,PP) MediaG(g,MM,NN,PP)

/*Curl arrays*/
#define CurlHx(MM,NN,PP)   CurlHxG(g,MM,NN,PP)
#define CurlHy(MM,NN,PP)   CurlHyG(g,MM,NN,PP)
#define CurlHz(MM,NN,PP)   CurlHzG(g,MM,NN,PP)


#define CurlEx(MM,NN,PP)   CurlExG(g,MM,NN,PP)
#define CurlEy(MM,NN,PP)   CurlEyG(g,MM,NN,PP)
#define CurlEz(MM,NN,PP)   CurlEzG(g,MM,NN,PP)

/*PML arrays*/
// #define SigmaX(MM,NN,PP) SigmaXG(g,MM,NN,PP)
// #define SigmaY(MM,NN,PP) SigmaYG(g,MM,NN,PP)
#define SigmaZD(PP)	  SigmaZDG(g,PP)
#define SigmaZH(PP)	  SigmaZHG(g,PP)

#define PMLHx_0(PP)   PMLHx_0G(g,PP) 
#define PMLHx_1(PP)   PMLHx_1G(g,PP) 
#define PMLHx_2(PP)   PMLHx_2G(g,PP) 

#define PMLHz_2(PP)   PMLHz_2G(g,PP)
#define PMLHz_3(PP)   PMLHz_3G(g,PP)

#define PMLDx_0(PP)   PMLDx_0G(g,PP)
#define PMLDx_1(PP)   PMLDx_1G(g,PP) 
#define PMLDx_2(PP)   PMLDx_2G(g,PP)

#define PMLDz_2(PP)   PMLDz_2G(g,PP)
#define PMLDz_3(PP)   PMLDz_3G(g,PP)
// #define PMLHy_0(MM,NN,PP)   PMLHy_0G(g,MM,NN,PP)
// #define PMLHz_0(MM,NN,PP)   PMLHz_0G(g,MM,NN,PP)

// #define PMLHy_1(MM,NN,PP)   PMLHy_1G(g,MM,NN,PP)
// #define PMLHz_1(MM,NN,PP)   PMLHz_1G(g,MM,NN,PP)

// #define PMLHy_2(MM,NN,PP)   PMLHy_2G(g,MM,NN,PP)

// #define PMLHx_3(MM,NN,PP)   PMLHx_3G(g,MM,NN,PP) 
// #define PMLHy_3(MM,NN,PP)   PMLHy_3G(g,MM,NN,PP)

// #define PMLHx_4(MM,NN,PP)   PMLHx_4G(g,MM,NN,PP) 
// #define PMLHy_4(MM,NN,PP)   PMLHy_4G(g,MM,NN,PP)
// #define PMLHz_4(MM,NN,PP)   PMLHz_4G(g,MM,NN,PP)

// #define PMLDy_0(MM,NN,PP)   PMLDy_0G(g,MM,NN,PP)
// #define PMLDz_0(MM,NN,PP)   PMLDz_0G(g,MM,NN,PP)

// #define PMLDy_1(MM,NN,PP)   PMLDy_1G(g,MM,NN,PP)
// #define PMLDz_1(MM,NN,PP)   PMLDz_1G(g,MM,NN,PP)
 
// #define PMLDy_2(MM,NN,PP)   PMLDy_2G(g,MM,NN,PP)

// #define PMLDx_3(MM,NN,PP)   PMLDx_3G(g,MM,NN,PP) 
// #define PMLDy_3(MM,NN,PP)   PMLDy_3G(g,MM,NN,PP)

// #define PMLDx_4(MM,NN,PP)   PMLDx_4G(g,MM,NN,PP) 
// #define PMLDy_4(MM,NN,PP)   PMLDy_4G(g,MM,NN,PP)
// #define PMLDz_4(MM,NN,PP)   PMLDz_4G(g,MM,NN,PP)

/*PML integration arrays*/
#define ICHx(MM,NN,PP)   ICHxG(g,MM,NN,PP)
#define ICHy(MM,NN,PP)   ICHyG(g,MM,NN,PP)
#define ICHz(MM,NN,PP)   ICHzG(g,MM,NN,PP)

#define ICEx(MM,NN,PP)   ICExG(g,MM,NN,PP)
#define ICEy(MM,NN,PP)   ICEyG(g,MM,NN,PP)
#define ICEz(MM,NN,PP)   ICEzG(g,MM,NN,PP)

 
/* Frequency Sensors */
#define IXSensor(MM)     		  IXSensorG(g,MM)

#define RXSensor(MM, NN, PP)      RXSensorG(g,MM, NN, PP)
#define RYSensor(MM, NN, PP)      RYSensorG(g,MM, NN, PP)
#define RZSensor(MM, NN, PP)      RZSensorG(g,MM, NN, PP)

#define RYHSensor(MM, NN, PP)     RYHSensorG(g,MM, NN, PP)
#define RXHSensor(MM, NN, PP)     RXHSensorG(g,MM, NN, PP)

#define TXSensor(MM, NN, PP)      TXSensorG(g,MM, NN, PP)
#define TYSensor(MM, NN, PP)      TYSensorG(g,MM, NN, PP)
#define TZSensor(MM, NN, PP)      TZSensorG(g,MM, NN, PP)


// #define TXHSensor(MM, NN, PP)     TXHSensorG(g,MM, NN, PP)

// #define TYHSensor(MM, NN, PP)     TYHSensorG(g,MM, NN, PP)


// #define IYSensor(MM)     		   IYSensorG(g,MM)
// #define TYSensor(MM)     		   TYSensorG(g,MM)

// #define RZSensor(MM)     		   RZSensorG(g,MM)
// #define IZSensor(MM)     		   IZSensorG(g,MM)
// #define TZSensor(MM)     		   TZSensorG(g,MM)


#define Pz(MM,NN,PP)   PzG(g,MM,NN,PP)
#define PzT(MM,NN,PP)  PzTG(g,MM,NN,PP)


#define KRef(MM)    KRefG(g,MM)		   

/*Time dependent sensors*/
#define RXTime(MM)  RXTimeG(g,MM)
#define TXTime(MM)  TXTimeG(g,MM)

/* Movie Sensor */
#define ExMovie(TT,MM,PP)     ExMovieG(g, TT,MM,PP) 
#define ExMovieXY(TT,MM,PP)   ExMovieXYG(g, TT,MM,PP) 


#define SizeX       SizeXG(g)
#define SizeY       SizeYG(g)
#define SizeZ       SizeZG(g)
#define Time        TimeG(g)
#define MaxTime     MaxTimeG(g)   
#define Cdtds       CdtdsG(g)
#define Type        TypeG(g)
#define Phix        PhixG(g)
#define Phiy        PhiyG(g)
#define Phiz        PhizG(g)

#endif
