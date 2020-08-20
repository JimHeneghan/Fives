#include "ezinc.h"
#include "p-macro.h"
#include <openacc.h>
#include <complex.h>
#include <accelmath.h>
static double cdtds, ppw = 0;

/* initialize source-function variables */
void ezIncInit(Grid *g){

  ppw = 30;
  cdtds = Cdtds;
  return;
}

/* calculate source function at given time and location */
void TFSF_E(double time, Grid *g) {
  double arg;
  double rick, sig, dt, ddx, t0, tau, Gauss, Gauss1, c0, ExSr, Eysr, Nlam;
  int mm, nn, pp;

  int tc = SizeX0*SizeY0*SizeZ0;
  /*complex*/ double *pCurlex = g->Curlex;
  /*complex*/ double *pCurley = g->Curley;
  /*complex*/ double *pCurlez = g->Curlez;
  /*complex*/ double *pex = g->ex;
  /*complex*/ double *pey = g->ey;
  /*complex*/ double *pez = g->ez;

  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;
  complex double *pixsensor = g->ixsensor;
  complex double *pkref = g->kref;

  double pi = 3.141592654;
  if (ppw <= 0) {
    fprintf(stderr,
       "ezInc: ezIncInit() must be called before ezInc.\n"
       "       Points per wavelength must be positive.\n");
    exit(-1);
  }
  c0 = 3e8;
  ddx = DDz;
  Nlam = 3000.0;


  arg = (time*pi)/(3*Nlam) - 10.0;
  arg = arg*arg;
  Gauss = exp(-0.5*arg)*cos((time*pi/Nlam) - 30.0);


  // arg = ((time - 0.5)*pi)/(3*Nlam) - 10.0;
  // arg = arg*arg;
  // Gauss1 = exp(-0.5*arg)*cos(((time - 0.5)*pi/Nlam) - 30.0);

  pp = 10 + PMLs;

  // #pragma acc enter data copyin(pex[0:tc], pey[0:tc], pez[0:tc], pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc])
  // {
    #pragma acc parallel loop async(1) independent collapse(2) present(pCurlex[0:tc], pey[0:tc],pez[0:tc])
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        pCurlEx(mm, nn, pp) = pCurlEx(mm, nn, pp) + Gauss/DDz;
      }
    }

    // #pragma acc parallel loop async(2) independent present(pCurlex[0:tc], pey[0:tc],pez[0:tc])
    // for (mm = 0; mm < SizeX0; mm++){
    //  pCurlEx(mm, SizeY0 - 1, pp) = (pEz(mm, 0, pp) - pEz(mm, SizeY0 - 1, pp))/DDy
    //   - (pEy(mm, SizeY0 - 1, pp + 1) - pEy(mm, SizeY0 -1, pp))/DDz + Gauss/DDz;
    // }

    #pragma acc parallel loop async(3) independent collapse(2) present(pCurley[0:tc], pex[0:tc],pez[0:tc])
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        pCurlEy(mm, nn, pp) = pCurlEy(mm, nn, pp) - Gauss/DDz;
      }
    }  
    // #pragma acc parallel loop async(4) independent present(pCurley[0:tc], pex[0:tc],pez[0:tc])
    // for (nn = 0; nn < SizeY0; nn++){
    //     pCurlEy(SizeX0 - 1, nn, pp) = (pEx(SizeX0 - 1, nn, pp + 1) - pEx(SizeX0 - 1, nn, pp))/DDz
    //      - (pEz(0, nn, pp) - pEz(SizeX0 - 1, nn, pp))/DDx- Gauss/DDz;
    // }
    // #pragma acc exit data copyout(pex[0:tc], pey[0:tc], pez[0:tc], pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc])
  // }
  // #pragma acc parallel loop independent present(pixsensor[:NFREQS], pkref[:NFREQS])
  // for(pp = 0; pp < NFREQS; pp++){
  //   pIXSensor(pp) = pIXSensor(pp) + creal(pow(pKRef(pp), (time))*Gauss*(pow(pKRef(pp), (time))*Gauss))/2;
  // }

    #pragma acc wait

    return;
}


void TFSF_H(double time, Grid *g) {
  double arg;
  double rick, sig, dt, ddx, t0, tau, Gauss, c0, del_t, Hysr, Hxsr, Nlam;
  int mm, nn, pp;

  double imp0 = 377.0;
  int tc = SizeX0*SizeY0*SizeZ0;
  /*complex*/ double *phx = g->hx;
  /*complex*/ double *phy = g->hy;
  /*complex*/ double *phz = g->hz;
  
  /*complex*/ double *pCurlhx = g->Curlhx;
  /*complex*/ double *pCurlhy = g->Curlhy;

  complex double *pixsensor = g->ixsensor;
  complex double *pkref = g->kref;  

  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;

  double pi = 3.141592654;
  if (ppw <= 0) {
    fprintf(stderr,
       "ezInc: ezIncInit() must be called before ezInc.\n"
       "       Points per wavelength must be positive.\n");
    exit(-1);
  }
  c0 = 3e8;
  ddx = DDz;
  Nlam = 3000.0;

  arg = ((time - 0.5)*pi)/(3*Nlam) - 10.0;
  arg = arg*arg;
  Gauss = exp(-0.5*arg)*cos(((time - 0.5)*pi/Nlam) - 30.0);

  pp = 10 + PMLs;
  // Hxsr = -1*Gauss;
  // Hysr =    Gauss;
  //printf("0 CurlHx(20, 20, 30) is %g \n", CurlHx(20, 20, 30));

  // #pragma acc enter data copyin(pCurlhx[0:tc], pCurlhy[0:tc])
  // {
  #pragma acc parallel loop async(1) independent collapse(2) present(pCurlhx[0:tc],phz[0:tc], phy[0:tc])  
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
      pCurlHx(mm, nn, pp) = pCurlHx(mm, nn, pp) + Gauss/DDz;
      }
    }
  // #pragma acc parallel loop async(2) independent present(pCurlhx[0:tc],phz[0:tc], phy[0:tc])
  //   for (mm = 0; mm < SizeX0; mm++){
  //       pCurlHx(mm, 0, pp) = pCurlHx(mm, nn, pp) + Gauss/DDz;
  //   }

    /*for periodic boundary may need to add in the TFSF term explicitly*/
  #pragma acc parallel loop async(3) independent collapse(2) present(pCurlhy[0:tc],phx[0:tc], phz[0:tc])
    for (mm = 0; mm < SizeX0; mm++){
        for (nn = 0; nn < SizeY0; nn++){
          pCurlHy(mm, nn, pp) = pCurlHy(mm, nn, pp)  + Gauss/DDz;
        }
    }

  // #pragma acc parallel loop async(4) independent present(pCurlhy[0:tc],phx[0:tc], phz[0:tc])
  //   for (nn = 0; nn < SizeY0; nn++){
  //       pCurlHy(0, nn, pp) = pCurlHy(mm, nn, pp)  + Gauss/DDz;
  //   }
  // #pragma acc exit data copyout(phx[0:tc],phy[0:tc],phz[0:tc],pCurlhx[0:tc], pCurlhy[0:tc], pCurlhz[0:tc])
  #pragma acc parallel loop independent present(pixsensor[:NFREQS], pkref[:NFREQS])
  for(pp = 0; pp < NFREQS; pp++){
    pIXSensor(pp) = pIXSensor(pp) + cpow(pKRef(pp), (time))*Gauss;
  }
    //printf("CurlHx(20, 20, 30) is %g \n \n", CurlHx(20, 20, 30));
  #pragma acc wait
    return;
}
