#include "fdtd-macro.h"
#include "p-macro.h"
#include <stdio.h>
#include <complex.h>
// #include <math.h>.
#include <accelmath.h>
#include <stdlib.h>
#include <openacc.h>


#define SQR(x) ((x)*(x))

/* update magnetic field */

void CompCurlE(Grid *g){
  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;
  // printf("started CurlE \n");
  /*Calculate the Hx field*/
  int mm, nn, pp;

  int tc = SizeX0*SizeY0*SizeZ0;
  /*complex*/ double *pCurlex = g->Curlex;
  /*complex*/ double *pCurley = g->Curley;
  /*complex*/ double *pCurlez = g->Curlez;

  /*complex*/ double *pex = g->ex;
  /*complex*/ double *pey = g->ey;
  /*complex*/ double *pez = g->ez;
  // #pragma acc parallel  
  // {
        // #pragma acc loop independent collapse(3)
  #pragma acc kernels loop independent collapse(3) present(pCurlex[0:tc], pey[0:tc],pez[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      // #pragma acc loop independent
      for (nn = 0; nn < SizeY0 - 1; nn++){
        // #pragma acc loop independent
        for (pp = 0; pp < SizeZ0 - 1; pp++){
          pCurlEx(mm, nn, pp) = (pEz(mm, nn + 1, pp) - pEz(mm, nn, pp))/DDy
            - (pEy(mm, nn, pp + 1) - pEy(mm, nn, pp))/DDz;
      }
     }
    }

  #pragma acc kernels loop independent collapse(2) present(pCurlex[0:tc], pey[0:tc],pez[0:tc])  
    for (mm = 0; mm < SizeX0; mm++){
      for (pp = 0; pp < SizeZ0 -1; pp++){
       pCurlEx(mm, SizeY0 - 1, pp) = (pEz(mm, 0, pp) - pEz(mm, SizeY0 - 1, pp))/DDy
        - (pEy(mm, SizeY0 - 1, pp + 1) - pEy(mm, SizeY0 -1, pp))/DDz;
      }
    }


    /* Calculate the Curl Ey field */
  #pragma acc kernels loop independent collapse(3) present(pCurley[0:tc], pex[0:tc],pez[0:tc]) 
    for (mm = 0; mm < SizeX0 - 1; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEy(mm, nn, pp) = (pEx(mm, nn, pp + 1) - pEx(mm, nn, pp))/DDz
          - (pEz(mm + 1, nn, pp) - pEz(mm, nn, pp))/DDx;
        }
      }
   }

  #pragma acc kernels loop independent collapse(2) present(pCurley[0:tc], pex[0:tc],pez[0:tc]) 
    for (nn = 0; nn < SizeY0; nn++){
       for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEy(SizeX0 - 1, nn, pp) = (pEx(SizeX0 - 1, nn, pp + 1) - pEx(SizeX0 - 1, nn, pp))/DDz
          - (pEz(0, nn, pp) - pEz(SizeX0 - 1, nn, pp))/DDx;
        }
      }

      /* Calculate the Curl of Ez field */ 
  #pragma acc kernels loop independent collapse(3) present(pCurlez[0:tc], pex[0:tc],pey[0:tc])     
    for (mm = 0; mm < SizeX0 - 1; mm++){
     for (nn = 0; nn < SizeY0 - 1; nn++){
       for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEz(mm, nn, pp) = (pEy(mm + 1, nn, pp) - pEy(mm, nn, pp))/DDx
          - (pEx(mm, nn + 1, pp) - pEx(mm, nn, pp))/DDy;
        }
      }
    }
    #pragma acc kernels loop independent collapse(2) present(pCurlez[0:tc], pex[0:tc],pey[0:tc])  
     for (nn = 0; nn < SizeY0 - 1; nn++){
       for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEz(SizeX0 - 1, nn, pp) = (pEy(0, nn, pp) - pEy(SizeX0 - 1, nn, pp))/DDx
          - (pEx(SizeX0 - 1, nn + 1, pp) - pEx(SizeX0 - 1, nn, pp))/DDy; 
        }
      }
  #pragma acc kernels loop independent collapse(2) present(pCurlez[0:tc], pex[0:tc],pey[0:tc]) 
   for (mm = 0; mm < SizeX0 - 1; mm++){
     for (pp = 0; pp < SizeZ0 -1; pp++){
       pCurlEz(mm, SizeY0  - 1, pp) = (pEy(mm + 1, SizeY0  - 1, pp) - pEy(mm, SizeY0  - 1, pp))/DDx
        - (pEx(mm, 0, pp) - pEx(mm, SizeY0  - 1, pp))/DDy;
     }
   }
  #pragma acc kernels loop independent present(pCurlez[0:tc], pex[0:tc],pey[0:tc]) 
   for (pp = 0; pp < SizeZ0 -1; pp++){
     pCurlEz(SizeX0 - 1, SizeY0 - 1, pp) = (pEy(0, SizeY0 - 1, pp) - pEy(SizeX0 - 1, SizeY0 - 1, pp))/DDx
      - (pEx(SizeX0 - 1, 0, pp) - pEx(SizeX0 - 1, SizeY0  - 1, pp))/DDy;
   }


  #pragma acc wait
  return;
} /*End Comp Curl E*/


void updateIntE(Grid *g){
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;
  /*complex*/ double *pCurlez = g->Curlez;
  /*complex*/ double *pCurlex = g->Curlex;
  /*complex*/ double *pCurley = g->Curley;
  /*complex*/ double *pICez = g->ICez;
  // #pragma acc enter data copyin(pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc]) 
  {
    #pragma acc kernels loop independent collapse(3) present(pICez[:tc], pCurlez[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
       for (nn = 0; nn < SizeY0; nn++){
         for (pp = 0; pp < SizeZ0 - 1; pp++){
          pICEz(mm, nn, pp) = pICEz(mm, nn, pp) + pCurlEz(mm, nn, pp);
         }
       }
    }
  }
   return;
}


void updateH(Grid *g, int t) {
  int mm, nn, pp,ff;
  int tc = SizeX0*SizeY0*SizeZ0;
  
  complex double *pkref = g->kref;
  complex double *prxhsensor = g->rxhsensor;
  complex double *pryhsensor = g->ryhsensor;

  int tfs = SizeX0*SizeY0*NFREQS;

  /*create new pointers that point to the struct pointers*/
  /*create new macro's for those pointers*/
  /*profit*/
  /*complex*/ double *pCurlex = g->Curlex;
  /*complex*/ double *pCurley = g->Curley;
  /*complex*/ double *pCurlez = g->Curlez;
  /*complex*/ double *pICez = g->ICez;

  /*complex*/ double *pPMLhx_1 = g->PMLhx_1;
  /*complex*/ double *pPMLhx_2 = g->PMLhx_2;

  /*complex*/ double *pPMLhz_2 = g->PMLhz_2;
  /*complex*/ double *pPMLhz_3 = g->PMLhz_3;
  /*complex*/ double *phx = g->hx;
  /*complex*/ double *phy = g->hy;
  /*complex*/ double *phz = g->hz;
    /*Calculate the Hx field*/
    #pragma acc parallel loop async independent collapse(3) present(pCurlex[0:tc], phx[0:tc], pPMLhx_1[0:SizeZ0], pPMLhx_2[0:SizeZ0])
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 0; pp < SizeZ0 - 1; pp++){
          pHx(mm, nn, pp) =  pHx(mm, nn, pp) * pPMLHx_1(pp)
            + pCurlEx(mm, nn, pp)*pPMLHx_2(pp);
        }
      }
    }

    /* Calculate the Hy field */    
  #pragma acc parallel loop async independent collapse(3) present(pCurley[0:tc],phy[0:tc],pPMLhx_1[0:SizeZ0],pPMLhx_2[0:SizeZ0]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 0; pp < SizeZ0 - 1; pp++){
         pHy(mm, nn, pp) =  pHy(mm, nn, pp) * pPMLHx_1(pp) 
          + pCurlEy(mm, nn, pp) * pPMLHx_2(pp);
        }
      }
    }

   /* Calculate the Hz field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlez[0:tc],phz[0:tc],pPMLhz_3[0:SizeZ0],pPMLhz_2[0:SizeZ0],pICez[:tc])     
   for (mm = 0; mm < SizeX0; mm++){
     for (nn = 0; nn < SizeY0; nn++){
       for (pp = 0; pp < SizeZ0 - 1; pp++){
         pHz(mm, nn, pp) = pHz(mm, nn, pp) 
          + pCurlEz(mm, nn, pp) * pPMLHz_2(pp)
            + pICEz(mm, nn, pp) * pPMLHz_3(pp); 
       }
     }
   }

  #pragma acc wait
   double arg;
  // #pragma acc parallel loop independent collapse(3) present(prxhsensor[:tfs], pkref[:NFREQS], phx[:tc])
  // for (mm = 0; mm < SizeX0 - 1; mm++){
  //   for (nn = 0; nn < SizeY0; nn++){  
  //     for(ff = 0; ff < NFREQS; ff++){
  //       arg = (pHx(mm, nn, 10 + PMLs) + pHx(mm + 1, nn, 10 + PMLs))/2;
  //       pRXHSensor(mm, nn, ff) = pRXHSensor(mm, nn, ff) + pow(pKRef(ff), (t+1))*arg;
  //     }
  //   }
  // }

  // #pragma acc parallel loop independent collapse(3) present(pryhsensor[:tfs], pkref[:NFREQS], phy[:tc])
  // for (mm = 0; mm < SizeX0; mm++){
  //   for (nn = 0; nn < SizeY0 - 1; nn++){  
  //     for(ff = 0; ff < NFREQS; ff++){
  //       arg = (pHy(mm, nn, 10 + PMLs) + pHy(mm, nn + 1, 10 + PMLs))/2;
  //       pRYHSensor(mm, nn, ff) = pRYHSensor(mm, nn, ff) + pow(pKRef(ff), (t+1))*arg;
  //     }
  //   }
  // }
  return;
}  /* end updateH() */


/* Calculate Curl Hx field */
void CompCurlH(Grid *g){
    
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;

  /*complex*/ double *pCurlhx = g->Curlhx;
  /*complex*/ double *pCurlhy = g->Curlhy;
  /*complex*/ double *pCurlhz = g->Curlhz;

  /*complex*/ double *phx = g->hx;
  /*complex*/ double *phy = g->hy;
  /*complex*/ double *phz = g->hz;
  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;

  #pragma acc kernels loop independent collapse(3) present(pCurlhx[0:tc], phy[0:tc],phz[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 1; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pCurlHx(mm, nn, pp) = (pHz(mm, nn, pp) - pHz(mm, nn - 1, pp))/DDy -
            (pHy(mm, nn, pp) - pHy(mm, nn, pp - 1))/DDz;
        }
      }
    }

  #pragma acc kernels loop independent collapse(2) present(pCurlhx[0:tc], phy[0:tc],phz[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (pp = 1; pp < SizeZ0; pp++){
        pCurlHx(mm, 0, pp) = (pHz(mm, 0, pp) - pHz(mm, SizeY0 - 1, pp))/DDy -
         (pHy(mm, 0, pp) - pHy(mm, 0, pp - 1))/DDz;
      }
    }

  /* Calculate Curl Hy field */
  #pragma acc kernels loop independent collapse(3) present(pCurlhy[0:tc], phx[0:tc],phz[0:tc])     
    for (mm = 1; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pCurlHy(mm, nn, pp) = (pHx(mm, nn, pp) - pHx(mm, nn, pp - 1))/DDz -
            (pHz(mm, nn, pp) - pHz(mm - 1, nn, pp))/DDx;
        }
      }
    }


  #pragma acc kernels loop independent collapse(2) present(pCurlhy[0:tc], phx[0:tc],phz[0:tc]) 
    for (nn = 0; nn < SizeY0; nn++){
      for (pp = 1; pp < SizeZ0; pp++){
        pCurlHy(0, nn, pp) = (pHx(0, nn, pp) - pHx(0, nn, pp - 1))/DDz -
          (pHz(0, nn, pp) - pHz(SizeX0 - 1, nn, pp))/DDx;
      }
    }

  /* Calculate Curl Hz field */
  #pragma acc kernels loop independent collapse(3) present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (mm = 1; mm < SizeX0; mm++){
      for (nn = 1; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pCurlHz(mm, nn, pp) = (pHy(mm, nn, pp) - pHy(mm - 1, nn, pp))/DDx -
            (pHx(mm, nn, pp) - pHx(mm, nn - 1, pp))/DDy;
        }
      }
    }

  #pragma acc kernels loop independent collapse(2) present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (nn = 1; nn < SizeY0; nn++){
      for (pp = 1; pp < SizeZ0; pp++){
        pCurlHz(0, nn, pp) = (pHy(0, nn, pp) - pHy(SizeX0 - 1, nn, pp))/DDx -
          (pHx(0, nn, pp) - pHx(0, nn - 1, pp))/DDy;
      }
    }

  #pragma acc kernels loop independent collapse(2) present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (mm = 1; mm < SizeX0; mm++){
      for (pp = 1; pp < SizeZ0; pp++){ 
        pCurlHz(mm, 0, pp) = (pHy(mm, 0, pp) - pHy(mm - 1, 0, pp))/DDx -
          (pHx(mm, 0, pp) - pHx(mm, SizeY0 - 1, pp))/DDy;
      }
    }

  #pragma acc kernels loop independent present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (pp = 1; pp < SizeZ0; pp++){
       pCurlHz(0, 0, pp) = (pHy(0, 0, pp) - pHy(SizeX0 - 1, 0, pp))/DDx - 
        (pHx(0, 0, pp) - pHx(0, SizeY0 - 1, pp))/DDy;
    }

  #pragma acc wait
  return;
} /*end CompCurlH()*/


/*update H integrations*/
void updateIntH(Grid *g){
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;

  /*complex*/ double *pCurlhx = g->Curlhx;
  /*complex*/ double *pCurlhy = g->Curlhy;
  /*complex*/ double *pCurlhz = g->Curlhz;

  /*complex*/ double *pIChz = g->IChz;

    #pragma acc kernels loop independent collapse(3) present(pIChz[:tc], pCurlhz[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn =0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pICHz(mm, nn, pp) = pICHz(mm, nn, pp) + pCurlHz(mm, nn, pp);
        }
      }
    }

  return;
}


/* update D field */
void updateD(Grid *g) {
  int mm, nn, pp, ff;
  /*complex*/ double CurlH;
  /*complex*/ double Psix = 1;//conj(Phix);
  /*complex*/ double Psiy = 1;//conj(Phiy);
  /*complex*/ double Psiz = 1;//conj(Phiz);

  int tc = SizeX0*SizeY0*SizeZ0;

  /*complex*/ double *pCurlhx = g->Curlhx;
  /*complex*/ double *pCurlhy = g->Curlhy;
  /*complex*/ double *pCurlhz = g->Curlhz;

  /*complex*/ double *pdx = g->dx;
  /*complex*/ double *pdy = g->dy;
  /*complex*/ double *pdz = g->dz;

  /*complex*/ double *pPMLdx_1 = g->PMLdx_1;
  /*complex*/ double *pPMLdx_2 = g->PMLdx_2;

  /*complex*/ double *pPMLdz_2 = g->PMLdz_2;
  /*complex*/ double *pPMLdz_3 = g->PMLdz_3;

  /*complex*/ double *pIChz = g->IChz;

    /* Calculate Dx field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlhx[0:tc], pdx[0:tc], pPMLdx_1[0:SizeZ0], pPMLdx_2[0:SizeZ0]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pDx(mm, nn, pp) = pDx(mm, nn, pp) * pPMLDx_1(pp) 
             + pCurlHx(mm, nn, pp) * pPMLDx_2(pp);
        }
      }
    }

  /* Calculate Dy field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlhy[0:tc], pdy[0:tc], pPMLdx_1[0:SizeZ0], pPMLdx_2[0:SizeZ0]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pDy(mm, nn, pp) = pDy(mm, nn, pp) * pPMLDx_1(pp) 
             + pCurlHy(mm, nn, pp) * pPMLDx_2(pp);
        }
      }
    }

  /* Calculate Dz field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlhz[0:tc], pdz[0:tc], pPMLdz_3[0:SizeZ0], pPMLdz_2[0:SizeZ0], pIChz[:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pDz(mm, nn, pp) = pDz(mm, nn, pp) 
            + pCurlHz(mm, nn, pp) * pPMLDz_2(pp)
              + pPMLDz_3(pp) * pICHz(mm, nn, pp);
        }
      }
    }

  #pragma acc wait

  return;
}  /* end updateD() */


void updateE(Grid *g, int t) {

  int mm, nn, pp, ff;
  /* Declaring Lorentz Equation Constants for hBN*/
  double w_nuXY, S_nu_XY, gam,  eps_inf,   A,  B,  C,  EXP1,  EXP2,  CS,  SN, C1XY, C2XY;
  double w_nuZ,  S_nu_Z,  gamZ, eps_infz,  Az, Bz, Cz, EXP1Z, EXP2Z, CSZ, SNZ, C1Z, C2Z;
  /* Declaring Drude Equation Constants for Ag*/
  double wp, gamAg, EXPAg, eps_inf_Ag, AAg, C1Ag, C2Ag;

  /* Declaring more generic constants*/
  eps_inf = 4.87;
  eps_infz = 2.95;
  double epsilon0 = 8.85418782e-12;
  double dt = DDz0/(2*3e8);
  double Epsr = 1.0;
  /* Redeclaring constants for use in update equations to not be part of a struct */
  int MaxTime0 = MaxTime;
  int PMLs0 = PMLs;
  int t0 = t;
  int tc = SizeX0*SizeY0*SizeZ0;
  int tfs = SizeX0*SizeY0*NFREQS;

  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;
  int PML_no = PMLs;

  /* Declaring array pointers to pointers to Struct:Grid type pointes*/
  /*complex*/ double *pex = g->ex;
  /*complex*/ double *pey = g->ey;
  /*complex*/ double *pez = g->ez;

  /*complex*/ double *pdx = g->dx; 
  /*complex*/ double *pdy = g->dy;
  /*complex*/ double *pdz = g->dz;

  /*complex*/ double *psxx = g->sxx;
  /*complex*/ double *psyy = g->syy;
  /*complex*/ double  *psz = g->sz;

  /*complex*/ double *ps1xx = g->s1xx;
  /*complex*/ double *ps1yy = g->s1yy;
  /*complex*/ double  *ps1z = g->s1z; 

  double *prxtime =  g->rxtime;
  double *ptxtime =  g->txtime;


  // /*complex*/ double *piysensor = g->iysensor;
  // /*complex*/ double *ptysensor = g->tysensor;

  // /*complex*/ double *przsensor = g->rzsensor;
  // /*complex*/ double *pizsensor = g->izsensor;
  // /*complex*/ double *ptzsensor = g->tzsensor;

  int *pmedia = g->media;
  double *pepsR   =  g->epsR;
  /*Calculate Ex Field*/
  #pragma acc parallel loop async independent collapse(3) present(pex[:tc], pdx[:tc], psxx[:tc], ps1xx[:tc], pmedia[:tc], pepsR[:tc], prxtime[:MT], ptxtime[:MT]) 
   for (mm = 0; mm < SizeX0; mm++){
     for (nn = 0; nn < SizeY0; nn++){
       for (pp = 1; pp < SizeZ0; pp++){
        /* If statements to indicate which update equation to use
        /* Media = 1 is Lorentz model hBN
        /* Media = 2 is non dispersive
        /* Media = 3 is Drude Model Ag */
        if (pMedia(mm, nn, pp) == 1){
          pEx(mm, nn, pp)  = (pDx(mm, nn, pp) - pSxx(mm, nn, pp))/eps_inf;
        } else if(pMedia(mm, nn, pp) == 2) {
          pEx(mm, nn, pp)  = pDx(mm, nn, pp)/pEpsR(mm, nn, pp);
        } else if(pMedia(mm, nn, pp) == 3) {
          pEx(mm, nn, pp)  = (pDx(mm, nn, pp) - pSxx(mm, nn, pp));
        }
       }
      }
      pRXTime(t) = (pEx(25, 25, PMLs + 5));
      pTXTime(t) = (pEx(25, 25, PMLs + 115));
    }
      /*Incident, REflection and Transmission sensors*/
      /*These could be modified to FFt sensors       */
      /*or stand aloneor FFt functions in future     */


  /*Calculate Ey Field*/ 
  #pragma acc parallel loop async independent collapse(3) present(pey[:tc], pdy[:tc], psyy[:tc], ps1yy[:tc], pmedia[:tc],pepsR[:tc])//, prysensor[:MT], piysensor[:MT], ptysensor[:MT]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          if (pMedia(mm, nn, pp) == 1){
            pEy(mm, nn, pp)  = (pDy(mm, nn, pp) -  pSyy(mm, nn, pp))/eps_inf;
          } else if(pMedia(mm, nn, pp) == 2) { 
            pEy(mm, nn, pp)  = pDy(mm, nn, pp)/pEpsR(mm, nn, pp);
          } else if(pMedia(mm, nn, pp) == 3) {
            pEy(mm, nn, pp)  = (pDy(mm, nn, pp) - pSyy(mm, nn, pp));// + pS1yy(mm, nn, pp))/eps_inf_Ag;
          }
        }
      }
    }

  /*Calculate Ez Field*/
  #pragma acc parallel loop async independent collapse(3) present(pez[:tc], pdz[:tc], psz[:tc], ps1z[:tc], pmedia[:tc], pepsR[:tc])//, przsensor[:MT], pizsensor[:MT], ptzsensor[:MT]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          if (pMedia(mm, nn, pp) == 1){
            pEz(mm, nn, pp) = (pDz(mm, nn, pp) - pSz(mm, nn, pp))/eps_infz;
          } else if(pMedia(mm, nn, pp) == 2) {
            pEz(mm, nn, pp)  = pDz(mm, nn, pp)/pEpsR(mm, nn, pp);//(mm, nn, pp);
          } else if(pMedia(mm, nn, pp) == 3) {
            pEz(mm, nn, pp)  = (pDz(mm, nn, pp) - pSz(mm, nn, pp));
        }
      }
    }
  }
  #pragma acc wait

  /*Lorentz term update equations*/
    /*could be done in the E* update loops but done here to increase 
      portability to an off diagnal medium*/
    /*XY hBN Lorentz constants*/  
    w_nuXY = 2.58276e14;
    gam = 6.02526055e12;
    S_nu_XY = 1.83;

    A = (gam/2);
    B = w_nuXY*(sqrt(1 - SQR(gam/(2*w_nuXY))));
    C = w_nuXY/(sqrt(1 - SQR(gam/(2*w_nuXY))));
    EXP1 = exp(-A*dt);
    EXP2 = exp(-2*A*dt);
    CS = cos(B*dt);
    SN = sin(B*dt);

    /*Z hBN Lorentz constants*/ 
    w_nuZ = 1.41e14;
    gamZ = 6.02526055e12;//2.3864e12;
    S_nu_Z= 0.61;

    Az = (gamZ/2);
    Bz = w_nuZ*(sqrt(1 - SQR(gamZ/(2*w_nuZ))));
    Cz = w_nuZ/(sqrt(1 - SQR(gamZ/(2*w_nuZ))));
    EXP1Z = exp(-Az*dt);
    EXP2Z = exp(-2*Az*dt);
    CSZ = cos(Bz*dt);
    SNZ = sin(Bz*dt);
    
   double *ps2xx = g->s2xx;
   double *ps2yy = g->s2yy;
   double  *ps2z = g->s2z;

   double *pe1x = g->e1x;
   double *pe1y = g->e1y;
   double *pe1z = g->e1z;   

  #pragma acc kernels loop independent collapse(3) present(pex[:tc], pey[:tc], pez[:tc], pe1x[:tc], pe1y[:tc], pe1z[:tc], psxx[:tc], psyy[:tc], psz[:tc], ps1xx[:tc], ps1yy[:tc], ps1z[:tc], ps2xx[:tc], ps2yy[:tc], ps2z[:tc], pmedia[:tc]) 
  for (mm = 0; mm < SizeX0; mm++){
     for (nn = 0; nn < SizeY0; nn++){
       for (pp = 1; pp < SizeZ0; pp++){
          if (pMedia(mm, nn, pp) == 1){
            
            /* epsilon xx term */
            pSxx(mm, nn, pp) = 2*EXP1*CS*pS1xx(mm, nn, pp) - EXP2*pS2xx(mm, nn, pp) + EXP1*SN*C*dt*S_nu_XY*pE1x(mm, nn, pp);

            // /* epsilon yy term */
            pSyy(mm, nn, pp) = 2*EXP1*CS*pS1yy(mm, nn, pp) - EXP2*pS2yy(mm, nn, pp) + EXP1*SN*C*dt*S_nu_XY*pE1y(mm, nn, pp);

            /* epsilon zz term */
            pSz(mm, nn, pp) = 2*EXP1Z*CSZ*pS1z(mm, nn, pp) - EXP2Z*pS2z(mm, nn, pp) + EXP1Z*SNZ*Cz*dt*S_nu_Z*pE1z(mm, nn, pp);

            /*Updating the recent momory fields*/
            pS2xx(mm, nn, pp) = pS1xx(mm, nn, pp);
            pS1xx(mm, nn, pp) = pSxx(mm, nn, pp);

            pS2yy(mm, nn, pp) = pS1yy(mm, nn, pp);
            pS1yy(mm, nn, pp) = pSyy(mm, nn, pp);

            pS2z(mm, nn, pp) = pS1z(mm, nn, pp);
            pS1z(mm, nn, pp) = pSz(mm, nn, pp);

            pE1x(mm, nn, pp) = pEx(mm, nn, pp);
            pE1y(mm, nn, pp) = pEy(mm, nn, pp);
            pE1z(mm, nn, pp) = pEz(mm, nn, pp);
          }
        }
      }
    }

    /*Drude model for Ag from Palik */

    // /*complex*/ double *ps2xx = g->s2xx;
    // /*complex*/ double *ps2yy = g->s2yy;
    // /*complex*/ double  *ps2z = g->s2z;

    // /*complex*/ double *pe1x = g->e1x;
    // /*complex*/ double *pe1y = g->e1y;
    // /*complex*/ double *pe1z = g->e1z;   

    gamAg = 9.79125662e13;
    wp = 1.15136316e16;
    EXPAg =exp(-gamAg*dt);
    AAg = SQR(wp)*dt/gamAg;
    /* Z Transfrom update arrays */
    #pragma acc kernels loop independent collapse(3) present(pex[:tc], pey[:tc], pez[:tc], psxx[:tc], psyy[:tc], psz[:tc], ps1xx[:tc], ps1yy[:tc], ps1z[:tc], ps2xx[:tc], ps2yy[:tc], ps2z[:tc], pmedia[:tc]) //pez[:tc], pe1x[:tc], pe1y[:tc], pe1z[:tc],
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp <SizeZ0; pp++) {
          if (pMedia(mm, nn, pp) == 3){

            pSxx(mm, nn, pp) = (1 + EXPAg)*pS1xx(mm, nn, pp) - EXPAg*pS2xx(mm, nn, pp) + AAg*(1 - EXPAg)*pEx(mm, nn, pp);

            pSyy(mm, nn, pp) = (1 + EXPAg)*pS1yy(mm, nn, pp) - EXPAg*pS2yy(mm, nn, pp) + AAg*(1 - EXPAg)*pEy(mm, nn, pp);

            pSz(mm, nn, pp) = (1 + EXPAg)*pS1z(mm, nn, pp) - EXPAg*pS2z(mm, nn, pp) + AAg*(1 - EXPAg)*pEz(mm, nn, pp);

            pS2xx(mm, nn, pp) = pS1xx(mm, nn, pp);
            pS1xx(mm, nn, pp) = pSxx(mm, nn, pp);

            pS2yy(mm, nn, pp) = pS1yy(mm, nn, pp);
            pS1yy(mm, nn, pp) = pSyy(mm, nn, pp);

            pS2z(mm, nn, pp) = pS1z(mm, nn, pp);
            pS1z(mm, nn, pp) = pSz(mm, nn, pp);

          }
        }
      }
    }
return;
}

void updateP(Grid *g, int t) {

  complex double RefEx, RefEy, RefEz, RefHx, RefHy, TraEx, TraEy, TraEz, TraHx, TraHy;

  complex double *pkref = g->kref;

  complex double *prxsensor  = g->rxsensor;
  complex double *prysensor  = g->rysensor;

  complex double *prxhsensor = g->rxhsensor;
  complex double *pryhsensor = g->ryhsensor;

  complex double *przsensor  = g->rzsensor;


  complex double *ptxsensor  = g->txsensor;
  complex double *ptysensor  = g->tysensor;
  complex double *ptzsensor  = g->tzsensor;


  complex double *ptxhsensor = g->txhsensor;
  complex double *ptyhsensor = g->tyhsensor;
  // complex double *ptyhsensor = g->tyhsensor;


  complex double *ppz        = g->pz;
  complex double *ppzT       = g->pzT;

  double *pex = g->ex;
  double *pey = g->ey;
  double *pez = g->ez;

  double *phx = g->hx;
  double *phy = g->hy;
  int tfs = SizeX0*SizeY0*NFREQS;
  int tc = SizeX0*SizeY0*SizeZ0;
  int mm, nn, ff;


  #pragma acc parallel loop independent collapse(3) present(prxsensor[:tfs], prysensor[:tfs], przsensor[:tfs], ptxsensor[:tfs], ptysensor[:tfs], ptzsensor[:tfs], pkref[:NFREQS], pex[:tc], pey[:tc], pez[:tc])
    for(mm = 0; mm < SizeX0; mm++){
      for(nn = 0; nn < SizeY0; nn++){
        for(ff = 0; ff < NFREQS; ff++){

          // printf("start P %d, %d, %d\n", mm, nn, ff);

          RefEx =  pEx(mm, nn, 5 + PMLs);
          RefEy =  pEy(mm, nn, 5 + PMLs);
          RefEz =  pEz(mm, nn, 5 + PMLs);
          // RefHx = (pHx(mm, nn, 5 + PMLs) + pHx(mm, nn, 4 + PMLs))/2; 
          // RefHy = (pHy(mm, nn, 5 + PMLs) + pHy(mm, nn, 4 + PMLs))/2; 

          TraEx =  pEx(mm, nn, 105 + PMLs);
          TraEy =  pEy(mm, nn, 105 + PMLs);
          TraEz =  pEz(mm, nn, 105 + PMLs);

          // TraHx = (pHx(mm, nn, 115 + PMLs) + pHx(mm, nn, 114 + PMLs))/2; 
          // TraHy = (pHy(mm, nn, 115 + PMLs) + pHy(mm, nn, 114 + PMLs))/2;

          pRXSensor(mm, nn, ff)  = pRXSensor(mm, nn, ff)  + cpow(pKRef(ff), (t))*RefEx;
          pRYSensor(mm, nn, ff)  = pRYSensor(mm, nn, ff)  + cpow(pKRef(ff), (t))*RefEy;
          pRZSensor(mm, nn, ff)  = pRZSensor(mm, nn, ff)  + cpow(pKRef(ff), (t))*RefEz;

          // pRXHSensor(mm, nn, ff) = pRXHSensor(mm, nn, ff) + cpow(pKRef(ff), (t))*RefHx;

          // pRYHSensor(mm, nn, ff) = pRYHSensor(mm, nn, ff) + cpow(pKRef(ff), (t))*RefHy;

          pTXSensor(mm, nn, ff)  = pTXSensor(mm, nn, ff)  + cpow(pKRef(ff), (t))*TraEx;
          pTYSensor(mm, nn, ff)  = pTYSensor(mm, nn, ff)  + cpow(pKRef(ff), (t))*TraEy;
          pTZSensor(mm, nn, ff)  = pTZSensor(mm, nn, ff)  + cpow(pKRef(ff), (t))*TraEz;

          // pTXHSensor(mm, nn, ff) = pTXHSensor(mm, nn, ff) + cpow(pKRef(ff), (t))*TraHx;

          // pTYHSensor(mm, nn, ff) = pTYHSensor(mm, nn, ff) + cpow(pKRef(ff), (t))*TraHy;

          // pPz(mm, nn, ff)  = pPz(mm, nn, ff) 
          //                    + creal(pRXSensor(mm, nn, ff)*conj(pRYHSensor(mm, nn, ff))
          //                         - pRYSensor(mm, nn, ff)*conj(pRXHSensor(mm, nn, ff)));

          // pPzT(mm, nn, ff) = pPzT(mm, nn, ff) 
          //                    + creal(pTXSensor(mm, nn, ff)*conj(pTYHSensor(mm, nn, ff))
          //                          - pTYSensor(mm, nn, ff)*conj(pTXHSensor(mm, nn, ff)));
          
        }
      }
    }

  return;
}

void Movie(Grid *g, int t) {
  int mm, nn, pp;
  int tcm = FRAMES*SizeX0*SizeZ0;
  int tcm2 = FRAMES*SizeX0*SizeY0;
  int tc = SizeX0*SizeY0*SizeZ0;

  double *pex = g->ex;
  double *pexMoviexy = g->exMoviexy;
  double *pexMovie   = g->exMovie;


  #pragma acc kernels loop independent collapse(2) present(pex[:tc], pexMovie[:tcm])
  for (mm = 0; mm < SizeX0; mm++){
    for (pp = 0; pp < SizeZ0; pp++){
      // if(t%10 == 0){
      
      pExMovie(t, mm, pp) = creal(pEx(mm, 20, pp));      
      // }
    }
  }

  #pragma acc kernels loop independent collapse(2) present(pex[:tc], pexMoviexy[:tcm2])
  for (mm = 0; mm < SizeX0; mm++){
    for (nn = 0; nn < SizeY0; nn++){
      // if(t%10 == 0){
      
      pExMovieXY(t, mm, nn) = creal(pEx(mm, nn, PMLs + 80));      
      // }
    }
  }

  return;
}