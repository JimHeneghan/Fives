#include "fdtd-macro.h"
#include "fdtd-alloc.h"
#include "p-macro.h"
#include <complex.h>
#include <math.h>
#include <openacc.h>
// #include <cuda_runtime_api.h>
void gridInit(Grid *g) {
  double imp0 = 377.0, ddx;
  
  int mm, nn, pp, i;
  double XCenter1, XCenter2, YCenter1, YCenter2, r2, XLocC, YLocC, rad, XLoc1, XLoc2, YLoc1, YLoc2, XLoc3, YLoc3, XLoc4, YLoc4, XLoc5, YLoc5, dist1, dist2, dist3, dist4, dist5;
  double epsr, a, PML_no, eps0, c0, dt, df, freq;
  Type = threeDGrid;   /*@ \label{grid3dhomoA} @*/
  

  /* Defining structural constants from p-macro.h */
  PML_no = PMLs;
  MaxTime = MT;   
  SizeX = SizeX0; 
  SizeY = SizeY0;
  SizeZ = SizeZ0;

  a = SizeX;
  c0 = 3e8;
  Cdtds = 1.0/2; // Courant number /*@ \label{grid3dhomoB} @*/
  ddx = DDz0;
  dt = Cdtds*ddx/c0;
  eps0 = 8.85418782e-12;
  df = (9e-6)/NFREQS;
  rad = 0.68e-6/DDx0; 
  XCenter1 = 0;
  YCenter1 = 0;
  // XCenter2 = 65;
  // YCenter2 = 37;
  pp = 1;
  r2 = rad*rad;

  /* Bloch boundary constants */
  Phix = cexp(I*0*SizeX*ddx);
  Phiy = cexp(I*0*SizeY*ddx);
  Phiz = cexp(I*0*SizeY*ddx);

  printf("SizeX is %d \t SizeY is %d \t SizeZ is %d, dt is %g \n", SizeX0, SizeY0, SizeZ0, dt);
  /* memory allocationi*/
  ALLOC_3D(g->hx,   SizeX, SizeY, SizeZ, /*complex*/ double); 
  ALLOC_3D(g->hy,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->hz,   SizeX, SizeY, SizeZ, /*complex*/ double);

  ALLOC_3D(g->ex,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->ey,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->ez,   SizeX, SizeY, SizeZ, /*complex*/ double);

  ALLOC_3D(g->e1x,  SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->e1y,  SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->e1z,  SizeX, SizeY, SizeZ, /*complex*/ double);

  ALLOC_3D(g->dx,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->dy,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->dz,   SizeX, SizeY, SizeZ, /*complex*/ double);

  /*Curl arrays*/
  ALLOC_3D(g->Curlhx,   SizeX, SizeY, SizeZ, /*complex*/ double); 
  ALLOC_3D(g->Curlhy,   SizeX, SizeY, SizeZ, /*complex*/ double); 
  ALLOC_3D(g->Curlhz,   SizeX, SizeY, SizeZ, /*complex*/ double); 

  ALLOC_3D(g->Curlex,   SizeX, SizeY, SizeZ, /*complex*/ double); 
  ALLOC_3D(g->Curley,   SizeX, SizeY, SizeZ, /*complex*/ double); 
  ALLOC_3D(g->Curlez,   SizeX, SizeY, SizeZ, /*complex*/ double); 

  /* PML memory allocationi*/
  ALLOC_1D(g->sigmazD,    SizeZ, /*complex*/ double);
  ALLOC_1D(g->sigmazH,    SizeZ, /*complex*/ double);

  ALLOC_1D(g->PMLhx_0,   SizeZ, /*complex*/ double);
  ALLOC_1D(g->PMLhx_1,   SizeZ, /*complex*/ double);
  ALLOC_1D(g->PMLhx_2,   SizeZ, /*complex*/ double);

  ALLOC_1D(g->PMLhz_2,   SizeZ, /*complex*/ double);
  ALLOC_1D(g->PMLhz_3,   SizeZ, /*complex*/ double); 

  ALLOC_1D(g->PMLdx_0,   SizeZ, /*complex*/ double);
  ALLOC_1D(g->PMLdx_1,   SizeZ, /*complex*/ double);
  ALLOC_1D(g->PMLdx_2,   SizeZ, /*complex*/ double);
 
  ALLOC_1D(g->PMLdz_2,   SizeZ, /*complex*/ double);
  ALLOC_1D(g->PMLdz_3,   SizeZ, /*complex*/ double);

  /*PML integration arrays*/
  ALLOC_3D(g->IChz,   SizeX, SizeY, SizeZ, /*complex*/ double); 
  ALLOC_3D(g->ICez,   SizeX, SizeY, SizeZ, /*complex*/ double); 

  /*memory allocation for z transform fields for X comp*/  

  // ALLOC_3D(g->ixx,   SizeX, SizeY, SizeZ, /*complex*/ double);

  ALLOC_3D(g->sxx,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->s1xx,  SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->s2xx,  SizeX, SizeY, SizeZ, /*complex*/ double);

  /*memory allocation for z transform fields for Y comp*/
  // ALLOC_3D(g->iyy,   SizeX, SizeY, SizeZ, /*complex*/ double);

  ALLOC_3D(g->syy,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->s1yy,  SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->s2yy,  SizeX, SizeY, SizeZ, /*complex*/ double);

  /*memory allocation for z transform fields for Z comp*/
  // ALLOC_3D(g->iz,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->sz,   SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->s1z,  SizeX, SizeY, SizeZ, /*complex*/ double);
  ALLOC_3D(g->s2z,  SizeX, SizeY, SizeZ, /*complex*/ double);
  
  ALLOC_3D(g->epsR, SizeX, SizeY, SizeZ,  double);
  ALLOC_3D(g->media, SizeX, SizeY, SizeZ,  int);

  /* Sensor array memory allocation */
  ALLOC_1D(g->ixsensor,  NFREQS, complex double);

  ALLOC_3D(g->rxsensor,  SizeX, SizeY, NFREQS, complex double);
  // ALLOC_3D(g->rxhsensor, SizeX, SizeY, NFREQS, complex double);

  ALLOC_3D(g->rysensor,  SizeX, SizeY, NFREQS, complex double);
  ALLOC_3D(g->rzsensor,  SizeX, SizeY, NFREQS, complex double);

  // ALLOC_3D(g->ryhsensor, SizeX, SizeY, NFREQS, complex double);

  ALLOC_3D(g->txsensor,  SizeX, SizeY, NFREQS, complex double);
  ALLOC_3D(g->tysensor,  SizeX, SizeY, NFREQS, complex double);
  ALLOC_3D(g->tzsensor,  SizeX, SizeY, NFREQS, complex double);

  // ALLOC_3D(g->txhsensor, SizeX, SizeY, NFREQS, complex double);

  // ALLOC_3D(g->tyhsensor, SizeX, SizeY, NFREQS, complex double);

  // ALLOC_1D(g->rzsensor,   200, complex double);
  ALLOC_1D(g->rxtime,   MaxTime, double);
  ALLOC_1D(g->txtime,   MaxTime, double);
  // ALLOC_1D(g->iysensor,   MaxTime, /*complex*/ double);

  // ALLOC_3D(g->pz,  SizeX, SizeY, NFREQS, complex double);
  // ALLOC_3D(g->pzT, SizeX, SizeY, NFREQS, complex double);


  ALLOC_1D(g->kref,   NFREQS, complex double);
  
  ALLOC_3D(g->exMovie,   FRAMES, SizeX, SizeZ, double);
  ALLOC_3D(g->exMoviexy,   FRAMES, SizeX, SizeY, double);


  for (i = 0; i < NFREQS; i++){
    // KRef(i) = cexp(-1*I*2*M_PI*dt*(30e12 + i*540e9)); //150 GHz being 1/200th of 30 THz, 30 THz being the difference btn 5 and 10 microns
    freq = c0/(1e-6 + i*df);
    KRef(i) = cexp(-1*I*2*M_PI*dt*freq);

  }

  /*Media = 1 is magnetized plasma
    Media = 2 is non dispersive
    Media = 3 is reversed biased magnetized plasma*/
  /* Arrays are zeroed when memory is allocated to them*/
  /* Arraays that need to be initialized to other values are reinitialized here*/
  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
       for (pp = 0; pp < SizeZ; pp++) {
          EpsR(mm, nn, pp) = 1.0;
          Media(mm, nn, pp) = 2;
	}


  /* PML conductivity definition only done for 1 dimension as only applied in the Z direction*/
  /* SF PML Sigma*/
  double arg;
  int pmpm;
  // for (pp = 0; pp < PML_no + 1; pp++) {
  //   pmpm = (SizeZ - pp - 1);
  //   arg = (PML_no - pp)/(1.0*PML_no);
  //   SigmaZ(pmpm) = (0.33)*arg*arg*arg;
  //   SigmaZ(pp) = (0.33)*arg*arg*arg;
  // }
  /* TF PML Sigma*/
  int kd = 0;
  int kh = 0;
  for (pp = 0; pp < 2*PML_no; pp++) {
    // pmpm = (SizeZ - 2*PML_no - 1 + pp);
    arg = (2*PML_no - pp)/(2.0*PML_no);
    // SigmaZ(pmpm) = (0.33)*arg*arg*arg;
    // SigmaZ(pp) = (0.33)*arg*arg*arg;
    if (pp%2 == 0){
      pmpm = (SizeZ - kh);
      SigmaZD(kh) = (eps0/(2*dt))*arg*arg*arg;
      kh++;
    }else{
      pmpm = (SizeZ - kd);
      SigmaZH(kd) = (eps0/(2*dt))*arg*arg*arg;
      kd++;
    }
  }
 printf("\n \n **************************** \n \n");
  /* TF PML Sigma*/
  kd = 0;
  kh = 0;
  for (pp = 0; pp < 2*PML_no; pp++) {
    // pmpm = (SizeZ - 2*PML_no - 1 + pp);
    arg = (2*PML_no - pp)/(2.0*PML_no);
    // SigmaZ(pmpm) = (0.33)*arg*arg*arg;
    // SigmaZ(pp) = (0.33)*arg*arg*arg;
    if (pp%2 == 0){
      pmpm = (SizeZ - kh);
      SigmaZH(pmpm) = (eps0/(2*dt))*arg*arg*arg;
      kh++;
    }else{
      pmpm = (SizeZ - kd);
      SigmaZD(pmpm) = (eps0/(2*dt))*arg*arg*arg;
      kd++;
    }
  }
   printf("\n \n **************************** \n \n");
       for (pp = 0; pp < SizeZ; pp++) {
      printf("pp is %d, SigmaZ D is %f\n",pp, SigmaZD(pp));
      printf("pp is %d, SigmaZ H is %f\n",pp, SigmaZH(pp));
  }

  /*Dx PML update constants*/
  for (pp = 0; pp < SizeZ; pp++) {
    PMLDx_0(pp) = (1/dt) + (SigmaZD(pp)/(2*eps0));
    PMLDx_1(pp) = ((1/dt) - (SigmaZD(pp)/(2*eps0)))/PMLDx_0(pp);
    PMLDx_2(pp) = (c0/PMLDx_0(pp));
    
  }
  /*Dy PML update constants*/
  /*Dz PML update constants*/ 
  for (pp = 0; pp < SizeZ; pp++) {
    PMLDz_2(pp) = c0*dt;
    PMLDz_3(pp) = c0*dt*dt*SigmaZD(pp)/eps0;
  }


  // for (pp = 0; pp < PML_no; pp++) {
  //   pmpm = (SizeZ - pp - 1);
  //   arg = (PML_no - pp - 0.5)/(1.0*PML_no);
  //   // SigmaZ(pmpm) = (0.33)*arg*arg*arg;
  //   // SigmaZ(pp) = (0.33)*arg*arg*arg;

  //   SigmaZ(pmpm) = (eps0/(2*dt))*arg*arg*arg;
  //   SigmaZ(pp) = (eps0/(2*dt))*arg*arg*arg;
  // }


  // for (pp = 0; pp < SizeZ - 1; pp++) {
  //   SigmaZH(pp) = (SigmaZD(pp) + SigmaZD(pp - 1))/2;
  // }


  
  /* PML constants for the Z boundary only PML*/
  /*Hx PML update constants*/
 
  for (pp = 0; pp < SizeZ; pp++) {
    PMLHx_0(pp) =  (1/dt) + (SigmaZH(pp)/(2*eps0));

    PMLHx_1(pp) = ((1/dt) - (SigmaZH(pp)/(2*eps0)))/PMLHx_0(pp);
    // printf("PMLHx_1 is %g \n", PMLHx_1(pp));
    /*assuming relative permiability of 1*/
    PMLHx_2(pp) = -1*(c0/PMLHx_0(pp));
    // printf("pp is %d, SigmaZ is %g , PMLHx_0 is %g \n", pp, SigmaZ(pp), PMLHx_0(pp));
    // printf("PMLHx_2 is %g \n", PMLHx_2(pp));
  }
  /*Hy PML update constants*/
  /*Not defined as identical Hx terms are reused*/ 


  /*Hz PML update constants*/ 
  for (pp = 0; pp < SizeZ; pp++) {
    /*assuming relative perm of 1*/
    PMLHz_2(pp) = -1*c0*dt;
    PMLHz_3(pp) = -1*c0*dt*dt*SigmaZH(pp)/eps0;
  }
printf("Hz PML assigned \n");

// /* putting in the hBN layer*/
//   for (mm = 0; mm < SizeX; mm++)
//       for (nn = 0; nn < SizeY; nn++) 
//         for (pp = (PML_no + 60); pp <(PML_no + 100); pp++) {
//           Media(mm, nn, pp)  = 1;

//         }


/* putting in the Ag layer*/
  for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY; nn++) 
        for (pp = (PML_no + 60); pp <(PML_no + 70); pp++) {
          Media(mm, nn, pp)  = 3;

        }

/* putting in the vacuum holes */
  for (mm = 0; mm < SizeX; mm++) /*@ \label{grid3dhomoE} @*/
      for (nn = 0; nn < SizeY; nn++)
       for (pp = (PML_no + 60); pp < (PML_no + 70); pp++) {

          XLoc1 = XCenter1 - mm;
          XLoc2 = XCenter1 - mm;
          XLoc3 = XCenter1 + a/2 - mm;
          XLoc4 = XCenter1 + a - mm;
          XLoc5 = XCenter1 + a - mm;

          YLoc1 = YCenter1 - nn;
          YLoc2 = YCenter1 + a*sqrt(3) - nn;
          YLoc3 = YCenter1 + a*sqrt(3)/2 - nn;
          YLoc4 = YCenter1 - nn;
          YLoc5 = YCenter1 + a*sqrt(3) - nn;

          dist1 = (int)(sqrt(XLoc1*XLoc1 + YLoc1*YLoc1));
          dist2 = (int)(sqrt(XLoc2*XLoc2 + YLoc2*YLoc2));
          dist3 = (int)(sqrt(XLoc3*XLoc3 + YLoc3*YLoc3));
          dist4 = (int)(sqrt(XLoc4*XLoc4 + YLoc4*YLoc4));
          dist5 = (int)(sqrt(XLoc5*XLoc5 + YLoc5*YLoc5));
          if((dist1< rad)||(dist2 < rad)||(dist3 < rad)||(dist4 < rad)||(dist5 < rad)){
            Media(mm, nn, pp) = 2;
            EpsR(mm, nn, pp) = 1.0;
          }
        }

  for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY; nn++) 
        for (pp = (PML_no + 70); pp <SizeZ0; pp++) {
          Media(mm, nn, pp)  = 2;
          EpsR(mm, nn, pp) = 11.7;
        }

/*media visualization */
  float temp;
  char filename[100];
  FILE *out;
  sprintf(filename,"Media/Media.dat");
  out = fopen(filename,"w");
  mm = 20;
  //nn = 5;
  for (pp = 0; pp < SizeZ; pp++) {
    for (nn = 0; nn < SizeY; nn++){
      temp = (float)Media(mm,nn, pp); // store data as a float
      if(pp >= PML_no+70){
        temp = 4;}
      if(pp < PML_no){
          temp = 0;}
      if(pp> (SizeZ - PML_no)){
          temp = 0;
      }
      
      fprintf(out, "%d \t %f \n", pp, temp);
    }
  }

  fclose(out);

  sprintf(filename,"Media/XY.dat");
  out = fopen(filename,"w");
  pp = PML_no + 65;
  for (nn = 0; nn < SizeY; nn++) {
    for (mm=0; mm<SizeX; mm++) {
      temp = (float)Media(mm,nn, pp); // store data as a float
      fprintf(out, "%f \n", temp);
    }
  }
  fclose(out);

  return;
}  /* end gridInit() */
