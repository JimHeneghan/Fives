#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"
#include "p-macro.h"
#include <complex.h>
#include <math.h>
#include <openacc.h>

/*code to record a sensor of the data*/
void sensorInit(Grid *g){
	char filename[100];
	FILE *out;
	sprintf(filename, "Ref/TraPoynt.txt");
  out = fopen(filename, "w");
	fprintf(out, "Trans Sensor \n");
	fclose(out);

  // sprintf(filename, "Inc/Inc.txt");
  // out = fopen(filename, "w");
  // fprintf(out, "Inc Sensor\n");
  // fclose(out);

  sprintf(filename, "Ref/ExRef.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ex Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/EyRef.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ey Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/EzRef.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ez Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/ExTra.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ex Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/EyTra.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ey Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/EzTra.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ez Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/freq.txt");
  out = fopen(filename, "w");
  fprintf(out, "freq \n");
  fclose(out);

  sprintf(filename, "Ref/ExRefTime.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ex Time Reflection Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/ExTraTime.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ex Time Transmission Sensor \n");
  fclose(out);

    sprintf(filename, "Ref/ExInc.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ex Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/EyInc.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ey Sensor \n");
  fclose(out);

  sprintf(filename, "Ref/EzInc.txt");
  out = fopen(filename, "w");
  fprintf(out, "Ez Sensor \n");
  fclose(out);

return;
}

// void IncSensor(Grid *g){
//         int t;
//         /*complex*/ double *pixsensor = g->ixsensor;
//         /*complex*/ double *piysensor = g->iysensor;
//         /*complex*/ double *pizsensor = g->izsensor;
//         #pragma acc exit data copyout(pixsensor[0:MT], piysensor[0:MT], pizsensor[0:MT])
//         double tempX, tempY, tempZ, itempX, itempY, itempZ;
//         char filename[100];
//         FILE *out;
//         sprintf(filename, "Inc/Inc.txt");
//         out = fopen(filename, "a");
//         for(t = 0; t < MT; t++){              
//                 /* print the time stamp and the Ex field right before the QWS*/
//                 tempX = IXSensor(t);
//                 // tempY = creal(IYSensor(t));
//                 // tempZ = creal(IZSensor(t));

//                 // itempX = cimag(IXSensor(t));
//                 // itempY = cimag(IYSensor(t));
//                 // itempZ = cimag(IZSensor(t));
//                 fprintf(out, "%g \n", tempX);
//                 // fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", tempX, itempX, tempY, itempY, tempZ, itempZ);//, cimag(Hz(80, 34, 1)), creal(Hz(33, 43, 1)), cimag(Hz(33, 43, 1)), creal(Hz(83, 8, 1)), cimag(Hz(83, 8, 1)), creal(Hz(81, 44, 1)), cimag(Hz(81, 44, 1)), creal(Hz(36, 44, 1)), cimag(Hz(36, 44, 1)));
//         }
//         fclose(out);

// }


void RefSensor(Grid *g){
  printf("Starting FT R op 1 \n");
  int ff, mm, nn;
  complex double *prxsensor = g->rxsensor;
  complex double *prysensor = g->rysensor;
  complex double *przsensor = g->rzsensor;

  complex double *ptxsensor = g->txsensor;
  complex double *ptysensor = g->tysensor;
  complex double *ptzsensor = g->tzsensor;

  complex double *pixsensor = g->ixsensor;

  double *prxtime =  g->rxtime;
  double *ptxtime =  g->txtime;


  int tfs = SizeX0*SizeY0*NFREQS;

  printf("Starting FT R op 2 \n");
  #pragma acc exit data copyout(prxsensor[0:tfs], prysensor[0:tfs], przsensor[0:tfs], pixsensor[0:NFREQS], ptxsensor[0:tfs], ptysensor[0:tfs], ptzsensor[0:tfs], prxtime[:MT], ptxtime[:MT])// prxhsensor[0:tfs], pryhsensor[0:tfs],
  printf("Starting FT R op 3 \n");
  double IncX, arg, freq, sum, poynt, df, c0;// tempY, tempZ, itempX, itempY, itempZ;
  char filename[100];
  double dt = DDz0/(2*3e8);
  c0 = 3e8;
  df = (9e-6)/NFREQS;
  FILE *out;

  /*****************************************************************************/
  /*************************    Ref Ex    **************************************/

  sprintf(filename, "Ref/ExRef.txt");
  out = fopen(filename, "a");

  printf("Starting FT R Ex \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(RXSensor(mm, nn, ff)), cimag(RXSensor(mm, nn, ff)));
      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Ref Ey    **************************************/

  sprintf(filename, "Ref/EyRef.txt");
  out = fopen(filename, "a");

  printf("Starting FT R Ey \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(RYSensor(mm, nn, ff)), cimag(RYSensor(mm, nn, ff)));      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Ref Ez    **************************************/

  sprintf(filename, "Ref/EzRef.txt");
  out = fopen(filename, "a");

  printf("Starting FT R Ez \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(RZSensor(mm, nn, ff)), cimag(RZSensor(mm, nn, ff)));
      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Tra Ex    **************************************/

  sprintf(filename, "Ref/ExTra.txt");
  out = fopen(filename, "a");

  printf("Starting FT T Ex \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(TXSensor(mm, nn, ff)), cimag(TXSensor(mm, nn, ff)));
      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Tra Ey    **************************************/

  sprintf(filename, "Ref/EyTra.txt");
  out = fopen(filename, "a");

  printf("Starting FT T Ey \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(TYSensor(mm, nn, ff)), cimag(TYSensor(mm, nn, ff)));      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Tra Ez    **************************************/

  sprintf(filename, "Ref/EzTra.txt");
  out = fopen(filename, "a");

  printf("Starting FT T Ez \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(TZSensor(mm, nn, ff)), cimag(TZSensor(mm, nn, ff)));
      }
    }
  }
  fclose(out);
  /*****************************************************************************/
  /*************************    freq    ****************************************/

  sprintf(filename, "Ref/freq.txt");
  out = fopen(filename, "a");

  printf("Starting FT  freq \n");

  for(ff = 0; ff < NFREQS; ff++){
    arg  = c0/(1e-6 + ff*df);
    fprintf(out, "%g \n", arg);
  }

  fclose(out);

  /*****************************************************************************/
  /*************************    Inc Ex    **************************************/

  sprintf(filename, "Ref/ExInc.txt");
  out = fopen(filename, "a");

  printf("Starting FT R Ex \n");

  for(ff = 0; ff < NFREQS; ff++){
    // printf("Starting FT Ex is %g \n", IXSensor(ff));
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(IXSensor(ff)), cimag(IXSensor(ff)));
      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Inc Ey    **************************************/

  sprintf(filename, "Ref/EyInc.txt");
  out = fopen(filename, "a");

  printf("Starting FT R Ey \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", creal(IXSensor(ff)), cimag(IXSensor(ff)));      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Inc Ez    **************************************/

  // sprintf(filename, "Ref/EzRef.txt");
  sprintf(filename, "Ref/EzInc.txt");
  out = fopen(filename, "a");

  printf("Starting FT R Ez \n");

  for(ff = 0; ff < NFREQS; ff++){
    for(mm = 0; mm < SizeX0; mm++){
      for(nn  = 0; nn < SizeY0; nn++){
        fprintf(out, "%g \t %g \n", 0, 0);
      }
    }
  }
  fclose(out);

  /*****************************************************************************/
  /*************************    Time Reflection    *****************************/

  sprintf(filename, "Ref/ExRefTime.txt");
  out = fopen(filename, "a");

  printf("Starting Time Reflection \n");

  for(ff = 0; ff < MT; ff++){
    fprintf(out, "%g \n", RXTime(ff));
  }

  fclose(out);

  /*****************************************************************************/
  /*************************    Time Transmission    ***************************/

  sprintf(filename, "Ref/ExTraTime.txt");
  out = fopen(filename, "a");

  printf("Starting Time Reflection \n");

  for(ff = 0; ff < MT; ff++){
    fprintf(out, "%g \n", TXTime(ff));
  }

  fclose(out);


  /*****************************************************************************/
  /*************************    Tra Poynt **************************************/

  // sprintf(filename, "Ref/TraPoynt.txt");
  // out = fopen(filename, "a");

  // printf("Starting FT R op \n");

  // for(ff = 0; ff < NFREQS; ff++){
  //   sum = 0;
  //   freq = (30e12 + df*ff);
  //   // printf("sum is %g, freq is %g \n", sum, freq);
  //   for(mm = 0; mm < SizeX0; mm++){
  //     for(nn  = 0; nn < SizeY0; nn++){
  //       arg  = PzT(mm, nn, ff);//(cabs(Pz(mm, nn, ff)))/cabs(IncX);
  //       // arg  = cabs(arg*arg);
  //       sum = sum + arg;
  //     }
  //   }
  //     sum = sum;
  //     fprintf(out, "%g \t %g \n", freq, sum);
  // }
  // fclose(out);

return;      
}


// void Transmission(Grid *g){
//   int t, mm, nn, ff;
//   complex double *ptxsensor = g->txsensor;
//   // complex double *ptysensor = g->tysensor;
//   // /*complex*/ double *ptzsensor = g->tzsensor;

//   double TranX, IncX, tempZ, itempX, itempY, itempZ, arg, freq, sum;

//   complex double *pixsensor = g->ixsensor;
//   int tfs = SizeX0*SizeY0*NFREQS;
//   #pragma acc exit data copyout(ptxsensor[0:tfs], pixsensor[0:NFREQS])
//   char filename[100];
//   FILE *out;
//   sprintf(filename, "Trans/Trans.txt");
//   out = fopen(filename, "a");
//   for(ff = 0; ff < NFREQS; ff++){
//     sum = 0;
//     freq = (30e12 + 50e9*ff);
//     // printf("sum is %g, freq is %g \n", sum, freq);
//     for(mm = 0; mm < SizeX0; mm++){
//       for(nn = 0; nn < SizeY0; nn++){
//         TranX = creal(TXSensor(mm, nn, ff));
//         IncX = creal(IXSensor(ff));
//         arg  = TranX/IncX;
//         arg  = cabs(arg*arg);
//         // freq = (30e12 + 6e11*ff);//3e8/((500.0 + 5*ff)*DDz0);
//         sum = sum + arg;
//         // fprintf(out, "%d \t %d \t %g \n", mm, nn, RXSensor(mm, nn, ff));

//       }
//     }
//       sum = sum/((SizeX0)*SizeY0);
//       fprintf(out, "%g \t %g \t %g \n", freq, sum, TranX);
//   }
//   fclose(out);

// 	return;      
// }




