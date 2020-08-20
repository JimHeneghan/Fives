#ifndef _EZINC_H
#define _EZINC_H

// #include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"
#include <complex.h>

/*Wiggly Gausian with phased sources*/
// void ezIncInit(Grid *g); 
//void ezInc(Grid *g, int time, double LocX, double LocY, double LocZ, double Kx, double Ky, double Kz);

/*Changing to a simpler source function*/
void ezIncInit(Grid *g); 
void TFSF_H(double time, Grid *g);
void TFSF_E(double time, Grid *g);


#endif
