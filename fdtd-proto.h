#ifndef _FDTD_PROTO_H
#define _FDTD_PROTO_H
#include <complex.h>
#include "fdtd-grid.h"

/* Function prototypes */
void abcInit(Grid *g);
void abc(Grid *g);

void gridInit(Grid *g);

void snapshot3dInit(Grid *g);
void snapshot3d(Grid *g);

void updateE(Grid *g, int t);
void updateH(Grid *g, int t);
void updateD(Grid *g);
void updateP(Grid *g, int t);

void CompCurlE(Grid *g);
void CompCurlH(Grid *g);

void updateIntE(Grid *g);
void updateIntH(Grid *g);


void sensorInit(Grid *g);
void Transmission(Grid *g);
void IncSensor(Grid *g);
void RefSensor(Grid *g);

void Movie(Grid *g, int t);

#endif
