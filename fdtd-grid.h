#ifndef _FDTD_GRID_H
#define _FDTD_GRID_H
// #include </*complex*/.h>
enum GRIDTYPE {oneDGrid, teZGrid, tmZGrid, threeDGrid};

struct Grid {
  /*complex*/ double *hx, *hy, *hz;
  double *chxh, *chxe;
  double *chyh, *chye;
  double *chzh, *chze;
  /*complex*/ double *ex, *ey, *ez;
  /*complex*/ double *e1x, *e1y, *e1z;
  double *cexe, *cexh;
  double *ceye, *ceyh;
  double *ceze, *cezh;
  int *media;

  /*complex*/ double *dx, *ixx, *ixy, *sxx, *s1xx, *s2xx, *sxy, *s1xy, *s2xy, *qxx, *q1xx, *q2xx, *qxy, *q1xy, *q2xy;
  /*complex*/ double *dy, *iyy, *iyx, *syy, *s1yy, *s2yy, *syx, *s1yx, *s2yx, *qyy, *q1yy, *q2yy, *qyx, *q1yx, *q2yx;
  /*complex*/ double *dz, *iz, *sz, *s1z, *s2z;

/*PML pointers*/
  /*complex*/ double *PMLhx_0, *PMLhy_0, *PMLhz_0;
  /*complex*/ double *PMLhx_1, *PMLhy_1, *PMLhz_1;
  /*complex*/ double *PMLhx_2, *PMLhy_2, *PMLhz_2;
  /*complex*/ double *PMLhx_3, *PMLhy_3, *PMLhz_3;
  /*complex*/ double *PMLhx_4, *PMLhy_4, *PMLhz_4;

  /*complex*/ double *PMLdx_0, *PMLdy_0, *PMLdz_0;
  /*complex*/ double *PMLdx_1, *PMLdy_1, *PMLdz_1;
  /*complex*/ double *PMLdx_2, *PMLdy_2, *PMLdz_2;
  /*complex*/ double *PMLdx_3, *PMLdy_3, *PMLdz_3;
  /*complex*/ double *PMLdx_4, *PMLdy_4, *PMLdz_4;

  /*complex*/ double *Curlhx, *Curlhy, *Curlhz;
  /*complex*/ double *Curlex, *Curley, *Curlez;

  /*complex*/ double *IChx, *IChy, *IChz;
  /*complex*/ double *ICex, *ICey, *ICez;

  /*complex*/ double *sigmax, *sigmay, *sigmazD, *sigmazH;

  complex double *rxsensor, *rxhsensor, *ixsensor, *txsensor, *txhsensor;
  complex double *rysensor, *ryhsensor, *iysensor, *tysensor, *tyhsensor;
  complex double *rzsensor, *izsensor, *tzsensor;

  double *rxtime, *txtime;

  complex double *kref, *pz, *pzT;

  /*Movie pointer*/
  double *exMovie, *exMoviexy;

  double *epsR, *sigma, *omega_p, *omega_b, *nu_c;
  int sizeX, sizeY, sizeZ;
  int time, maxTime;
  int type;
  double cdtds;
  /*complex*/ double phix, phiy, phiz;
};

typedef struct Grid Grid;

#endif

