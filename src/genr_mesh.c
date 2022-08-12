/*******************************************************************************
* genr_mesh.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [GPLv3 license]

*******************************************************************************/

#include "define.h"
#include "genr_mesh.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifdef OMP
#include <omp.h>
#define OMP_ATOMIC      _Pragma("omp atomic")
#else
#define OMP_ATOMIC
#endif

const char *powspec_assign_names[] = {
  "NGP",
  "CIC",
  "TSC",
  "PCS"
};


/*============================================================================*\
                 Implementation of particle assignment schemes
\*============================================================================*/

/******************************************************************************
Function `ngp`:
  Assign particles to the mesh with the Nearest-Grid-Point scheme.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the sample;
  * `Ng`:       number of grid cells per box side;
  * `bmin`:     lower boundaries of the box;
  * `Lbox`:     side length of the box;
  * `rho`:      the mesh.
******************************************************************************/
static void ngp(const DATA *data, const size_t ndata, const int Ng,
    const double *bmin, const double *Lbox, FFT_REAL *rho) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    /* Boundary check has been performed already. */
    double x = (data[i].x[0] - bmin[0]) * Ng / Lbox[0];
    double y = (data[i].x[1] - bmin[1]) * Ng / Lbox[1];
    double z = (data[i].x[2] - bmin[2]) * Ng / Lbox[2];
    int x0 = (int) x;
    int y0 = (int) y;
    int z0 = (int) z;

    if (x - x0 >= 0.5) x0 = (x0 == Ng - 1) ? 0 : x0 + 1;
    if (y - y0 >= 0.5) y0 = (y0 == Ng - 1) ? 0 : y0 + 1;
    if (z - z0 >= 0.5) z0 = (z0 == Ng - 1) ? 0 : z0 + 1;

OMP_ATOMIC
    rho[IDX(Ng,x0,y0,z0)] += data[i].w;
  }
}

/******************************************************************************
Function `cic`:
  Assign particles to the mesh with the Could-In-Cell scheme.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the sample;
  * `Ng`:       number of grid cells per box side;
  * `bmin`:     lower boundaries of the box;
  * `Lbox`:     side length of the box;
  * `rho`:      the mesh.
******************************************************************************/
static void cic(const DATA *data, const size_t ndata, const int Ng,
    const double *bmin, const double *Lbox, FFT_REAL *rho) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    /* Boundary check has been performed already. */
    double x = (data[i].x[0] - bmin[0]) * Ng / Lbox[0];
    double y = (data[i].x[1] - bmin[1]) * Ng / Lbox[1];
    double z = (data[i].x[2] - bmin[2]) * Ng / Lbox[2];
    int x0 = (int) x;
    int y0 = (int) y;
    int z0 = (int) z;

    /* Weights for neighbours. */
    double wx1 = x - x0;
    double wx0 = 1 - wx1;
    double wy1 = y - y0;
    double wy0 = 1 - wy1;
    double wz1 = z - z0;
    double wz0 = 1 - wz1;

    int x1 = (x0 == Ng - 1) ? 0 : x0 + 1;
    int y1 = (y0 == Ng - 1) ? 0 : y0 + 1;
    int z1 = (z0 == Ng - 1) ? 0 : z0 + 1;

    wx0 *= data[i].w;
    wx1 *= data[i].w;

OMP_ATOMIC
    rho[IDX(Ng,x0,y0,z0)] += wx0 * wy0 * wz0;
OMP_ATOMIC
    rho[IDX(Ng,x0,y0,z1)] += wx0 * wy0 * wz1;
OMP_ATOMIC
    rho[IDX(Ng,x0,y1,z0)] += wx0 * wy1 * wz0;
OMP_ATOMIC
    rho[IDX(Ng,x0,y1,z1)] += wx0 * wy1 * wz1;
OMP_ATOMIC
    rho[IDX(Ng,x1,y0,z0)] += wx1 * wy0 * wz0;
OMP_ATOMIC
    rho[IDX(Ng,x1,y0,z1)] += wx1 * wy0 * wz1;
OMP_ATOMIC
    rho[IDX(Ng,x1,y1,z0)] += wx1 * wy1 * wz0;
OMP_ATOMIC
    rho[IDX(Ng,x1,y1,z1)] += wx1 * wy1 * wz1;
  }
}

/******************************************************************************
Function `tsc`:
  Assign particles to the mesh with the Triangular Shaped Cloud scheme.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the sample;
  * `Ng`:       number of grid cells per box side;
  * `bmin`:     lower boundaries of the box;
  * `Lbox`:     side length of the box;
  * `rho`:      the mesh.
******************************************************************************/
static void tsc(const DATA *data, const size_t ndata, const int Ng,
    const double *bmin, const double *Lbox, FFT_REAL *rho) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    double w0[3], w1[3], w2[3];
    int x0[3], x1[3], x2[3];
    for (int k = 0; k < 3; k++) {       /* let the compiler expand the loop */
      double x = (data[i].x[k] - bmin[k]) * Ng / Lbox[k];
      x0[k] = (int) x;
      double d = x - x0[k];

      /* keep x0,x1,x2 left-to-right to reduce cache miss */
      if (d < 0.5) {
        x1[k] = x0[k];
        x0[k] = (x1[k] == 0) ? Ng - 1 : x1[k] - 1;
        x2[k] = (x1[k] == Ng - 1) ? 0 : x1[k] + 1;
        w0[k] = 0.5 - d;
        w0[k] *= w0[k] * 0.5;
        w1[k] = 0.75 - d * d;
      }
      else {
        x1[k] = (x0[k] == Ng - 1) ? 0 : x0[k] + 1;
        x2[k] = (x1[k] == Ng - 1) ? 0 : x1[k] + 1;
        d = 1 - d;      /* replace one multiplication by subtraction */
        w0[k] = 0.5 + d;
        w0[k] *= w0[k] * 0.5;
        w1[k] = 0.75 - d * d;
      }
      w2[k] = 1 - w0[k] - w1[k];
    }
    w0[0] *= data[i].w;
    w1[0] *= data[i].w;
    w2[0] *= data[i].w;

OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x0[2])] += w0[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x1[2])] += w0[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x2[2])] += w0[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x0[2])] += w0[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x1[2])] += w0[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x2[2])] += w0[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x0[2])] += w0[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x1[2])] += w0[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x2[2])] += w0[0] * w2[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x0[2])] += w1[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x1[2])] += w1[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x2[2])] += w1[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x0[2])] += w1[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x1[2])] += w1[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x2[2])] += w1[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x0[2])] += w1[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x1[2])] += w1[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x2[2])] += w1[0] * w2[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x0[2])] += w2[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x1[2])] += w2[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x2[2])] += w2[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x0[2])] += w2[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x1[2])] += w2[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x2[2])] += w2[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x0[2])] += w2[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x1[2])] += w2[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x2[2])] += w2[0] * w2[1] * w2[2];
  }
}

/******************************************************************************
Function `pcs`:
  Assign particles to the mesh with the Piecewise Cubic Spline scheme.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the sample;
  * `Ng`:       number of grid cells per box side;
  * `bmin`:     lower boundaries of the box;
  * `Lbox`:     side length of the box;
  * `rho`:      the mesh.
******************************************************************************/
static void pcs(const DATA *data, const size_t ndata, const int Ng,
    const double *bmin, const double *Lbox, FFT_REAL *rho) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    double w0[3], w1[3], w2[3], w3[3];
    int x0[3], x1[3], x2[3], x3[3];
    for (int k = 0; k < 3; k++) {
      double x = (data[i].x[k] - bmin[k]) * Ng / Lbox[k];
      x1[k] = (int) x;          /* ensure the left-to-right order of xn */

      x0[k] = (x1[k] == 0) ? Ng - 1 : x1[k] - 1;
      x2[k] = (x1[k] == Ng - 1) ? 0 : x1[k] + 1;
      x3[k] = (x2[k] == Ng - 1) ? 0 : x2[k] + 1;

      double d = x - x1[k];
      double d2 = d * d;

      /* Division by 6 done later for all nodes. */
      w3[k] = d2 * d;
      w2[k] = 1 + 3 * (d + d2 - w3[k]);
      w1[k] = 4 - 6 * d2 + 3 * w3[k];
      w0[k] = 6 - w1[k] - w2[k] - w3[k];
    }

    /* Divide by 6^3: 0x1.2f684bda12f68p-8. */
    double norm = data[i].w * 0x1.2f684bda12f68p-8;
    w0[0] *= norm;
    w1[0] *= norm;
    w2[0] *= norm;
    w3[0] *= norm;

OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x0[2])] += w0[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x1[2])] += w0[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x2[2])] += w0[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x0[1],x3[2])] += w0[0] * w0[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x0[2])] += w0[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x1[2])] += w0[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x2[2])] += w0[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x1[1],x3[2])] += w0[0] * w1[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x0[2])] += w0[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x1[2])] += w0[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x2[2])] += w0[0] * w2[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x2[1],x3[2])] += w0[0] * w2[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x3[1],x0[2])] += w0[0] * w3[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x3[1],x1[2])] += w0[0] * w3[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x3[1],x2[2])] += w0[0] * w3[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x0[0],x3[1],x3[2])] += w0[0] * w3[1] * w3[2];

OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x0[2])] += w1[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x1[2])] += w1[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x2[2])] += w1[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x0[1],x3[2])] += w1[0] * w0[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x0[2])] += w1[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x1[2])] += w1[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x2[2])] += w1[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x1[1],x3[2])] += w1[0] * w1[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x0[2])] += w1[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x1[2])] += w1[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x2[2])] += w1[0] * w2[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x2[1],x3[2])] += w1[0] * w2[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x3[1],x0[2])] += w1[0] * w3[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x3[1],x1[2])] += w1[0] * w3[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x3[1],x2[2])] += w1[0] * w3[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x1[0],x3[1],x3[2])] += w1[0] * w3[1] * w3[2];

OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x0[2])] += w2[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x1[2])] += w2[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x2[2])] += w2[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x0[1],x3[2])] += w2[0] * w0[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x0[2])] += w2[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x1[2])] += w2[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x2[2])] += w2[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x1[1],x3[2])] += w2[0] * w1[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x0[2])] += w2[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x1[2])] += w2[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x2[2])] += w2[0] * w2[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x2[1],x3[2])] += w2[0] * w2[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x3[1],x0[2])] += w2[0] * w3[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x3[1],x1[2])] += w2[0] * w3[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x3[1],x2[2])] += w2[0] * w3[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x2[0],x3[1],x3[2])] += w2[0] * w3[1] * w3[2];

OMP_ATOMIC
    rho[IDX(Ng,x3[0],x0[1],x0[2])] += w3[0] * w0[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x0[1],x1[2])] += w3[0] * w0[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x0[1],x2[2])] += w3[0] * w0[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x0[1],x3[2])] += w3[0] * w0[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x1[1],x0[2])] += w3[0] * w1[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x1[1],x1[2])] += w3[0] * w1[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x1[1],x2[2])] += w3[0] * w1[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x1[1],x3[2])] += w3[0] * w1[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x2[1],x0[2])] += w3[0] * w2[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x2[1],x1[2])] += w3[0] * w2[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x2[1],x2[2])] += w3[0] * w2[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x2[1],x3[2])] += w3[0] * w2[1] * w3[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x3[1],x0[2])] += w3[0] * w3[1] * w0[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x3[1],x1[2])] += w3[0] * w3[1] * w1[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x3[1],x2[2])] += w3[0] * w3[1] * w2[2];
OMP_ATOMIC
    rho[IDX(Ng,x3[0],x3[1],x3[2])] += w3[0] * w3[1] * w3[2];
  }
}


/*============================================================================*\
                         Functions for defining the box
\*============================================================================*/

/******************************************************************************
Function `get_coord_bound`:
  Obtain lower and upper limits of coordinates from the catalogs.
Arguments:
  * `cat`:      structure for storing the catalogs;
  * `bmin`:     lower limit of the coordinates;
  * `bmax`:     upper limit of the coordinates.
******************************************************************************/
static void get_coord_bound(const CATA *cat, double min[3], double max[3]) {
  min[0] = min[1] = min[2] = DBL_MAX;
  max[0] = max[1] = max[2] = -DBL_MAX;
  
  for (int i = 0; i < cat->num; i++) {
#ifdef OMP
#pragma omp parallel
    {
      /* Thread-private minimum/maximum trackers. */
      double pmin[3];
      double pmax[3];
      pmin[0] = pmin[1] = pmin[2] = DBL_MAX;
      pmax[0] = pmax[1] = pmax[2] = -DBL_MAX;

#pragma omp for
#endif
      for (size_t n = 0; n < cat->ndata[i]; n++) {
#ifdef OMP
        if (pmin[0] > cat->data[i][n].x[0]) pmin[0] = cat->data[i][n].x[0];
        if (pmax[0] < cat->data[i][n].x[0]) pmax[0] = cat->data[i][n].x[0];
        if (pmin[1] > cat->data[i][n].x[1]) pmin[1] = cat->data[i][n].x[1];
        if (pmax[1] < cat->data[i][n].x[1]) pmax[1] = cat->data[i][n].x[1];
        if (pmin[2] > cat->data[i][n].x[2]) pmin[2] = cat->data[i][n].x[2];
        if (pmax[2] < cat->data[i][n].x[2]) pmax[2] = cat->data[i][n].x[2];
#else
        if (min[0] > cat->data[i][n].x[0]) min[0] = cat->data[i][n].x[0];
        if (max[0] < cat->data[i][n].x[0]) max[0] = cat->data[i][n].x[0];
        if (min[1] > cat->data[i][n].x[1]) min[1] = cat->data[i][n].x[1];
        if (max[1] < cat->data[i][n].x[1]) max[1] = cat->data[i][n].x[1];
        if (min[2] > cat->data[i][n].x[2]) min[2] = cat->data[i][n].x[2];
        if (max[2] < cat->data[i][n].x[2]) max[2] = cat->data[i][n].x[2];
#endif
      }
#ifdef OMP
#pragma omp for
#endif
      for (size_t n = 0; n < cat->nrand[i]; n++) {
#ifdef OMP
        if (pmin[0] > cat->rand[i][n].x[0]) pmin[0] = cat->rand[i][n].x[0];
        if (pmax[0] < cat->rand[i][n].x[0]) pmax[0] = cat->rand[i][n].x[0];
        if (pmin[1] > cat->rand[i][n].x[1]) pmin[1] = cat->rand[i][n].x[1];
        if (pmax[1] < cat->rand[i][n].x[1]) pmax[1] = cat->rand[i][n].x[1];
        if (pmin[2] > cat->rand[i][n].x[2]) pmin[2] = cat->rand[i][n].x[2];
        if (pmax[2] < cat->rand[i][n].x[2]) pmax[2] = cat->rand[i][n].x[2];
#else
        if (min[0] > cat->rand[i][n].x[0]) min[0] = cat->rand[i][n].x[0];
        if (max[0] < cat->rand[i][n].x[0]) max[0] = cat->rand[i][n].x[0];
        if (min[1] > cat->rand[i][n].x[1]) min[1] = cat->rand[i][n].x[1];
        if (max[1] < cat->rand[i][n].x[1]) max[1] = cat->rand[i][n].x[1];
        if (min[2] > cat->rand[i][n].x[2]) min[2] = cat->rand[i][n].x[2];
        if (max[2] < cat->rand[i][n].x[2]) max[2] = cat->rand[i][n].x[2];
#endif
      }
#ifdef OMP
#pragma omp critical
      /* Gather the information from all threads. */
      if (min[0] > pmin[0]) min[0] = pmin[0];
      if (max[0] < pmax[0]) max[0] = pmax[0];
      if (min[1] > pmin[1]) min[1] = pmin[1];
      if (max[1] < pmax[1]) max[1] = pmax[1];
      if (min[2] > pmin[2]) min[2] = pmin[2];
      if (max[2] < pmax[2]) max[2] = pmax[2];
    }
#endif
  }
}

/******************************************************************************
Function `def_box`:
  Define the comoving box for the catalog to be placed in.
Arguments:
  * `cat`:      structure for storing the catalogs;
  * `issim`:    indicate whether the catalog is from a simulation box;
  * `bsize`:    box size set via configurations;
  * `bpad`:     box padding length set via configurations;
  * `size`:     the evaluated box size, if `bsize` is not given;
  * `bmin`:     the evaluated lower boundary of the box;
  * `max`:      the upper boundary of coordinates;
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int def_box(const CATA *cat, const bool issim, const double *bsize,
    const double *bpad, double size[3], double bmin[3], double max[3],
    const int verb) {
  double min[3];
  
  /* Obtain the boundaries of the input catalogs. */
  get_coord_bound(cat, min, max);
  const char c[3] = {'x', 'y', 'z'};

  for (int i = 0; i < 3; i++) {
    if (min[i] > max[i]) {
      P_ERR("invalid %c coordinate value %lf %lf in the catalogs\n", c[i], min[i], max[i]);
      return POWSPEC_ERR_CATA;
    }
    /* Verify coordinates of the simulation box. */
    if (issim) {
      if (min[i] < 0) {
        P_ERR("%c coordinate below 0: %lf\n", c[i], min[i]);
        return POWSPEC_ERR_CATA;
      }
      if (max[i] >= bsize[i]) {
        P_ERR("%c coordinate not smaller than BOX_SIZE: %lf\n", c[i], max[i]);
        return POWSPEC_ERR_CATA;
      }
    }
  }
  if (issim) {
    size[0] = bsize[0];
    size[1] = bsize[1];
    size[2] = bsize[2];
    bmin[0] = bmin[1] = bmin[2] = 0;
    return 0;
  }

  /* `BOX_SIZE` is set. */
  if (bsize) {
    for (int i = 0; i < 3; i++) {
      if (max[i] - min[i] > bsize[i]) {
        P_ERR("BOX_SIZE is too small for the %c coordinates,"
            " should be at least " OFMT_DBL "\n", c[i], max[i] - min[i]);
        return POWSPEC_ERR_CATA;
      }
      bmin[i] = (max[i] + min[i] - bsize[i]) * 0.5;
      size[i] = bsize[i];
    }

    if (verb) {
      printf("  Box size: [%lg, %lg, %lg]\n"
          "  Box boundaries: [[%lg,%lg], [%lg,%lg], [%lg,%lg]]\n",
          size[0], size[1], size[2], bmin[0], bmin[0] + size[0],
          bmin[1], bmin[1] + size[1], bmin[2], bmin[2] + size[2]);
    }
  }
  /* Determine box size automatically. */
  else {
    for (int i = 0; i < 3; i++) {
      size[i] = ceil((max[i] - min[i]) * (1 + bpad[i]) / POWSPEC_BOX_CEIL)
        * POWSPEC_BOX_CEIL;
      bmin[i] = (max[i] + min[i] - size[i]) * 0.5;
    }

    if (verb) {
      printf("  Box size: [%lg, %lg, %lg]\n"
          "  Box boundaries: [[%lg,%lg], [%lg,%lg], [%lg,%lg]]\n",
          size[0], size[1], size[2], bmin[0], bmin[0] + size[0],
          bmin[1], bmin[1] + size[1], bmin[2], bmin[2] + size[2]);
    }
  }

  return 0;
}

/******************************************************************************
Function `shift_cat`:
  Shift the catalog near the boundary if necessary, for periodic confitions.
  This is to avoid boundary checking for the interlaced grid, which is
    generally more expensive if done with particle assignment.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the data sample;
  * `mesh`:     structure for storing the meshes and FFT plans.
******************************************************************************/
static inline void shift_cat(DATA *data, const size_t ndata, MESH *mesh) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t n = 0; n < ndata; n++) {
    if (data[n].x[0] >= mesh->smin[0] + mesh->bsize[0])
      data[n].x[0] -= mesh->bsize[0];
    if (data[n].x[1] >= mesh->smin[1] + mesh->bsize[1])
      data[n].x[1] -= mesh->bsize[1];
    if (data[n].x[2] >= mesh->smin[2] + mesh->bsize[2])
      data[n].x[2] -= mesh->bsize[2];
  }
}


/*============================================================================*\
                   Functions for handling the mesh structure
\*============================================================================*/

/******************************************************************************
Function `mesh_destroy`:
  Deconstruct the interface for meshes.
Arguments:
  * `mesh`:     the structure for storing the meshes and FFT plans.
******************************************************************************/
void mesh_destroy(MESH *mesh) {
  if (!mesh) return;
  for (int i = 0; i < mesh->num; i++) {
    if (mesh->Fr && mesh->Fr[i]) FFT_FREE(mesh->Fr[i]);
    if (mesh->Frl && mesh->Frl[i]) FFT_FREE(mesh->Frl[i]);
    if (mesh->Fk0 && mesh->Fk0[i]) FFT_FREE(mesh->Fk0[i]);
    if (mesh->Fkl && mesh->Fkl[i]) FFT_FREE(mesh->Fkl[i]);
  }
  if (mesh->Fr) free(mesh->Fr);
  if (mesh->Frl) free(mesh->Frl);
  if (mesh->Fk0) free(mesh->Fk0);
  if (mesh->Fkl) free(mesh->Fkl);
  if (mesh->Fka) FFT_FREE(mesh->Fka);
  if (mesh->alias) free(mesh->alias);

  if (mesh->r2c) FFT_DESTROY(mesh->r2c);
  if (mesh->c2r) FFT_DESTROY(mesh->c2r);
  if (mesh->fft_init) {
#ifdef OMP
    FFT_CLEAN_OMP();
#else
    FFT_CLEAN();
#endif
  }
  free(mesh);
}

/******************************************************************************
Function `mesh_init`:
  Initialise the interface of FFT meshes and plans.
Arguments:
  * `conf`:     the structure for storing all configurations.
Return:
  Address of the interface.
******************************************************************************/
static MESH *mesh_init(const CONF *conf) {
  MESH *mesh = calloc(1, sizeof *mesh);
  if (!mesh) return NULL;

  mesh->issim = conf->issim;
  mesh->has_randoms = conf->has_randoms;
  mesh->intlace = conf->intlace;
  mesh->assign = conf->assign;
  mesh->fft_init = false;
  mesh->alias = NULL;
  mesh->Fr = mesh->Frl = NULL;
  mesh->Fk0 = mesh->Fkl = NULL;
  mesh->Fka = NULL;
  mesh->r2c = mesh->c2r = NULL;

  mesh->Ng = conf->gsize;
  mesh->Ngk = (mesh->Ng >> 1) + 1;
  mesh->Ntot = (size_t) mesh->Ng * mesh->Ng * mesh->Ng;
  /* Number of complex elements for R2C FFT. */
  mesh->Ncmplx = (size_t) mesh->Ng * mesh->Ng * mesh->Ngk;
  mesh->num = conf->ndata;

  /* Allocate memory for the density fields. */
  if (!(mesh->Fr = malloc(mesh->num * sizeof(FFT_REAL *)))) {
    free(mesh); return NULL;
  }
  //if (!mesh->issim || mesh->intlace) {
  if (mesh->has_randoms || mesh->intlace) {
    if (!(mesh->Frl = malloc(mesh->num * sizeof(FFT_REAL *)))) {
      mesh_destroy(mesh); return NULL;
    }
  }
  if (!(mesh->Fk0 = malloc(mesh->num * sizeof(FFT_CMPLX *)))) {
    mesh_destroy(mesh); return NULL;
  }
  if (mesh->intlace || conf->poles[conf->npole - 1]) {
    if (!(mesh->Fkl = malloc(mesh->num * sizeof(FFT_CMPLX *)))) {
      mesh_destroy(mesh); return NULL;
    }
  }

  /* Record the total allocated memory. */
  size_t size = 0;

  for (int i = 0; i < mesh->num; i++) {
    if (!(mesh->Fr[i] = FFT_MALLOC(mesh->Ntot * sizeof(FFT_REAL)))) {
      mesh_destroy(mesh); return NULL;
    }
    size += mesh->Ntot * sizeof(FFT_REAL);
    //if (!mesh->issim || mesh->intlace) {
    if (mesh->has_randoms || mesh->intlace) {
      if (!(mesh->Frl[i] = FFT_MALLOC(mesh->Ntot * sizeof(FFT_REAL)))) {
        mesh_destroy(mesh); return NULL;
      }
      size += mesh->Ntot * sizeof(FFT_REAL);
    }
    if (!(mesh->Fk0[i] = FFT_MALLOC(mesh->Ncmplx * sizeof(FFT_CMPLX)))) {
      mesh_destroy(mesh); return NULL;
    }
    size += mesh->Ncmplx * sizeof(FFT_CMPLX);
    if (mesh->intlace || conf->poles[conf->npole - 1]) {
      if (!(mesh->Fkl[i] = FFT_MALLOC(mesh->Ncmplx * sizeof(FFT_CMPLX)))) {
        mesh_destroy(mesh); return NULL;
      }
      size += mesh->Ncmplx * sizeof(FFT_CMPLX);
    }
  }
  if (!(mesh->alias = malloc(mesh->Ncmplx * sizeof(FFT_REAL)))) {
    mesh_destroy(mesh); return NULL;
  }
  size += mesh->Ncmplx * sizeof(FFT_REAL);
  //if (!mesh->issim && conf->poles[conf->npole - 1]) {
  if (mesh->has_randoms && conf->poles[conf->npole - 1]) {
    if (!(mesh->Fka = FFT_MALLOC(mesh->Ncmplx * sizeof(FFT_CMPLX)))) {
      mesh_destroy(mesh); return NULL;
    }
    size += mesh->Ncmplx * sizeof(FFT_CMPLX);
  }

  if (conf->verbose) {
    printf("  ~ %.3g Gb memory allocated for the density fields\n",
        size / 0x1p+30);
  }

#ifdef OMP
  if (FFT_INIT_OMP() == 0) {
    mesh_destroy(mesh);
    return NULL;
  }
  FFT_PLAN_OMP(omp_get_max_threads());
#endif

  mesh->r2c = FFT_PLAN_R2C(mesh->Ng, mesh->Ng, mesh->Ng,
      mesh->Fr[0], mesh->Fk0[0], FFT_FLAG);
  if (mesh->intlace) {
    mesh->c2r = FFT_PLAN_C2R(mesh->Ng, mesh->Ng, mesh->Ng,
        mesh->Fk0[0], mesh->Fr[0], FFT_FLAG);
  }
  mesh->fft_init = true;

  return mesh;
}


/*============================================================================*\
                    Functions for distributing data to mesh
\*============================================================================*/

/******************************************************************************
Function `to_mesh`:
  Distribute a sample to the mesh.
Arguments:
  * `data`:     structure for the data sample;
  * `ndata`:    number of objects in the sample;
  * `mesh`:     structure for all the meshes;
  * `min`:      lower boundary of the box;
  * `rho`:      pointer to the density field to be generated.
******************************************************************************/
static inline void to_mesh(const DATA *data, const size_t ndata, MESH *mesh,
    const double min[3], FFT_REAL *rho) {
  switch (mesh->assign) {
    case POWSPEC_ASSIGN_NGP:
      ngp(data, ndata, mesh->Ng, min, mesh->bsize, rho);
      break;
    case POWSPEC_ASSIGN_CIC:
      cic(data, ndata, mesh->Ng, min, mesh->bsize, rho);
      break;
    case POWSPEC_ASSIGN_TSC:
      tsc(data, ndata, mesh->Ng, min, mesh->bsize, rho);
      break;
    case POWSPEC_ASSIGN_PCS:
      pcs(data, ndata, mesh->Ng, min, mesh->bsize, rho);
      break;
    default:
      P_WRN("unrecognised particle assignment scheme: %d\n", mesh->assign);
      break;
  }
}

/******************************************************************************
Function `gen_dens`:
  Generate density fields for the catalogs.
Arguments:
  * `cat`:      structure for all the information of the catalogs;
  * `mesh`:     structure for all the meshes;
  * `verb`:     indicate whether to show detailed standard outputs.
******************************************************************************/
static void gen_dens(CATA *cat, MESH *mesh, const int verb) {
  for (int i = 0; i < cat->num; i++) {
    memset(mesh->Fr[i], 0, mesh->Ntot * sizeof(FFT_REAL));
    to_mesh(cat->data[i], cat->ndata[i], mesh, mesh->min, mesh->Fr[i]);
    //if (!mesh->issim) {
    if (mesh->has_randoms) {
      memset(mesh->Frl[i], 0, mesh->Ntot * sizeof(FFT_REAL));
      to_mesh(cat->rand[i], cat->nrand[i], mesh, mesh->min, mesh->Frl[i]);

      /* Subtract the random field from the data field. */
#ifdef OMP
#pragma omp parallel for
#endif
      for (size_t n = 0; n < mesh->Ntot; n++)
        mesh->Fr[i][n] -= mesh->Frl[i][n] * cat->alpha[i];

      /* Temporary save the field to the Fourier space mesh. */
      if (mesh->intlace) {
        memcpy(mesh->Fk0[i], mesh->Fr[i], mesh->Ntot * sizeof(FFT_REAL));
      }
    }

    /* Now deal with the interlaced mesh. */
    if (mesh->intlace) {
      /* Shift the lower boundary by half the cell size. */
      for (int k = 0; k < 3; k++)
        mesh->smin[k] = mesh->min[k] - 0.5 * mesh->bsize[k] / mesh->Ng;

      /* Boundary check for the catalog. */
      bool catshift = false;
      for (int k = 0; k < 3; k++) {
        if (mesh->smin[k] + mesh->bsize[k] <= mesh->max[k]) {
          catshift = true;
          break;
        }
      }
      if (catshift) shift_cat(cat->data[i], cat->ndata[i], mesh);

      /* Save the second field for the data to `Frl`: direct for simulations */
      memset(mesh->Frl[i], 0, mesh->Ntot * sizeof(FFT_REAL));
      to_mesh(cat->data[i], cat->ndata[i], mesh, mesh->smin, mesh->Frl[i]);

      //if (!mesh->issim) {
      if (mesh->has_randoms) {
        /* Boundary check and generate mesh for the random catalog. */
        if (catshift) shift_cat(cat->rand[i], cat->nrand[i], mesh);
        memset(mesh->Fr[i], 0, mesh->Ntot * sizeof(FFT_REAL));
        to_mesh(cat->rand[i], cat->nrand[i], mesh, mesh->smin, mesh->Fr[i]);

        /* Subtract the random field, and move back the first field. */
#ifdef OMP
#pragma omp parallel for
#endif
        for (size_t n = 0; n < mesh->Ntot; n++)
          mesh->Frl[i][n] -= mesh->Fr[i][n] * cat->alpha[i];
        memcpy(mesh->Fr[i], mesh->Fk0[i], mesh->Ntot * sizeof(FFT_REAL));
      }
    }
    if (verb) {
      if (cat->num == 2)
        printf("  Density field generated with %s for catalog %d\n",
            powspec_assign_names[mesh->assign], i);
      else
        printf("  Density field generated with %s for the catalog\n",
            powspec_assign_names[mesh->assign]);
    }
  }
}


/*============================================================================*\
                      Interface for generating the meshes
\*============================================================================*/

/******************************************************************************
Function `genr_mesh`:
  Interface for mesh generation.
Arguments:
  * `conf`:     the structure for configurations;
  * `cat`:      the structure for storing information from catalogs.
Return:
  Address of the structure for all meshes and FFT plans; NULL on error.
******************************************************************************/
MESH *genr_mesh(const CONF *conf, CATA *cat) {
  printf("Generating meshes for FFT ...");
  if (!conf) {
    P_ERR("configuration parameters not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

  if (!cat) {
    P_ERR("catalogs not read\n");
    return NULL;
  }
  printf("TEST\n"); fflush(stdout);
  /* Initialise the meshes. */
  MESH *mesh = mesh_init(conf);
  if (!mesh) {
    P_ERR("failed to initalise the meshes\n");
    return NULL;
  }
  

  /* Define the box. */
  if (def_box(cat, conf->issim, conf->bsize, conf->bpad, mesh->bsize,
        mesh->min, mesh->max, conf->verbose))
    return NULL;
  
  /* Generate the density fields. */
  gen_dens(cat, mesh, conf->verbose);

  /* Compute normalisation factor and shot noise for simulation boxes. */
  if (conf->issim) {
    for (int i = 0; i < cat->num; i++) {
      double vol = mesh->bsize[0] * mesh->bsize[1] * mesh->bsize[2];
      cat->shot[i] = vol / cat->wdata[i];
      cat->norm[i] = cat->wdata[i] * cat->wdata[i] / vol;
      //printf("%lf %lf %lf %lf", cat->shot[i], cat->norm[i], cat->num, cat->wdata[i]);
    }
  }

  /* Reuse `mesh->smin` for the coordinate on mesh of the lowest box corner. */
  for (int i = 0; i < 3; i++)
    mesh->smin[i] = mesh->min[i] * mesh->Ng / mesh->bsize[i];

  /* Release memory for the catalogs (meta data is still necessary). */
  for (int i = 0; i < cat->num; i++) {
    if (cat->data[i]) free(cat->data[i]);
    if (cat->rand[i]) free(cat->rand[i]);
  }
  free(cat->data); cat->data = NULL;
  free(cat->rand); cat->rand = NULL;

  printf(FMT_DONE);
  return mesh;
}

