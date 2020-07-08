/*******************************************************************************
* cnvt_coord.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include "cspline.h"
#include "legauss.h"
#include "cnvt_coord.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                 Data structure for cubic spline interpolation
\*============================================================================*/

typedef struct {
  size_t nsp;           /* number of sample points, excluding (0,0)      */
  double *z;            /* redshifts                                     */
  double *d;            /* radial comoving distances                     */
  double *ypp;          /* second derivative for spline interpolation    */
} COORD_CNVT;


/*============================================================================*\
             Functions for reading and interpolating sample points
\*============================================================================*/

/******************************************************************************
Function `cnvt_init`:
  Initialise the structure for coordinate conversion.
Return:
  Address of the structure.
******************************************************************************/
static COORD_CNVT *cnvt_init(void) {
  COORD_CNVT *cnvt = malloc(sizeof *cnvt);
  if (!cnvt) return NULL;
  cnvt->nsp = 0;
  cnvt->z = cnvt->d = cnvt->ypp = NULL;
  return cnvt;
}

/******************************************************************************
Function `cnvt_destroy`:
  Deconstruct the structure for coordinate conversion.
Arguments:
  * `cnvt`:     the structure to be deconstrcuted.
******************************************************************************/
static void cnvt_destroy(COORD_CNVT *cnvt) {
  if (!cnvt) return;
  if (cnvt->z) free(cnvt->z);
  if (cnvt->d) free(cnvt->d);
  if (cnvt->ypp) free(cnvt->ypp);
  free(cnvt);
}

/******************************************************************************
Function `read_sample`:
  Read sample points from a file.
Arguments:
  * `fname`:    name of the file for coordinate conversion;
  * `cnvt`:     the structure for coordinate conversion;
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_sample(const char *fname, COORD_CNVT *cnvt, const int verb) {
  if (!fname || !cnvt) return POWSPEC_ERR_ARG;
  if (read_ascii_simple(fname, &cnvt->z, &cnvt->d, &cnvt->nsp, verb))
    return POWSPEC_ERR_FILE;
  /* TODO: sort the sample points with probably TimSort. */
  return 0;
}

/******************************************************************************
Function `sample_interp`:
  Interpolate the coordinate conversion samples with binary search.
Arguments:
  * `cnvt`:     data structure for storing sample points;
  * `z`:        the redshift to be converted to a comoving distance.
Return:
  The radial comoving distance on success; HUGE_VAL on error.
******************************************************************************/
static inline double sample_interp(const COORD_CNVT *cnvt, const double z) {
  size_t i, l, u;
  i = l = 0;
  u = cnvt->nsp - 1;
  if (z < cnvt->z[l] || z >= cnvt->z[u]) return HUGE_VAL;

  while (l <= u) {
    i = (l + u) >> 1;
    if (cnvt->z[i + 1] <= z) l = i + 1;
    else if (cnvt->z[i] > z) u = i - 1;
    else break;
  }

  return cspline_eval(cnvt->z, cnvt->d, cnvt->ypp, z, i);
}

/******************************************************************************
Function `cnvt_cata_interp`:
  Apply coordinate conversion for a catalog with interpolation.
Arguments:
  * `cnvt`:     data structure for storing sample points;
  * `data`:     data strucutre for the coordinates;
  * `ndata`:    number of objects.
******************************************************************************/
static void cnvt_cata_interp(const COORD_CNVT *cnvt, DATA *data,
    const size_t ndata) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    double ra = data[i].x[0] * DEGREE_2_RAD;
    double dec = data[i].x[1] * DEGREE_2_RAD;
    double dist = sample_interp(cnvt, data[i].x[2]);

    data[i].x[0] = dist * cos(dec) * cos(ra);
    data[i].x[1] = dist * cos(dec) * sin(ra);
    data[i].x[2] = dist * sin(dec);
  }
}


/*============================================================================*\
                         Functions for the integration
\*============================================================================*/

/******************************************************************************
Function `cnvt_z_sample`:
  Sample redshifts uniformly in the redshift range of the catalogs.
Arguments:
  * `conf`:     the structure for storing configurations;
  * `cat`:      the structure for storing information from catalogs;
  * `num`:      number of sample points for the redshifts.
Return:
  Address of the sampled redshift array.
******************************************************************************/
static double *cnvt_z_sample(const CONF *conf, const CATA *cat, const int num) {
  if (num < 2) return NULL;
  double zmax = -DBL_MAX;
  double zmin = DBL_MAX;

#ifdef OMP
  const int nomp = omp_get_max_threads();
  double *pmax;                 /* thread-private maximum redshift */
  double *pmin;                 /* thread-private minimum redshift */
  if (!(pmax = calloc(nomp, sizeof(double)))) {
    P_ERR("failed to allocate memory for thread-private redshift.\n");
    return NULL;
  }
  if (!(pmin = calloc(nomp, sizeof(double)))) {
    P_ERR("failed to allocate memory for thread-private redshift.\n");
    free(pmax); return NULL;
  }
  for (int i = 0; i < nomp; i++) {
    pmax[i] = -DBL_MAX;
    pmin[i] = DBL_MAX;
  }
#endif

  /* Determine the minimum and maximum redshift of the samples. */
  for (int i = 0; i < cat->num; i++) {
    bool iscnvt = (conf->dcnvt) ? conf->dcnvt[i] : DEFAULT_CONVERT;
    /* Check the data catalog. */
    if (iscnvt) {
#ifdef OMP
#pragma omp parallel num_threads(nomp)
      {
        const int tid = omp_get_thread_num();
#pragma omp for
#endif
        for (size_t n = 0; n < cat->ndata[i]; n++) {
          double z = cat->data[i][n].x[2];
          if (z < 0) {
#ifdef OMP
#pragma omp critical
            {
#endif
              P_ERR("invalid negative redshift in the data catalog:\n"
                  "(%g, %g, %g)\n", cat->data[i][n].x[0],
                  cat->data[i][n].x[1], cat->data[i][n].x[2]);
#ifdef OMP
              exit(POWSPEC_ERR_CATA);
            }
#else
            return NULL;
#endif
          }
#ifdef OMP
          if (pmax[tid] < z) pmax[tid] = z;
          if (pmin[tid] > z) pmin[tid] = z;
#else
          if (zmax < z) zmax = z;
          if (zmin > z) zmin = z;
#endif
        }
#ifdef OMP
      }
#endif
    }

    iscnvt = (conf->rcnvt) ? conf->rcnvt[i] : DEFAULT_CONVERT;
    /* Check the random catalog. */
    if (iscnvt) {
#ifdef OMP
#pragma omp parallel num_threads(nomp)
      {
        const int tid = omp_get_thread_num();
#pragma omp for
#endif
        for (size_t n = 0; n < cat->nrand[i]; n++) {
          double z = cat->rand[i][n].x[2];
          if (z < 0) {
#ifdef OMP
#pragma omp critical
            {
#endif
              P_ERR("invalid negative redshift in the random catalog:\n"
                  "(%g, %g, %g)\n", cat->rand[i][n].x[0],
                  cat->rand[i][n].x[1], cat->rand[i][n].x[2]);
#ifdef OMP
              exit(POWSPEC_ERR_CATA);
            }
#else
            return NULL;
#endif
          }
#ifdef OMP
          if (pmax[tid] < z) pmax[tid] = z;
          if (pmin[tid] > z) pmin[tid] = z;
#else
          if (zmax < z) zmax = z;
          if (zmin > z) zmin = z;
#endif
        }
#ifdef OMP
      }
#endif
    }
  }
#ifdef OMP
  /* Gather the largest redshift from all threads. */
  for (int i = 0; i < nomp; i++) {
    if (zmax < pmax[i]) zmax = pmax[i];
    if (zmin > pmin[i]) zmin = pmin[i];
  }
  free(pmax);
  free(pmin);
#endif

  if (zmin > zmax) {
    P_ERR("invalid redshift value in the catalogs.\n");
    return NULL;
  }

  /* Generate sample points. */
  double *zsp = malloc(num * sizeof(double));
  if (!zsp) {
    P_ERR("failed to allocate memory for sample points of redshifts.\n");
    return NULL;
  }
  for (int i = 0; i < num; i++)
    zsp[i] = zmin + i * (zmax - zmin) / (num - 1);

  return zsp;
}

/******************************************************************************
Function `cnvt_integrand`:
  Integrand for the redshift to radial comoving distance conversion.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `z`:        the redshift to be converted to comoving distance;
Return:
  The integrand for comoving distance integration.
******************************************************************************/
static inline double cnvt_integrand(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const double z) {
  double z1 = z + 1;
  double z2 = z1 * z1;
  double d = Omega_m * z2 * z1;
  if (Omega_k) d += Omega_k * z2;
  if (widx) d += Omega_L * pow(z1, widx);
  else d += Omega_L;
  d = SPEED_OF_LIGHT * 0.01 / sqrt(d);
  return d;
}

/******************************************************************************
Function `cnvt_legauss`:
  Convert redshift to comoving distance using the Legendre-Gauss integration.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `order`:    order of the integration;
  * `z`:        redshift to be converted to comoving distance;
Return:
  The radial comoving distance on success; HUGE_VAL on error.
******************************************************************************/
static inline double cnvt_legauss(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const int order, const double z) {
  /* Variable transformation for integration from 0 to z. */
  double zp = z * 0.5;
  double sum = 0;
  int i;
  for (i = LEGAUSS_IDX(order);
      i < LEGAUSS_IDX(order) + LEGAUSS_LEN_NONZERO(order); i++) {
    /* Look up the abscissas and weights. */
    double x = legauss_x[i];
    double w = legauss_w[i];

    /* Integrate for both positive and negative abscissas. */
    double z1 = zp * (1 + x);
    double z2 = zp * (1 - x);
    sum += w * (cnvt_integrand(Omega_m, Omega_L, Omega_k, widx, z1) +
        cnvt_integrand(Omega_m, Omega_L, Omega_k, widx, z2));
  }
  /* For odd orders, there is also the abscissas x = 0. */
  if (order & 1)
    sum += legauss_w[i] * cnvt_integrand(Omega_m, Omega_L, Omega_k, widx, zp);

  return sum * zp;
}

/******************************************************************************
Function `cnvt_legauss_order`:
  Check the order of Legendre-Gauss integration given the desired precision.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `err`:      maximum allowed error for the integration;
  * `ztest`:    an array of redshifts to be tested;
  * `num`:      the length of the redshift array.
Return:
  The order on success; INT_MAX on error.
******************************************************************************/
static int cnvt_legauss_order(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const double err,
    const double *ztest, const int num) {
  int order = 0;
#ifdef OMP
#pragma omp parallel
  {
    int priv = 0;
#pragma omp for
#endif
    for (int i = 0; i < num; i++) {
      double oint, nint = 0;
      int n = LEGAUSS_MIN_ORDER - 1;
      do {
        if (n > LEGAUSS_MAX_ORDER) {
          n = INT_MAX;
          break;
        }
        oint = nint;
        nint = cnvt_legauss(Omega_m, Omega_L, Omega_k, widx, ++n, ztest[i]);
      }
      while (fabs(nint - oint) > nint * err);
#ifdef OMP
      if (priv < n) priv = n;
#else
      if (order < n) order = n;
#endif
    }
#ifdef OMP
#pragma omp critical
    if (order < priv) order = priv;
  }
#endif
  return order;
}

/******************************************************************************
Function `cnvt_cata_integr`:
  Apply coordinate conversion for a catalog with integration.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `Omega_k`:  density parameter of curvature at z = 0;
  * `widx`:     power index related to the dark energy equation of state;
  * `order`:    order of the integration;
  * `data`:     data strucutre for the coordinates;
  * `ndata`:    number of objects.
******************************************************************************/
static void cnvt_cata_integr(const double Omega_m, const double Omega_L,
    const double Omega_k, const double widx, const int order, DATA *data,
    const size_t ndata) {
#ifdef OMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < ndata; i++) {
    double ra = data[i].x[0] * DEGREE_2_RAD;
    double dec = data[i].x[1] * DEGREE_2_RAD;
    double dist =
      cnvt_legauss(Omega_m, Omega_L, Omega_k, widx, order, data[i].x[2]);

    data[i].x[0] = dist * cos(dec) * cos(ra);
    data[i].x[1] = dist * cos(dec) * sin(ra);
    data[i].x[2] = dist * sin(dec);
  }
}


/*============================================================================*\
                      Interfaces for coordinate conversion
\*============================================================================*/

/******************************************************************************
Function `cnvt_coord_interp`:
  Apply coordinate conversion for catalogs with interpolation.
Arguments:
  * `conf`:     the structure for configurations;
  * `cat`:      the structure for storing information from catalogs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cnvt_coord_interp(const CONF *conf, CATA *cat) {
  /* Create the structure for interpolation. */
  COORD_CNVT *cnvt = cnvt_init();
  if (!cnvt) {
    P_ERR("failed to allocate memory for coordinate conversion.\n");
    return POWSPEC_ERR_MEMORY;
  }

  /* Read samples from a file. */
  if (read_sample(conf->fcdst, cnvt, conf->verbose)) {
    cnvt_destroy(cnvt); return POWSPEC_ERR_FILE;
  }

  /* Creat the second derivative for interpolation. */
  cnvt->ypp = cspline_ypp(cnvt->z, cnvt->d, cnvt->nsp);
  if (!cnvt->ypp) {
    P_ERR("failed to interpolate the sample points.\n");
    return POWSPEC_ERR_CNVT;
  }

  /* Coordinate conversion. */
  for (int i = 0; i < cat->num; i++) {
    bool iscnvt = (conf->dcnvt) ? conf->dcnvt[i] : DEFAULT_CONVERT;
    if (iscnvt) {
      /* Conversion for the data catalog. */
      cnvt_cata_interp(cnvt, cat->data[i], cat->ndata[i]);
      if (conf->verbose) {
        if (cat->num == 1) printf("  Conversion done for the data catalog\n");
        else printf("  Conversion done for the data catalog %d\n", i + 1);
      }
    }

    iscnvt = (conf->rcnvt) ? conf->rcnvt[i] : DEFAULT_CONVERT;
    if (iscnvt) {
      /* Conversion for the random catalog. */
      cnvt_cata_interp(cnvt, cat->rand[i], cat->nrand[i]);
      if (conf->verbose) {
        if (cat->num == 1)
          printf("  Conversion done for the random catalog\n");
        else printf("  Conversion done for the random catalog %d\n", i + 1);
      }
    }
  }

  cnvt_destroy(cnvt);
  return 0;
}

/******************************************************************************
Function `cnvt_coord_integr`:
  Apply coordinate conversion for catalogs with integration.
Arguments:
  * `conf`:     the structure for configurations;
  * `cat`:      the structure for storing information from catalogs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cnvt_coord_integr(const CONF *conf, CATA *cat) {
  /* Sample redshifts for convergency test. */
  double *zsp = cnvt_z_sample(conf, cat, POWSPEC_INT_NUM_ZSP);
  if (!zsp) return POWSPEC_ERR_CATA;

  /* Pre-compute the power index for Lambda. */
  double widx = (conf->eos_w == -1) ? 0 : 3 * (1 + conf->eos_w);

  /* Obtain the highest order from all redshift samples. */
  int order = cnvt_legauss_order(conf->omega_m, conf->omega_l, conf->omega_k,
      widx, conf->ecdst, zsp, POWSPEC_INT_NUM_ZSP);
  if (order == INT_MAX) {
    P_ERR("failed to perform the convergency test for integrations.\n");
    return POWSPEC_ERR_CNVT;
  }

  free(zsp);
  if (conf->verbose)
    printf("  Legendre-Gauss order %d chosen for the integration error %g\n",
        order, conf->ecdst);

  /* Coordinate conversion. */
  for (int i = 0; i < cat->num; i++) {
    bool iscnvt = (conf->dcnvt) ? conf->dcnvt[i] : DEFAULT_CONVERT;
    if (iscnvt) {
      /* Conversion for the data catalog. */
      cnvt_cata_integr(conf->omega_m, conf->omega_l, conf->omega_k, widx,
          order, cat->data[i], cat->ndata[i]);
      if (conf->verbose) {
        if (cat->num == 1) printf("  Conversion done for the data catalog\n");
        else printf("  Conversion done for the data catalog %d\n", i + 1);
      }
    }

    iscnvt = (conf->rcnvt) ? conf->rcnvt[i] : DEFAULT_CONVERT;
    if (iscnvt) {
      /* Conversion for the random catalog. */
      cnvt_cata_integr(conf->omega_m, conf->omega_l, conf->omega_k, widx,
          order, cat->rand[i], cat->nrand[i]);
      if (conf->verbose) {
        if (cat->num == 1)
          printf("  Conversion done for the random catalog\n");
        else printf("  Conversion done for the random catalog %d\n", i + 1);
      }
    }
  }

  return 0;
}

/******************************************************************************
Function `cnvt_coord`:
  Interface for coordinate conversion.
Arguments:
  * `conf`:     the structure for configurations;
  * `cat`:      the structure for storing information from catalogs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cnvt_coord(const CONF *conf, CATA *cat) {
  if (!conf) {
    P_ERR("configuration parameters not loaded.\n");
    return POWSPEC_ERR_CONF;
  }
  if (!conf->cnvt) return 0;

  printf("Converting coordinates ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  if (!cat) {
    P_ERR("catalogs not read.\n");
    return POWSPEC_ERR_CATA;
  }

  /* Coordinate conversion. */
  if (conf->fcdst) {    /* binary search for samples from file */
    if (cnvt_coord_interp(conf, cat)) return POWSPEC_ERR_CNVT;
  }
  else {                /* Legendre-Gauss integration */
    if (cnvt_coord_integr(conf, cat)) return POWSPEC_ERR_CNVT;
  }

  printf(FMT_DONE);
  return 0;
}

