/*******************************************************************************
* multipole.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [GPLv3 license]

*******************************************************************************/

#include "define.h"
#include "multipole.h"
#include "load_conf.h"
#include "genr_mesh.h"
#include "spherical.h"
#include "legpoly.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                       Functions for pre-computed values
\*============================================================================*/

/******************************************************************************
  The alias correction is performed following Jing (2005):
  https://doi.org/10.1086/427087 (arXiv:astro-ph/0409240)
*******************************************************************************/

/******************************************************************************
Function `alias_corr`:
  The alias correction term for each dimension of the wave number.
Arguments:
  * `intlace`:  indicate whether interlaced grids are used;
  * `assign`:   the particle assignment scheme;
  * `pik`:      PI * k / (2 k_N), where k_N is the Nyquist frequncy.
Return:
  The (inverse of) alias correction factor.
******************************************************************************/
static inline double alias_corr(const bool intlace,
    const powspec_assign_t assign, const double pik) {
  double fac;
  if (intlace) {        /* correct only for the window function */
    if (pik > 0.01) fac = pik / sin(pik);
    /* Taylor expansion for small x:
       x / sin(x) = 1 + x^2 / 6 + 7 x^4 / 360 + 31 x^6 / 15120 + O(x^8)
     */
    else {
      double tmp = pik * pik;
      fac = 1 + 0x1.5555555555555p-3 * tmp + 0x1.3e93e93e93e94p-6 * tmp * tmp
        + 0x1.0cbb766210cbbp-9 * tmp * tmp * tmp;
    }
    fac *= fac;
    switch (assign) {
      case POWSPEC_ASSIGN_NGP: return fac;
      case POWSPEC_ASSIGN_CIC: return fac * fac;
      case POWSPEC_ASSIGN_TSC: return fac * fac * fac;
      case POWSPEC_ASSIGN_PCS: fac *= fac; return fac * fac;
      default:
        P_WRN("alias correction not available for particle assignment scheme: "
            "%d\n", assign);
        return 1;
    }
  }
  else {                /* correct for both window function and shot noise */
    if (pik > 0.01) fac = sin(pik);
    /* Taylor expansion for small x:
       sin(x) = x - x^3 / 6 + x^5 / 120 - x^7 / 5040 + O(x^9)
     */
    else {
      double tmp = pik *pik;
      fac = (1 - 0x1.5555555555555p-3 * tmp + 0x1.1111111111111p-7 * tmp * tmp
          - 0x1.a01a01a01a01ap-13 * tmp * tmp * tmp) * pik;
    }
    fac *= fac;
    switch (assign) {
      case POWSPEC_ASSIGN_NGP: return 1;
      case POWSPEC_ASSIGN_CIC:
        return 1 / (1 - 0x1.5555555555555p-1 * fac);
      case POWSPEC_ASSIGN_TSC:
        return 1 / (1 - fac + 0x1.1111111111111p-3 * fac * fac);
      case POWSPEC_ASSIGN_PCS:
        /* Not given in Jing (2005), should be:
           1 - 4/3 * sin(pik)^2 + 2/5 * sin(pik)^4 - 4/315 * sin(pik)^6
         */
        return 1 / (1 - 0x1.5555555555555p+0 * fac + 0.4 * fac * fac
            - 0x1.a01a01a01a01ap-7 * fac * fac * fac);
      default:
        P_WRN("alias correction not available for particle assignment scheme: "
            "%d\n", assign);
        return 1;
    }
  }
}

/******************************************************************************
Function `powspec_precomp`:
  Pre-compute alias corrections, number of modes, and the mean wave number
    value inside each bin.
Arguments:
  * `mesh`:     the structure for all meshes and FFT plans;
  * `pk`:       the structure for power spectra evaluation;
  * `verb`:     indicate whether to show detailed standard outputs.
******************************************************************************/
static void powspec_precomp(MESH *mesh, PK *pk, const int verb) {
  double fac, vec[3];
  fac = M_PI / mesh->Ng;
  for (int k = 0; k < 3; k++) vec[k] = 2 * M_PI / mesh->bsize[k];
#ifdef OMP
#pragma omp parallel
  {
    /* Thread-private array for mode counting. */
    double *kcnt = calloc(pk->nbin, sizeof(double));
    size_t *ccnt = calloc(pk->nbin, sizeof(size_t));
    if (!kcnt || !ccnt) {
#pragma omp critical
      {
        P_ERR("failed to allocate memory for counting modes\n");
        exit(POWSPEC_ERR_MEMORY);
      }
    }
#pragma omp for
#endif
    for (int i = 0; i < mesh->Ng; i++) {
      size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
      double ki = (i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng;
      double alias_i = alias_corr(mesh->intlace, mesh->assign, ki * fac);
      ki *= vec[0];
      double kmod_i = ki * ki;
      for (int j = 0; j < mesh->Ng; j++) {
        size_t idx_j = idx_i + j * mesh->Ngk;
        double kj = (j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng;
        double alias_j = alias_corr(mesh->intlace, mesh->assign, kj * fac);
        kj *= vec[1];
        double kmod_j = kmod_i + kj * kj;

        /* Expand k loop and treat the 0 and Nyquist frequencies separately. */
        /* k = 0 */
        size_t idx = idx_j;
        double kmod = kmod_j;
        if (pk->log) kmod = 0.5 * log10(kmod);
        else kmod = sqrt(kmod);

        /* Compare sqrt instead of squared value, for round errors. */
        if (kmod < pk->kedge[0] || kmod >= pk->kedge[pk->nbin]) {
          mesh->alias[idx] = 0;         /* mark cells that are not used */
        }
        else {
          /* Find the wave number bin. */
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
          if (kbin < 0 || kbin >= pk->nbin) {   /* in case of round error */
            mesh->alias[idx] = 0;
          }
          else {                /* this vector is inside the desired range */
            /* Record alias correction. */
            mesh->alias[idx] = alias_i * alias_j;
#ifdef OMP
            kcnt[kbin] += kmod;
            ccnt[kbin] += 1;
#else
            pk->km[kbin] += kmod;
            pk->cnt[kbin] += 1;
#endif
          }
        }

        /* 0 < k < kny: consider also the omitted negative plane. */
        for (int k = 1; k <= ((mesh->Ng - 1) >> 1); k++) {
          idx = idx_j + k;
          double alias_k = alias_corr(mesh->intlace, mesh->assign, k * fac);
          double kk = k * vec[2];

          /* `kmod` can also be pre-computed in principle,
             but the gain is small, compared to the complexity of FFT */
          kmod = kmod_j + kk * kk;
          if (pk->log) kmod = 0.5 * log10(kmod);
          else kmod = sqrt(kmod);

          /* Compare sqrt instead of squared value, for round errors. */
          if (kmod < pk->kedge[0] || kmod >= pk->kedge[pk->nbin]) {
            mesh->alias[idx] = 0;         /* mark cells that are not used */
          }
          else {
            /* Find the wave number bin. */
            int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
            if (kbin < 0 || kbin >= pk->nbin) {   /* in case of round error */
              mesh->alias[idx] = 0;
            }
            else {
              /* Record alias correction. */
              mesh->alias[idx] = alias_i * alias_j * alias_k;

              /* Count wave numbers in each bin, with the other half-space. */
#ifdef OMP
              kcnt[kbin] += 2 * kmod;
              ccnt[kbin] += 2;
#else
              pk->km[kbin] += 2 * kmod;
              pk->cnt[kbin] += 2;
#endif
            }
          }
        }

        /* k = kny */
        if ((mesh->Ng & 1) == 0) {
          int k = mesh->Ng >> 1;
          idx = idx_j + k;
          double alias_k = alias_corr(mesh->intlace, mesh->assign, k * fac);
          double kk = k * vec[2];

          kmod = kmod_j + kk * kk;
          if (pk->log) kmod = 0.5 * log10(kmod);
          else kmod = sqrt(kmod);

          if (kmod < pk->kedge[0] || kmod >= pk->kedge[pk->nbin]) {
            mesh->alias[idx] = 0;
          }
          else {
            int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
            if (kbin < 0 || kbin >= pk->nbin) mesh->alias[idx] = 0;
            else {
              mesh->alias[idx] = alias_i * alias_j * alias_k;
#ifdef OMP
              kcnt[kbin] += kmod;
              ccnt[kbin] += 1;
#else
              pk->km[kbin] += kmod;
              pk->cnt[kbin] += 1;
#endif
            }
          }
        }
      }
    }
#ifdef OMP
    /* Gather the counts from threads. */
#pragma omp critical
    for (int k = 0; k < pk->nbin; k++) {
      pk->km[k] += kcnt[k];
      pk->cnt[k] += ccnt[k];
    }
    free(kcnt);
    free(ccnt);
  }
#endif

  for (int i = 0; i < pk->nbin; i++)
    if (pk->cnt[i]) pk->km[i] /= pk->cnt[i];
  if (verb) printf("  Alias corrections and wave numbers are pre-computed\n");
}


/*============================================================================*\
                 Funtions for the power spectra data structure
\*============================================================================*/

/******************************************************************************
Function `powspec_destroy`:
  Deconstruct the interface for power spectra.
Arguments:
  * `pk`:       structure for the power spectra.
******************************************************************************/
void powspec_destroy(PK *pk) {
  if (!pk) return;
  if (pk->poles) free(pk->poles);
  if (pk->kedge) free(pk->kedge);
  if (pk->k) free(pk->k);
  if (pk->km) free(pk->km);
  if (pk->cnt) free(pk->cnt);
  if (pk->lcnt) free(pk->lcnt);
  for (int i = 0; i < 2; i++) {
    if (pk->pl[i]) {
      for (int j = 0; j < pk->nl; j++)
        if (pk->pl[i][j]) free(pk->pl[i][j]);
      free(pk->pl[i]);
    }
  }
  if (pk->xpl) {
    for (int j = 0; j < pk->nl; j++)
      if (pk->xpl[j]) free (pk->xpl[j]);
    free(pk->xpl);
  }
#ifdef OMP
  if (pk->pcnt) free(pk->pcnt);
  if (pk->plcnt) free(pk->plcnt);
#endif
  free(pk);
}

/******************************************************************************
Function `powspec_init`:
  Initialise the interface of power spectra.
Arguments:
  * `conf`:     the structure for storing all configurations;
  * `mesh`:     the structure for meshes and FFT plans.
Return:
  Address of the interface.
******************************************************************************/
static PK *powspec_init(const CONF *conf, const MESH *mesh) {
  PK *pk = calloc(1, sizeof *pk);
  if (!pk) return NULL;

  pk->issim = conf->issim; pk->log = conf->logscale;
  if (pk->issim) {
    pk->los[0] = conf->los[0];
    pk->los[1] = conf->los[1];
    pk->los[2] = conf->los[2];
  }
  pk->isauto[0] = conf->isauto[0]; pk->isauto[1] = conf->isauto[1];
  pk->iscross = conf->iscross;
  pk->nl = conf->npole;
  pk->dk = conf->kbin;
  pk->poles = NULL;
  pk->k = pk->km = pk->kedge = pk->lcnt = NULL;
  pk->pl[0] = pk->pl[1] = pk->xpl = NULL;
  pk->cnt = NULL;
#ifdef OMP
  pk->pcnt = NULL;
  pk->plcnt = NULL;
#endif

  if (!(pk->poles = malloc(pk->nl * sizeof(int)))) {
    free(pk); return NULL;
  }
  memcpy(pk->poles, conf->poles, pk->nl * sizeof(int));

  /* Find the maximum wave number to be considered. */
  double bmax = (mesh->bsize[0] > mesh->bsize[1]) ?
    mesh->bsize[0] : mesh->bsize[1];
  if (bmax < mesh->bsize[2]) bmax = mesh->bsize[2];
  double kny = M_PI * mesh->Ng / bmax;          /* Nyquist frequency */
  if (pk->log) kny = log10(kny);                /* for logarhtim scale */
  double kmax = kny;
  if (POWSPEC_KMAX_ISSET(conf->kmax) && kmax > conf->kmax) kmax = conf->kmax;

  /* Compute the number of wave number bins given the range and bin size. */
  double nbin = round((kmax - conf->kmin) / pk->dk);
  if (nbin >= INT_MAX) {
    P_ERR("too many wave number bins due to the small bin size: "
        OFMT_DBL "\n", pk->dk);
    powspec_destroy(pk); return NULL;
  }
  pk->nbin = (int) nbin;
  if (conf->kmin + pk->dk * pk->nbin > kny) pk->nbin -= 1;
  if (pk->nbin < 1) {
    P_ERR("not enough k bins given the Nyquist frequency and the bin size: "
        OFMT_DBL "\n", pk->dk);
    powspec_destroy(pk); return NULL;
  }

  if (!(pk->kedge = malloc((pk->nbin + 1) * sizeof(double)))) {
    powspec_destroy(pk); return NULL;
  }
  if (!(pk->k = malloc(pk->nbin * sizeof(double)))) {
    powspec_destroy(pk); return NULL;
  }
  if (!(pk->km = calloc(pk->nbin, sizeof(double)))) {
    powspec_destroy(pk); return NULL;
  }
  if (!(pk->cnt = calloc(pk->nbin, sizeof(size_t)))) {
    powspec_destroy(pk); return NULL;
  }
  /* Mode counts for multipoles of simulation boxes. */
  if (pk->issim) {
    if (!(pk->lcnt = malloc(sizeof(double) * pk->nbin * pk->nl))) {
      powspec_destroy(pk); return NULL;
    }
  }
#ifdef OMP
  pk->nomp = omp_get_max_threads();
  /* Array for thread-private counting. */
  if (pk->issim) {
    if (!(pk->plcnt = malloc(sizeof(double) * pk->nomp * pk->nl * pk->nbin))) {
      powspec_destroy(pk); return NULL;
    }
    pk->pcnt = malloc(sizeof(double) * pk->nomp * pk->nl * pk->nbin);
  }
  else pk->pcnt = malloc(pk->nomp * pk->nbin * sizeof(double));
  if (!pk->pcnt) {
    powspec_destroy(pk); return NULL;
  }
#endif

  /* Define the wave number bins. */
  for (int i = 0; i <= pk->nbin; i++) pk->kedge[i] = conf->kmin + pk->dk * i;
  for (int i = 0; i < pk->nbin; i++)
    pk->k[i] = (pk->kedge[i] + pk->kedge[i + 1]) * 0.5;

  /* Allocate memoy for multipoles. */
  for (int i = 0; i < 2; i++) {
    if (pk->isauto[i]) {
      if (!(pk->pl[i] = malloc(pk->nl * sizeof(double *)))) {
        powspec_destroy(pk); return NULL;
      }
      for (int j = 0; j < pk->nl; j++) {
        if (!(pk->pl[i][j] = calloc(pk->nbin, sizeof(double)))) {
          powspec_destroy(pk); return NULL;
        }
      }
    }
  }
  if (pk->iscross) {
    if (!(pk->xpl = malloc(pk->nl * sizeof(double *)))) {
      powspec_destroy(pk); return NULL;
    }
    for (int j = 0; j < pk->nl; j++) {
      if (!(pk->xpl[j] = calloc(pk->nbin, sizeof(double)))) {
        powspec_destroy(pk); return NULL;
      }
    }
  }

  return pk;
}


/*============================================================================*\
             Functions for generating Fourier space density fields
\*============================================================================*/

/******************************************************************************
Function `dens_k0`:
  Compute the Fourier space density field monopole (ell = 0).
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `verb`:     indicate whether to show detailed standard outputs.
******************************************************************************/
static void dens_k0(MESH *mesh, const int verb) {
  for (int n = 0; n < mesh->num; n++) {
    if (verb) {
      if (mesh->num != 2) printf("  Executing 1 FFT for l = 0 ...");
      else printf("  Executing 1 FFT for l = 0 with catalog %d ...", n + 1);
      fflush(stdout);
    }

    /* Fourier transform the original field. */
    FFT_EXEC_R2C(mesh->r2c, mesh->Fr[n], mesh->Fk0[n]);

    if (verb) {
      if (mesh->num != 2) printf("\r  Done with computing 1 FFT for l = 0\n");
      else printf("\r  Done with computing 1 FFT for l = 0 with catalog %d\n",
          n + 1);
    }

    if (mesh->intlace) {
      if (verb) {
        printf("  Executing 2 FFTs for grid interlacing ...");
        fflush(stdout);
      }

      /* Fourier transform the shifted field for intelacing. */
      FFT_EXEC_R2C(mesh->r2c, mesh->Frl[n], mesh->Fkl[n]);

      /* Combine the two fields in Fourier space. */
      double fac = M_PI / mesh->Ng;
#ifdef OMP
#pragma omp parallel for
#endif
      for (int i = 0; i < mesh->Ng; i++) {
        size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
        double ki = fac * ((i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng);
        for (int j = 0; j < mesh->Ng; j++) {
          size_t idx_j = idx_i + j * mesh->Ngk;
          double kj = fac * ((j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng);
          for (int k = 0; k < mesh->Ngk; k++) {
            size_t idx = idx_j + k;
            double kk = fac * k;

            double ksum = ki + kj + kk;
            /* The compiler may optimise the sin & cos calls with -O2. */
            double cosks = cos(ksum);
            double sinks = sin(ksum);

            mesh->Fk0[n][idx][0] = 0.5 * (mesh->Fk0[n][idx][0] +
                cosks * mesh->Fkl[n][idx][0] - sinks * mesh->Fkl[n][idx][1]);
            mesh->Fk0[n][idx][1] = 0.5 * (mesh->Fk0[n][idx][1] +
                sinks * mesh->Fkl[n][idx][0] + cosks * mesh->Fkl[n][idx][1]);
          }
        }
      }

      /* Copy Fk0 to Fkl, as FFTW C2R overwrites the input. */
      memcpy(mesh->Fkl[n], mesh->Fk0[n], mesh->Ncmplx * sizeof(FFT_CMPLX));

      /* Backward transform for the interlaced density field. */
      FFT_EXEC_C2R(mesh->c2r, mesh->Fkl[n], mesh->Fr[n]);

      /* Re-normalise the density field. */
#ifdef OMP
#pragma omp parallel for
#endif
      for (size_t i = 0; i < mesh->Ntot; i++)
        mesh->Fr[n][i] /= mesh->Ntot;

      if (verb) printf("\r  Done with computing 2 FFTs for grid interlacing\n");
    }
  }
}

/******************************************************************************
Function `dens_k<POWSPEC_ELL>`:
  Compute the Fourier space density field multipole (ell = <POWSPEC_ELL>).
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `verb`:     indicate whether to show detailed standard outputs.
******************************************************************************/

#ifdef POWSPEC_ELL
#undef POWSPEC_ELL
#endif
#define POWSPEC_ELL     1
extern inline double YlmR_l1(const int, const double, const double,
    const double, const double);
#include "mp_template.c"

#ifdef POWSPEC_ELL
#undef POWSPEC_ELL
#endif
#define POWSPEC_ELL     2
extern inline double YlmR_l2(const int, const double, const double,
    const double, const double);
#include "mp_template.c"

#ifdef POWSPEC_ELL
#undef POWSPEC_ELL
#endif
#define POWSPEC_ELL     3
extern inline double YlmR_l3(const int, const double, const double,
    const double, const double);
#include "mp_template.c"

#ifdef POWSPEC_ELL
#undef POWSPEC_ELL
#endif
#define POWSPEC_ELL     4
extern inline double YlmR_l4(const int, const double, const double,
    const double, const double);
#include "mp_template.c"

#ifdef POWSPEC_ELL
#undef POWSPEC_ELL
#endif
#define POWSPEC_ELL     5
extern inline double YlmR_l5(const int, const double, const double,
    const double, const double);
#include "mp_template.c"

#ifdef POWSPEC_ELL
#undef POWSPEC_ELL
#endif
#define POWSPEC_ELL     6
extern inline double YlmR_l6(const int, const double, const double,
    const double, const double);
#include "mp_template.c"


/*============================================================================*\
                   Functions for counting Fourier space modes
\*============================================================================*/

/******************************************************************************
Function `count_mode_lin`:
  Count modes of the Fourier space density field with k bins in linear scale,
    for non-periodic data.
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `Fk1`:      the first Fourier space density field;
  * `Fk2`:      the second Fourier space density field;
  * `pk`:       the structure for storing all power spectra;
  * `pl`:       the array for saving the resulting power spectrum.
******************************************************************************/
static void count_mode_lin(const MESH *mesh, const FFT_CMPLX *Fk1,
    const FFT_CMPLX *Fk2, PK *pk, double *pl) {
  double vec[3];
  for (int k = 0; k < 3; k++) vec[k] = 2 * M_PI / mesh->bsize[k];
#ifdef OMP
#pragma omp parallel num_threads(pk->nomp)
  {
    const int tid = omp_get_thread_num();
    double *pcnt = pk->pcnt + tid * pk->nbin;
    memset(pcnt, 0, pk->nbin * sizeof(double));
#pragma omp for
#endif
    for (int i = 0; i < mesh->Ng; i++) {
      size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
      double ki = vec[0] * ((i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng);
      double kmod_i = ki * ki;
      for (int j = 0; j < mesh->Ng; j++) {
        size_t idx_j = idx_i + j * mesh->Ngk;
        double kj = vec[1] * ((j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng);
        double kmod_j = kmod_i + kj * kj;

        /* Expand k loop and treat the 0 and Nyquist frequencies separately. */
        /* k = 0 */
        size_t idx = idx_j;
        if (mesh->alias[idx]) {
          double kmod = sqrt(kmod_j);
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
#ifdef OMP
          pcnt[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx];
#else
          pl[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx];
#endif
        }

        /* 0 < k < kny */
        for (int k = 1; k <= ((mesh->Ng - 1) >> 1); k++) {
          idx = idx_j + k;
          if (!mesh->alias[idx]) continue;

          double kk = vec[2] * k;
          double kmod = sqrt(kmod_j + kk * kk);
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
          /* Count twice for both half-planes. */
#ifdef OMP
          pcnt[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx] * 2;
#else
          pl[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx] * 2;
#endif
        }

        /* k = kny */
        if ((mesh->Ng & 1) == 0) {
          int k = mesh->Ng >> 1;
          idx = idx_j + k;
          if (mesh->alias[idx]) {
            double kk = vec[2] * k;
            double kmod = sqrt(kmod_j + kk * kk);
            int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
#ifdef OMP
            pcnt[kbin] += (Fk1[idx][0] * Fk2[idx][0]
                + Fk1[idx][1] * Fk2[idx][1]) * mesh->alias[idx];
#else
            pl[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
              * mesh->alias[idx];
#endif
          }
        }
      }
    }
#ifdef OMP
#pragma omp for
    for (int i = 0; i < pk->nbin; i++) {
      for (int j = 0; j < pk->nomp; j++)
        pl[i] += pk->pcnt[j * pk->nbin + i];
    }
  }
#endif
}

/******************************************************************************
Function `count_mode_log`:
  Count modes of the Fourier space density field with k bins in log scale,
    for non-periodic data.
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `Fk1`:      the first Fourier space density field;
  * `Fk2`:      the second Fourier space density field;
  * `pk`:       the structure for storing all power spectra;
  * `pl`:       the array for saving the resulting power spectrum.
******************************************************************************/
static void count_mode_log(const MESH *mesh, const FFT_CMPLX *Fk1,
    const FFT_CMPLX *Fk2, PK *pk, double *pl) {
  double vec[3];
  for (int k = 0; k < 3; k++) vec[k] = 2 * M_PI / mesh->bsize[k];
#ifdef OMP
#pragma omp parallel num_threads(pk->nomp)
  {
    const int tid = omp_get_thread_num();
    double *pcnt = pk->pcnt + tid * pk->nbin;
    memset(pcnt, 0, pk->nbin * sizeof(double));
#pragma omp for
#endif
    for (int i = 0; i < mesh->Ng; i++) {
      size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
      double ki = vec[0] * ((i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng);
      double kmod_i = ki * ki;
      for (int j = 0; j < mesh->Ng; j++) {
        size_t idx_j = idx_i + j * mesh->Ngk;
        double kj = vec[1] * ((j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng);
        double kmod_j = kmod_i + kj * kj;

        /* Expand k loop and treat the 0 and Nyquist frequencies separately. */
        /* k = 0 */
        size_t idx = idx_j;
        if (mesh->alias[idx]) {
          double kmod = 0.5 * log10(kmod_j);
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
#ifdef OMP
          pcnt[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx];
#else
          pl[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx];
#endif
        }

        /* 0 < k < kny */
        for (int k = 1; k <= ((mesh->Ng - 1) >> 1); k++) {
          size_t idx = idx_j + k;
          if (!mesh->alias[idx]) continue;

          double kk = vec[2] * k;
          double kmod = 0.5 * log10(kmod_j + kk * kk);
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
          /* Count twice for both half-planes. */
#ifdef OMP
          pcnt[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx] * 2;
#else
          pl[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx] * 2;
#endif
        }

        /* k = kny */
        if ((mesh->Ng & 1) == 0) {
          int k = mesh->Ng >> 1;
          idx = idx_j + k;
          if (mesh->alias[idx]) {
            double kk = vec[2] * k;
            double kmod = 0.5 * log10(kmod_j + kk * kk);
            int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);
#ifdef OMP
            pcnt[kbin] += (Fk1[idx][0] * Fk2[idx][0]
                + Fk1[idx][1] * Fk2[idx][1]) * mesh->alias[idx];
#else
            pl[kbin] += (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
              * mesh->alias[idx];
#endif
          }
        }
      }
    }
#ifdef OMP
#pragma omp for
    for (int i = 0; i < pk->nbin; i++) {
      for (int j = 0; j < pk->nomp; j++)
        pl[i] += pk->pcnt[i + j * pk->nbin];
    }
  }
#endif
}

/******************************************************************************
Functions `legpoly`:
  Compute the Legengre polynomial.
Arguments:
  * `ell`:      the order;
  * `x`:        the variable.
Return:
  P_ell (x)
******************************************************************************/
extern inline double legpoly(const int, const double);

/******************************************************************************
Function `count_mode_sim_lin`:
  Count modes of the Fourier space density field with k bins in linear scale,
    for meshes from cubic simulation boxes.
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `Fk1`:      the first Fourier space density field;
  * `Fk2`:      the second Fourier space density field;
  * `pk`:       the structure for storing all power spectra;
  * `pl`:       the array for saving the resulting power spectra.
******************************************************************************/
static void count_mode_sim_lin(const MESH *mesh, const FFT_CMPLX *Fk1,
    const FFT_CMPLX *Fk2, PK *pk, double **pl) {
  double vec[3];
  for (int k = 0; k < 3; k++) vec[k] = 2 * M_PI / mesh->bsize[k];
  memset(pk->lcnt, 0, sizeof(double) * pk->nl * pk->nbin);
#ifdef OMP
#pragma omp parallel num_threads(pk->nomp)
  {
    const int tid = omp_get_thread_num();
    double *pcnt = pk->pcnt + (size_t) tid * pk->nl * pk->nbin;
    memset(pcnt, 0, sizeof(double) * pk->nl * pk->nbin);
    double *lcnt = pk->plcnt + (size_t) tid * pk->nl * pk->nbin;
    memset(lcnt, 0, sizeof(double) * pk->nl * pk->nbin);
#pragma omp for
#endif
    for (int i = 0; i < mesh->Ng; i++) {
      size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
      double ki = vec[0] * ((i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng);
      double kmod_i = ki * ki;
      for (int j = 0; j < mesh->Ng; j++) {
        size_t idx_j = idx_i + j * mesh->Ngk;
        double kj = vec[1] * ((j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng);
        double kmod_j = kmod_i + kj * kj;

        /* Expand k loop and treat the 0 and Nyquist frequencies separately. */
        /* k = 0 */
        size_t idx = idx_j;
        /* Skip i = j = k = 0. */
        if (mesh->alias[idx] && kmod_j) {
          double kmod = sqrt(kmod_j);
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);

          double mu = (ki * pk->los[0] + kj * pk->los[1]) / kmod;
          double p = (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx];
          double num = 1;

          for (int l = 0; l < pk->nl; l++) {
            double leg = legpoly(pk->poles[l], mu);
#ifdef OMP
            pcnt[l * pk->nbin + kbin] += p * leg;
            lcnt[l * pk->nbin + kbin] += num * leg;
#else
            pl[l][kbin] += p * leg;
            pk->lcnt[l * pk->nbin + kbin] += num * leg;
#endif
          }
        }

        /* 0 < k < kny */
        for (int k = 1; k <= ((mesh->Ng - 1) >> 1); k++) {
          idx = idx_j + k;
          if (!mesh->alias[idx]) continue;

          double kk = vec[2] * k;
          double kmod = sqrt(kmod_j + kk * kk);
          int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);

          double mu = (ki * pk->los[0] + kj * pk->los[1] + kk * pk->los[2])
            / kmod;
          double p = (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx] * 2;
          double num = 2;

          for (int l = 0; l < pk->nl; l++) {
            /* Results for the two half-planes are canceled for odd mu. */
            if (pk->poles[l] & 1) continue;
            double leg = legpoly(pk->poles[l], mu);
#ifdef OMP
            pcnt[l * pk->nbin + kbin] += p * leg;
            lcnt[l * pk->nbin + kbin] += num * leg;
#else
            pl[l][kbin] += p * leg;
            pk->lcnt[l * pk->nbin + kbin] += num * leg;
#endif
          }
        }

        /* k = kny */
        if ((mesh->Ng & 1) == 0) {
          int k = mesh->Ng >> 1;
          idx = idx_j + k;
          if (mesh->alias[idx]) {
            double kk = vec[2] * k;
            double kmod = sqrt(kmod_j + kk * kk);
            int kbin = (int) ((kmod - pk->kedge[0]) / pk->dk);

            double mu = (ki * pk->los[0] + kj * pk->los[1] + kk * pk->los[2])
              / kmod;
            double p = (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
              * mesh->alias[idx];
            double num = 1;

            for (int l = 0; l < pk->nl; l++) {
              double leg = legpoly(pk->poles[l], mu);
#ifdef OMP
              pcnt[l * pk->nbin + kbin] += p * leg;
              lcnt[l * pk->nbin + kbin] += num * leg;
#else
              pl[l][kbin] += p * leg;
              pk->lcnt[l * pk->nbin + kbin] += num * leg;
#endif
            }
          }
        }
      }
    }
#ifdef OMP
#pragma omp for
    for (int i = 0; i < pk->nbin; i++) {
      for (int l = 0; l < pk->nl; l++) {
        for (int j = 0; j < pk->nomp; j++) {
          pl[l][i] += pk->pcnt[i + pk->nbin * (l + pk->nl * j)];
          pk->lcnt[l * pk->nbin + i] +=
            pk->plcnt[i + pk->nbin * (l + pk->nl * j)];
        }
      }
    }
  }
#endif
}

/******************************************************************************
Function `count_mode_sim_log`:
  Count modes of the Fourier space density field with k bins in log scale,
    for meshes from cubic simulation boxes.
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `Fk1`:      the first Fourier space density field;
  * `Fk2`:      the second Fourier space density field;
  * `pk`:       the structure for storing all power spectra;
  * `pl`:       the array for saving the resulting power spectra.
******************************************************************************/
static void count_mode_sim_log(const MESH *mesh, const FFT_CMPLX *Fk1,
    const FFT_CMPLX *Fk2, PK *pk, double **pl) {
  double vec[3];
  for (int k = 0; k < 3; k++) vec[k] = 2 * M_PI / mesh->bsize[k];
  memset(pk->lcnt, 0, sizeof(double) * pk->nl * pk->nbin);
#ifdef OMP
#pragma omp parallel num_threads(pk->nomp)
  {
    const int tid = omp_get_thread_num();
    double *pcnt = pk->pcnt + (size_t) tid * pk->nl * pk->nbin;
    memset(pcnt, 0, sizeof(double) * pk->nl * pk->nbin);
    double *lcnt = pk->plcnt + (size_t) tid * pk->nl * pk->nbin;
    memset(lcnt, 0, sizeof(double) * pk->nl * pk->nbin);
#pragma omp for
#endif
    for (int i = 0; i < mesh->Ng; i++) {
      size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
      double ki = vec[0] * ((i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng);
      double kmod_i = ki * ki;
      for (int j = 0; j < mesh->Ng; j++) {
        size_t idx_j = idx_i + j * mesh->Ngk;
        double kj = vec[1] * ((j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng);
        double kmod_j = kmod_i + kj * kj;

        /* Expand k loop and treat the 0 and Nyquist frequencies separately. */
        /* k = 0 */
        size_t idx = idx_j;
        /* Skip i = j = k = 0. */
        if (mesh->alias[idx] && kmod_j) {
          double kmod = sqrt(kmod_j);
          int kbin = (int) ((0.5 * log10(kmod_j) - pk->kedge[0]) / pk->dk);

          double mu = (ki * pk->los[0] + kj * pk->los[1]) / kmod;
          double p = (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx];
          double num = 1;

          for (int l = 0; l < pk->nl; l++) {
            double leg = legpoly(pk->poles[l], mu);
#ifdef OMP
            pcnt[l * pk->nbin + kbin] += p * leg;
            lcnt[l * pk->nbin + kbin] += num * leg;
#else
            pl[l][kbin] += p * leg;
            pk->lcnt[l * pk->nbin + kbin] += num * leg;
#endif
          }
        }

        /* 0 < k < kny */
        for (int k = 1; k <= ((mesh->Ng - 1) >> 1); k++) {
          idx = idx_j + k;
          if (!mesh->alias[idx]) continue;

          double kk = vec[2] * k;
          double kmod = kmod_j + kk * kk;
          int kbin = (int) ((0.5 * log10(kmod) - pk->kedge[0]) / pk->dk);

          double mu = (ki * pk->los[0] + kj * pk->los[1] + kk * pk->los[2])
            / sqrt(kmod);
          double p = (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
            * mesh->alias[idx] * 2;
          double num = 2;

          for (int l = 0; l < pk->nl; l++) {
            /* Results for the two half-planes are canceled for odd mu. */
            if (pk->poles[l] & 1) continue;
            double leg = legpoly(pk->poles[l], mu);
#ifdef OMP
            pcnt[l * pk->nbin + kbin] += p * leg;
            lcnt[l * pk->nbin + kbin] += num * leg;
#else
            pl[l][kbin] += p * leg;
            pk->lcnt[l * pk->nbin + kbin] += num * leg;
#endif
          }
        }

        /* k = kny */
        if ((mesh->Ng & 1) == 0) {
          int k = mesh->Ng >> 1;
          idx = idx_j + k;
          if (mesh->alias[idx]) {
            double kk = vec[2] * k;
            double kmod = kmod_j + kk * kk;
            int kbin = (int) ((0.5 * log10(kmod) - pk->kedge[0]) / pk->dk);

            double mu = (ki * pk->los[0] + kj * pk->los[1] + kk * pk->los[2])
              / sqrt(kmod);
            double p = (Fk1[idx][0] * Fk2[idx][0] + Fk1[idx][1] * Fk2[idx][1])
              * mesh->alias[idx];
            double num = 1;

            for (int l = 0; l < pk->nl; l++) {
              double leg = legpoly(pk->poles[l], mu);
#ifdef OMP
              pcnt[l * pk->nbin + kbin] += p * leg;
              lcnt[l * pk->nbin + kbin] += num * leg;
#else
              pl[l][kbin] += p * leg;
              pk->lcnt[l * pk->nbin + kbin] += num * leg;
#endif
            }
          }
        }
      }
    }
#ifdef OMP
#pragma omp for
    for (int i = 0; i < pk->nbin; i++) {
      for (int l = 0; l < pk->nl; l++) {
        for (int j = 0; j < pk->nomp; j++) {
          pl[l][i] += pk->pcnt[i + pk->nbin * (l + pk->nl * j)];
          pk->lcnt[l * pk->nbin + i] +=
            pk->plcnt[i + pk->nbin * (l + pk->nl * j)];
        }
      }
    }
  }
#endif
}

/******************************************************************************
Function `count_mode`:
  Count modes of the Fourier space density field, and normalise power spectra.
Arguments:
  * `cat`:      the structure for information of catalogs;
  * `mesh`:     the structure for meshes and FFT plans;
  * `pk`:       the structure for storing all power spectra;
  * `idx`:      index of the catalog (-1 for cross correlation);
  * `ell`:      order of the multipole;
  * `verb`:     indicate whether to show detailed standard outputs.
******************************************************************************/
static void count_mode(const CATA *cat, const MESH *mesh, PK *pk,
    const int idx, const int ell, const int verb) {
  /* Deal with simulation boxes separately. */
  if (pk->issim) {
    if (idx < 0) {      /* cross power spectrum */
      if (pk->log)
        count_mode_sim_log(mesh, mesh->Fk0[0], mesh->Fk0[1], pk, pk->xpl);
      else
        count_mode_sim_lin(mesh, mesh->Fk0[0], mesh->Fk0[1], pk, pk->xpl);
      /* Normalisation. */
      for (int l = 0; l < pk->nl; l++) {
        for (int i = 0; i < pk->nbin; i++) {
          if (pk->cnt[i])
            pk->xpl[l][i] /= sqrt(cat->norm[0] * cat->norm[1]) * pk->cnt[i];
          pk->xpl[l][i] *= 2 * pk->poles[l] + 1;
        }
      }
      if (verb) printf("  Cross power spectra computed\n");
    }
    else {              /* auto power spectrum */
      if (pk->log) {
        count_mode_sim_log(mesh, mesh->Fk0[idx], mesh->Fk0[idx], pk,
            pk->pl[idx]);
      }
      else {
        count_mode_sim_lin(mesh, mesh->Fk0[idx], mesh->Fk0[idx], pk,
            pk->pl[idx]);
      }
      /* Normalisation and shot noise correction. */
      for (int l = 0; l < pk->nl; l++) {
        for (int i = 0; i < pk->nbin; i++) {
          if (pk->cnt[i]) {
            pk->pl[idx][l][i] /= cat->norm[idx] * pk->cnt[i];
            /* Shot noise for multipoles. */
            pk->pl[idx][l][i] -= cat->shot[idx]
              * pk->lcnt[i + l * pk->nbin] / pk->cnt[i];
          }
          pk->pl[idx][l][i] *= 2 * pk->poles[l] + 1;
        }
      }
    }

    return;
  }

  /* Function pointer for different binning schemes. */
  void (*cntfunc) (const MESH *, const FFT_CMPLX *, const FFT_CMPLX *,
      PK *, double *);
  if (pk->log) cntfunc = count_mode_log;
  else cntfunc = count_mode_lin;

  if (ell == 0) {       /* monopole */
    if (idx < 0) {      /* cross power spectrum */
      cntfunc(mesh, mesh->Fk0[0], mesh->Fk0[1], pk, pk->xpl[0]);
      /* Normalisation */
      for (int i = 0; i < pk->nbin; i++) {
        if (pk->cnt[i]) {
          double norm = sqrt(cat->norm[0] * cat->norm[1]) * pk->cnt[i];
          pk->xpl[0][i] /= norm;
        }
      }
      if (verb) printf("  Cross power spectrum monopole computed\n");
    }
    else {              /* auto power spectrum */
      cntfunc(mesh, mesh->Fk0[idx], mesh->Fk0[idx], pk, pk->pl[idx][0]);
      /* Normalisation and shot noise subtraction. */
      for (int i = 0; i < pk->nbin; i++) {
        if (pk->cnt[i]) {
          pk->pl[idx][0][i] = pk->pl[idx][0][i] / pk->cnt[i] - cat->shot[idx];
          pk->pl[idx][0][i] /= cat->norm[idx];
        }
      }
      if (verb) {
        if (cat->num == 2)
          printf("  Auto power spectrum monopole for catalog %d computed\n",
              idx + 1);
        else printf("  Auto power spectrum monopole computed\n");
      }
    }
  }
  else {                /* multipole with ell > 0 */
    /* Find the index for the ell. */
    int pid = 0;
    for (int i = 0; i < pk->nl; i++) {
      if (pk->poles[i] == ell) {
        pid = i;
        break;
      }
    }
    if (idx < 0) {      /* cross power spectrum */
      cntfunc(mesh, mesh->Fk0[0], mesh->Fkl[1], pk, pk->xpl[pid]);
      cntfunc(mesh, mesh->Fk0[1], mesh->Fkl[0], pk, pk->xpl[pid]);
      /* Normalisation */
      for (int i = 0; i < pk->nbin; i++) {
        if (pk->cnt[i]) {
          double norm = sqrt(cat->norm[0] * cat->norm[1]) * pk->cnt[i];
          pk->xpl[pid][i] *= 2 * M_PI / norm;
        }
      }
      if (verb) printf("  Cross power spectrum with l = %d computed\n", ell);
    }
    else {              /* auto power spectrum */
      cntfunc(mesh, mesh->Fk0[idx], mesh->Fkl[idx], pk, pk->pl[idx][pid]);
      /* Normalisation */
      for (int i = 0; i < pk->nbin; i++) {
        if (pk->cnt[i]) {
          double norm = cat->norm[idx] * pk->cnt[i];
          pk->pl[idx][pid][i] *= 4 * M_PI / norm;
        }
      }
      if (verb) {
        if (cat->num == 2)
          printf("  Auto power spectrum with l = %d for catalog %d computed\n",
              ell, idx + 1);
        else printf("  Auto power spectrum with l = %d computed\n", ell);
      }
    }
  }
}


/*============================================================================*\
                    Interface for power spectrum evaluation
\*============================================================================*/

/******************************************************************************
Function `powspec`:
  Compute power spectra given the meshes.
Arguments:
  * `conf`:     the structure for storing all configurations;
  * `cat`:      the structure for storing information from catalogs;
  * `mesh`:     the structure for all meshes and FFT plans.
Return:
  Address of the structure for power spectra on success; NULL on error.
******************************************************************************/
PK *powspec(const CONF *conf, const CATA *cat, MESH *mesh) {
  printf("Evaluating power spectra ...");
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
  if (!mesh) {
    P_ERR("meshes not generated\n");
    return NULL;
  }

  /* Function pointer for density field multipole generation. */
  typedef void (*densfunc) (MESH *, const int);
  densfunc densFFT[POWSPEC_MAX_ELL + 1];
  densFFT[0] = dens_k0;
  densFFT[1] = dens_k1;
  densFFT[2] = dens_k2;
  densFFT[3] = dens_k3;
  densFFT[4] = dens_k4;
  densFFT[5] = dens_k5;
  densFFT[6] = dens_k6;

  /* Initialise the structure for power spectra. */
  PK *pk = powspec_init(conf, mesh);
  if (!pk) {
    P_ERR("failed to initialise the power spectra\n");
    return NULL;
  }

  /* Pre-compute number counts and aliasing corrections. */
  powspec_precomp(mesh, pk, conf->verbose);

  /* Generate the density field monopole, and the interlaced field. */
  densFFT[0](mesh, conf->verbose);

  if (pk->issim) {
    for (int i = 0; i < mesh->num; i++) {
      if (conf->isauto[i]) count_mode(cat, mesh, pk, i, 0, conf->verbose);
    }
    if (conf->iscross) count_mode(cat, mesh, pk, -1, 0, conf->verbose);
  }
  else {
    /* Count modes for monopole. */
    if (pk->poles[0] == 0) {
      for (int i = 0; i < mesh->num; i++) {
        if (conf->isauto[i]) count_mode(cat, mesh, pk, i, 0, conf->verbose);
      }
      if (conf->iscross) count_mode(cat, mesh, pk, -1, 0, conf->verbose);
    }

    /* Deal with the other multipoles. */
    for (int n = 1; n < pk->nl; n++) {
      int l = pk->poles[n];
      densFFT[l](mesh, conf->verbose);

      for (int i = 0; i < mesh->num; i++) {
        if (conf->isauto[i]) count_mode(cat, mesh, pk, i, l, conf->verbose);
      }
      if (conf->iscross) count_mode(cat, mesh, pk, -1, l, conf->verbose);
    }
  }

  /* Release memory for the meshes (meta data is still necessary). */
  for (int i = 0; i < mesh->num; i++) {
    FFT_FREE(mesh->Fr[i]);
    if (mesh->Frl && mesh->Frl[i]) FFT_FREE(mesh->Frl[i]);
    FFT_FREE(mesh->Fk0[i]);
    if (mesh->Fkl && mesh->Fkl[i]) FFT_FREE(mesh->Fkl[i]);
  }
  free(mesh->Fr);
  if (mesh->Frl) free(mesh->Frl);
  free(mesh->Fk0);
  if (mesh->Fkl) free(mesh->Fkl);
  free(mesh->alias);
  if (mesh->Fka) FFT_FREE(mesh->Fka);
  mesh->Fr = mesh->Frl = NULL;
  mesh->Fk0 = mesh->Fkl = NULL;
  mesh->alias = NULL;
  mesh->Fka = NULL;

  /* Rescale logarithm bins. */
  if (pk->log) {
    for (int i = 0; i < pk->nbin; i++) {
      pk->k[i] = pow(10, pk->k[i]);
      pk->km[i] = pow(10, pk->k[i]);
      pk->kedge[i] = pow(10, pk->k[i]);
    }
    pk->kedge[pk->nbin] = pow(10, pk->kedge[pk->nbin]);
  }

  printf(FMT_DONE);
  return pk;
}

