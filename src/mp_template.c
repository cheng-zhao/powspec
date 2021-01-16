/*******************************************************************************
* mp_template.h: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

*******************************************************************************/

/*******************************************************************************
  Implementation of the power spectra multipole estimator based on spherical
  harmonic expansion decomposition of the Legender polynomials.
  ref: Hand et al. (2017): http://dx.doi.org/10.1088/1475-7516/2017/07/002
       (arXiv:1704.02357)
*******************************************************************************/

/* A template function implemented using C macros */
#ifdef POWSPEC_ELL

#define CAT_FNAME(a,b)  a##b
#define CAT(a,b)        CAT_FNAME(a,b)  /* to avoid multiple arguments */
#define POWSPEC_FUNC(x) CAT(dens_k, x)
#define SPH_HAR(x)      CAT(YlmR_l, x)

/******************************************************************************
Function `dens_k<POWSPEC_ELL>`:
  Compute the Fourier space density field multipole (ell = <POWSPEC_ELL>).
Arguments:
  * `mesh`:     the structure for meshes and FFT plans;
  * `verb`:     indicate whether to show detailed standard outputs.
******************************************************************************/
static void POWSPEC_FUNC(POWSPEC_ELL) (MESH *mesh, const int verb) {
  double fac[3];
  for (int i = 0; i < 3; i++) fac[i] = 1 / mesh->bsize[i];

  for (int n = 0; n < mesh->num; n++) {
    if (verb) {
      if (mesh->num != 2)
        printf("  Executing %d FFTs for l = %d ...",
            2 * POWSPEC_ELL + 1, POWSPEC_ELL);
      else
        printf("  Executing %d FFTs for l = %d with catalog %d ...",
            2 * POWSPEC_ELL + 1, POWSPEC_ELL, n + 1);
      fflush(stdout);
    }

    memset(mesh->Fkl[n], 0, mesh->Ncmplx * sizeof(FFT_CMPLX));
    for (int m = - POWSPEC_ELL; m <= POWSPEC_ELL; m++) {
#ifdef OMP
#pragma omp parallel for
#endif
      /* Weight the density field by spherical harmonics. */
      for (int i = 0; i < mesh->Ng; i++) {
        size_t idx_i = (size_t) i * mesh->Ng * mesh->Ng;
        double ri = (i + mesh->smin[0]) * mesh->bsize[0];
        double r1 = ri * ri;
        for (int j = 0; j < mesh->Ng; j++) {
          size_t idx_j = idx_i + j * mesh->Ng;
          if (idx_j == 0) {
            memcpy(mesh->Frl[n], mesh->Fr[n], mesh->Ng * sizeof(FFT_REAL));
            continue;
          }
          double rj = (j + mesh->smin[1]) * mesh->bsize[1];
          double r2 = r1 + rj * rj;
          for (int k = 0; k < mesh->Ng; k++) {
            size_t idx = idx_j + k;
            double rk = (k + mesh->smin[2]) * mesh->bsize[2];
            double r3 = sqrt(r2 + rk * rk);
            double r22 = sqrt(r2);

            double cost = rk / r3;
            double sint = r22 / r3;
            double cosp = ri / r22;
            double sinp = rj / r22;
            mesh->Frl[n][idx] = mesh->Fr[n][idx] *
              SPH_HAR(POWSPEC_ELL) (m, cost, sint, cosp, sinp);
          }
        }
      }

      /* Fourier transform the weighted field. */
      FFT_EXEC_R2C(mesh->r2c, mesh->Frl[n], mesh->Fka);

      /* Weight the Fourier space density field by spherical harmonics. */
#ifdef OMP
#pragma omp parallel for
#endif
      for (int i = 0; i < mesh->Ng; i++) {
        size_t idx_i = (size_t) i * mesh->Ng * mesh->Ngk;
        double ki = fac[0] * ((i <= (mesh->Ng >> 1)) ? i : i - mesh->Ng);
        double k1 = ki * ki;
        for (int j = 0; j < mesh->Ng; j++) {
          size_t idx_j = idx_i + j * mesh->Ngk;
          double kj = fac[1] * ((j <= (mesh->Ng >> 1)) ? j : j - mesh->Ng);
          double k2 = k1 + kj * kj;
          for (int k = 0; k < mesh->Ngk; k++) {
            size_t idx = idx_j + k;

            /* Skip modes that are not necessary. */
            if (!mesh->alias[idx]) continue;

            double kk = fac[2] * k;
            double k3 = sqrt(k2 + kk * kk);

            if (k2 == 0 || k3 == 0) {
              mesh->Fkl[n][idx][0] += mesh->Fka[idx][0];
              mesh->Fkl[n][idx][1] += mesh->Fka[idx][1];
            }
            else {
              double k22 = sqrt(k2);
              double cost = kk / k3;
              double sint = k22 / k3;
              double cosp = ki / k22;
              double sinp = kj / k22;
              double sh = SPH_HAR(POWSPEC_ELL) (m, cost, sint, cosp, sinp);
              mesh->Fkl[n][idx][0] += mesh->Fka[idx][0] * sh;
              mesh->Fkl[n][idx][1] += mesh->Fka[idx][1] * sh;
            }
          }
        }
      }
    }

    if (verb) {
      if (mesh->num != 2)
        printf("\r  Done with computing %d FFTs for l = %d\n",
            2 * POWSPEC_ELL + 1, POWSPEC_ELL);
      else
        printf("\r  Done with computing %d FFTs for l = %d with catalog %d\n",
          2 * POWSPEC_ELL + 1, POWSPEC_ELL, n + 1);
    }
  }
}

#endif

