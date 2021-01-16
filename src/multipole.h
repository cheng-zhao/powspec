/*******************************************************************************
* multipole.h: this file is part of the powspec program.

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

#ifndef _MULTIPOLE_H_
#define _MULTIPOLE_H_

#include "fftw_define.h"
#include "genr_mesh.h"
#include <stdio.h>
#include <stdbool.h>

/*============================================================================*\
                                Data structures
\*============================================================================*/

typedef struct {
  bool issim;           /* indicate if the inputs are from simulations  */
  bool log;             /* indicate if using logarithm wave number bins */
  bool isauto[2];       /* indicate if auto power spectra are required  */
  bool iscross;         /* indicate if cross power spectrum is required */
  int nl;               /* number of multipoles to be computed          */
  int nbin;             /* number of k bins for the power spectra       */
  int nmu;              /* number of mu bins for the 2-D power spectra  */
  int *poles;           /* multipoles to be computed                    */
  double los[3];        /* line-of-sight vector for simulation boxes    */
  double dk;            /* width of the wave number bins                */
  double *kedge;        /* edges of the pre-defined wave number bins    */
  double *k;            /* centre of the pre-defined wave number bins   */
  double *km;           /* mean wave number of all grids inside the bin */
  size_t *cnt;          /* number of grids inside the bin               */
  double *lcnt;         /* number of weighted grids for multipoles      */
  double **pl[2];       /* auto power spectrum multipoles               */
  double **xpl;         /* cross power spectrum multipoles              */
#ifdef OMP              /* thread-private arrays for countings          */
  int nomp;             /* number of OpenMP threads                     */
  double *pcnt;         /* for summing the Fourier space densities      */
  double *plcnt;        /* for counting the number of multipole modes   */
#endif
} PK;


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
PK *powspec(const CONF *conf, const CATA *cat, MESH *mesh);

/******************************************************************************
Function `powspec_destroy`:
  Deconstruct the interface for power spectra.
Arguments:
  * `pk`:       structure for the power spectra.
******************************************************************************/
void powspec_destroy(PK *pk);

#endif
