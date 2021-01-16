/*******************************************************************************
* genr_mesh.h: this file is part of the powspec program.

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

#ifndef _GENR_MESH_H_
#define _GENR_MESH_H_

#include "fftw_define.h"
#include "read_cata.h"

/*============================================================================*\
                                Data structures
\*============================================================================*/

/* Particle assignment schemes. */
typedef enum {
  POWSPEC_ASSIGN_NGP = 0,
  POWSPEC_ASSIGN_CIC = 1,
  POWSPEC_ASSIGN_TSC = 2,
  POWSPEC_ASSIGN_PCS = 3
} powspec_assign_t;

/* Name of particle assignment schemes. */
extern const char *powspec_assign_names[];

/* Data structure for storing all the meshes and FFT plans. */
typedef struct {
  int num;                      /* number of catalogs / density fields   */
  int Ng;                       /* number of grid cells per side         */
  int Ngk;                      /* number of cells on the last dimension */
  size_t Ntot;                  /* total number of grids                 */
  size_t Ncmplx;                /* number of Fourier space grids         */
  double min[3];                /* lower limit of the coordinates        */
  double max[3];                /* upper limit of the coordinates        */
  double smin[3];               /* lower boundary for the interlaced box */
  double bsize[3];              /* box size                              */
  bool issim;                   /* whether this is for a simulation box  */
  bool intlace;                 /* whether grid interlacing is enabled   */
  bool fft_init;                /* whether FFTW has been initialised     */
  powspec_assign_t assign;      /* particle assignment scheme            */
  FFT_PLAN r2c;                 /* real-to-complex FFT                   */
  FFT_PLAN c2r;                 /* complex-to-real FFT                   */
  FFT_REAL **Fr;                /* the real-space density field          */
  FFT_REAL **Frl;               /* weighted field by spherical harmonics */
  FFT_REAL *alias;              /* alias corrections on grids            */
  FFT_CMPLX **Fk0;              /* Fourier-space field for monopole      */
  FFT_CMPLX **Fkl;              /* Fourier-space field for multipoles    */
  FFT_CMPLX *Fka;               /* auxiliary Fourier-space field         */
} MESH;


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
MESH *genr_mesh(const CONF *conf, CATA *cat);

/******************************************************************************
Function `mesh_destroy`:
  Deconstruct the interface for meshes.
Arguments:
  * `mesh`:     the structure for storing the meshes and FFT plans.
******************************************************************************/
void mesh_destroy(MESH *mesh);

#endif
