/*******************************************************************************
* fftw_define.h: this file is part of the powspec program.

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

#ifndef _FFTW_DEFINE_H_
#define _FFTW_DEFINE_H_

#include <fftw3.h>

/* Macros for different FFTW precisions. */
#ifdef SINGLE_PREC
typedef float                   FFT_REAL;
typedef fftwf_complex           FFT_CMPLX;
typedef fftwf_plan              FFT_PLAN;
#define FFT_MALLOC(...)         fftwf_malloc(__VA_ARGS__)
#define FFT_INIT_OMP(...)       fftwf_init_threads(__VA_ARGS__)
#define FFT_PLAN_R2C(...)       fftwf_plan_dft_r2c_3d(__VA_ARGS__)
#define FFT_PLAN_C2R(...)       fftwf_plan_dft_c2r_3d(__VA_ARGS__)
#define FFT_PLAN_OMP(...)       fftwf_plan_with_nthreads(__VA_ARGS__)
#define FFT_EXEC_R2C(...)       fftwf_execute_dft_r2c(__VA_ARGS__)
#define FFT_EXEC_C2R(...)       fftwf_execute_dft_c2r(__VA_ARGS__)
#define FFT_FREE(...)           fftwf_free(__VA_ARGS__)
#define FFT_CLEAN(...)          fftwf_cleanup(__VA_ARGS__)
#define FFT_CLEAN_OMP(...)      fftwf_cleanup_threads(__VA_ARGS__)
#define FFT_DESTROY(...)        fftwf_destroy_plan(__VA_ARGS__)
#else
typedef double                  FFT_REAL;
typedef fftw_complex            FFT_CMPLX;
typedef fftw_plan               FFT_PLAN;
#define FFT_MALLOC(...)         fftw_malloc(__VA_ARGS__)
#define FFT_INIT_OMP(...)       fftw_init_threads(__VA_ARGS__)
#define FFT_PLAN_R2C(...)       fftw_plan_dft_r2c_3d(__VA_ARGS__)
#define FFT_PLAN_C2R(...)       fftw_plan_dft_c2r_3d(__VA_ARGS__)
#define FFT_PLAN_OMP(...)       fftw_plan_with_nthreads(__VA_ARGS__)
#define FFT_EXEC_R2C(...)       fftw_execute_dft_r2c(__VA_ARGS__)
#define FFT_EXEC_C2R(...)       fftw_execute_dft_c2r(__VA_ARGS__)
#define FFT_FREE(...)           fftw_free(__VA_ARGS__)
#define FFT_CLEAN(...)          fftw_cleanup(__VA_ARGS__)
#define FFT_CLEAN_OMP(...)      fftw_cleanup_threads(__VA_ARGS__)
#define FFT_DESTROY(...)        fftw_destroy_plan(__VA_ARGS__)
#endif

#define FFT_FLAG        FFTW_ESTIMATE   /* or FFTW_MEASURE */

#endif
