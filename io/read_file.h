/*******************************************************************************
* read_file.h: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef _READ_FILE_H_
#define _READ_FILE_H_

#include "read_cata.h"
#include <stdio.h>
#include <stdbool.h>

typedef enum {
  POWSPEC_FFMT_ASCII = 0,
  POWSPEC_FFMT_FITS = 1
} powspec_ffmt_t;

extern const char *powspec_ffmt_names[];

/*============================================================================*\
                           Interfaces for ASCII files
\*============================================================================*/

/******************************************************************************
Function `read_ascii_simple`:
  Read the first two columns of an ASCII file as double arrays.
Arguments:
  * `fname`:    filename of the input catalog;
  * `x`:        array for the first column;
  * `y`:        array for the second column;
  * `num`:      number of lines read successfully;
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_simple(const char *fname, double **x, double **y, size_t *num,
    const int verb);

/******************************************************************************
Function `read_ascii_data`:
  Read an ASCII file for the positions and weights.
Arguments:
  * `fname`:    filename of the input catalog;
  * `skip`:     number of lines to be skipped before reading positions;
  * `comment`:  character indicating the beginning of a comment line;
  * `fmtr`:     formatter string for `sscanf`;
  * `pos`:      columns of the positions;
  * `wc`:       completeness weight;
  * `wfkp`:     FKP weight;
  * `nz`:       radial number density distribution;
  * `sel`:      data selection criteria;
  * `wt`:       indicate whether to compute weights;
  * `data`:     address of the structure for storing positions;
  * `num`:      number of lines read successfully;
  * `sumw`:     sum of completeness weights;
  * `I12`:      the I12 factor;
  * `I22`:      the I22 factor;
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_data(const char *fname, const size_t skip, const char comment,
    const char *fmtr, char *const *pos, const char *wc, const char *wfkp,
    const char *nz, const char *sel, const bool wt, DATA **data, size_t *num,
    double *sumw, double *I12, double *I22, const int verb);

#endif
