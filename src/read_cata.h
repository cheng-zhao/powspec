/*******************************************************************************
* read_cata.h: this file is part of the powspec program.

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

#ifndef _READ_CATA_H_
#define _READ_CATA_H_

#include "load_conf.h"
#include <stdio.h>

/*============================================================================*\
              Data structure for storing information from catalogs
\*============================================================================*/

/* Data structure for the corrdinates and weight. */
typedef struct {
  double x[3];
  double w;
} DATA;

/* Data structure for the catalogs. */
typedef struct {
  int num;              /* number of data catalogs to be read       */
  DATA **data;          /* structure for the data catalogs          */
  DATA **rand;          /* structure for the random catalogs        */
  size_t *ndata;        /* number of objects in the data catalogs   */
  size_t *nrand;        /* number of objects in the random catalogs */
  double *wdata;        /* total completeness weight of the data    */
  double *wrand;        /* total completeness weight of the random  */
  double *alpha;        /* alpha value for each data-random sample  */
  double *shot;         /* shot noise term for each sample          */
  double *norm;         /* normalisation term for each sample       */
} CATA;


/*============================================================================*\
                         Interface for reading catalogs
\*============================================================================*/

/******************************************************************************
Function `read_cata`:
  Read the catalogs.
Arguments:
  * `conf`:     the structure for all configurations.
Return:
  The structure for storing information of the catalogs.
******************************************************************************/
CATA *read_cata(const CONF *conf);

/******************************************************************************
Function `cata_destroy`:
  Release memory allocated for the catalogs.
Arguments:
  * `cat`:      the structure for storing information from the catalogs.
******************************************************************************/
void cata_destroy(CATA *cat);

#endif
