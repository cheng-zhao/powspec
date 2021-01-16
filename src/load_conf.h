/*******************************************************************************
* load_conf.h: this file is part of the powspec program.

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

#ifndef _LOAD_CONF_H_
#define _LOAD_CONF_H_

#include <stdbool.h>

/*============================================================================*\
                   Data structure for storing configurations
\*============================================================================*/

typedef struct {
  char *fconf;          /* Name of the configuration file. */
  int ndata;            /* Number of data catalogs. */
  char **dfname;        /* DATA_CATALOG    */
  char **rfname;        /* RAND_CATALOG    */
  int *dftype;          /* DATA_FORMAT     */
  int *rftype;          /* RAND_FORMAT     */
  long *dskip;          /* DATA_SKIP       */
  long *rskip;          /* RAND_SKIP       */
  char *dcmt;           /* DATA_COMMENT    */
  char *rcmt;           /* RAND_COMMENT    */
  char **dfmtr;         /* DATA_FORMATTER  */
  char **rfmtr;         /* RAND_FORMATTER  */
  char **dpos;          /* DATA_POSITION   */
  char **rpos;          /* RAND_POSITION   */
  char **dwcomp;        /* DATA_WT_COMP    */
  char **rwcomp;        /* RAND_WT_COMP    */
  char **dwfkp;         /* DATA_WT_FKP     */
  char **rwfkp;         /* RAND_WT_FKP     */
  char **dnz;           /* DATA_NZ         */
  char **rnz;           /* RAND_NZ         */
  char **dsel;          /* DATA_SELECTION  */
  char **rsel;          /* RAND_SELECTION  */
  bool has_asc[2];      /* Indicate whether there is an ASCII file. */
  bool *dcnvt;          /* DATA_CONVERT    */
  bool *rcnvt;          /* RAND_CONVERT    */
  bool cnvt;            /* Indicate if coordinate conversion is required. */
  double omega_m;       /* OMEGA_M         */
  double omega_l;       /* OMEGA_LAMBDA    */
  double omega_k;       /* 1 - OMEGA_M - OMEGA_LAMBDA */
  double eos_w;         /* DE_EOS_W        */
  double ecdst;         /* CMVDST_ERR      */
  char *fcdst;          /* Z_CMVDST_CONV   */
  bool issim;           /* CUBIC_SIM       */
  double *los;          /* LINE_OF_SIGHT   */
  double *bsize;        /* BOX_SIZE        */
  double *bpad;         /* BOX_PAD         */
  int gsize;            /* GRID_SIZE       */
  int assign;           /* PARTICLE_ASSIGN */
  bool intlace;         /* GRID_INTERLACE  */
  int *poles;           /* MULTIPOLE       */
  int npole;            /* Number of multipoles to be computed. */
  double kmin;          /* KMIN            */
  double kmax;          /* KMAX            */
  bool logscale;        /* LOG_SCALE       */
  double kbin;          /* BIN_SIZE        */
  char **oauto;         /* OUTPUT_AUTO     */
  char *ocross;         /* OUTPUT_CROSS    */
  bool isauto[2];       /* Indicate whether to compute auto power spectra. */
  bool iscross;         /* Indicate whether to compute cross power spectra. */
  bool oheader;         /* OUTPUT_HEADER   */
  int ovwrite;          /* OVERWRITE       */
  bool verbose;         /* VERBOSE         */
} CONF;


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv);

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf);

#endif
