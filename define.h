/*******************************************************************************
* define.h: this file is part of the powspec program.

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

#ifndef _DEFINE_H_
#define _DEFINE_H_

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/
#define SPEED_OF_LIGHT  299792.458
#ifndef M_PI
#define M_PI            0x1.921fb54442d18p+1    /* PI */
#endif
#define DEGREE_2_RAD    0x1.1df46a2529d39p-6    /* M_PI / 180 */
#define DOUBLE_EPSILON  1e-16   /* ~ machine epsilon for double numbers */
#define DOUBLE_TOL      1e-8    /* tolerance for double number comparison */

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
/* Default value for unset parameters. */
#define DEFAULT_CONF_FILE               "powspec.conf"
#define DEFAULT_FILE_FORMAT             POWSPEC_FFMT_ASCII
#define DEFAULT_FILE_SKIP               0
#define DEFAULT_FILE_CMT                '\0'
#define DEFAULT_WT_COMP                 "1"
#define DEFAULT_WT_FKP                  "1"
#define DEFAULT_CONVERT                 false
#define DEFAULT_CMVDST_ERR              1e-8
#define DEFAULT_BOX_PAD                 0.02
#define DEFAULT_PARTICLE_ASSIGN         POWSPEC_ASSIGN_TSC
#define DEFAULT_GRID_INTERLACE          false
#define DEFAULT_KMIN                    0.0
#define DEFAULT_LOG_SCALE               false
#define DEFAULT_HEADER                  true
#define DEFAULT_OVERWRITE               0
#define DEFAULT_VERBOSE                 true

#define POWSPEC_MAX_GSIZE               65536
#define POWSPEC_KMAX_UNSET_VAL          -1
#define POWSPEC_KMAX_ISSET(x)           ((x) > 0)
#define POWSPEC_MAX_MU_BIN              65536

/* Priority of parameters from different sources. */
#define POWSPEC_PRIOR_CMD               5
#define POWSPEC_PRIOR_FILE              1

/*============================================================================*\
                            Definitions for file IO
\*============================================================================*/
#define POWSPEC_PATH_SEP        '/'     /* separator for file paths     */
#define POWSPEC_FILE_CHUNK      1048576 /* chunk size for ASCII file IO */
#define POWSPEC_MAX_CHUNK       INT_MAX /* maximum allowed chunk size   */
/* Initial number of objects allocated for the catalogs.        */
#define POWSPEC_DATA_INIT_NUM   128
/* Maximum number of objects stored for each thread.            */
#define POWSPEC_DATA_THREAD_NUM 1024
/* Comment symbol for the coordinate conversion file.           */
#define POWSPEC_READ_COMMENT    '#'
/* Comment symbol for the output files.                         */
#define POWSPEC_SAVE_COMMENT    '#'

/*============================================================================*\
                            Other runtime constants
\*============================================================================*/
/* Number of redshifts for convergency tests of integrations    */
#define POWSPEC_INT_NUM_ZSP     128
/* Ceiling the adaptive box size if it is not a multiple of POWSPEC_BOX_CEIL */
#define POWSPEC_BOX_CEIL        10
#define POWSPEC_MAX_ELL         6       /* maximum ell for multipoles */

/*============================================================================*\
                     Definitions for the format of outputs
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m"               /* Red "Exit"        */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"  /* Red "FAIL"        */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL "%.10lg"             /* Output format for double parameters */

/*============================================================================*\
                          Definitions for error codes
\*============================================================================*/
#define POWSPEC_ERR_MEMORY      (-1)
#define POWSPEC_ERR_ARG         (-2)
#define POWSPEC_ERR_FILE        (-3)
#define POWSPEC_ERR_CFG         (-4)
#define POWSPEC_ERR_AST         (-5)
#define POWSPEC_ERR_ASCII       (-6)
#define POWSPEC_ERR_CONF        (-10)
#define POWSPEC_ERR_CATA        (-11)
#define POWSPEC_ERR_CNVT        (-12)
#define POWSPEC_ERR_MESH        (-13)
#define POWSPEC_ERR_PK          (-14)
#define POWSPEC_ERR_SAVE        (-15)
#define POWSPEC_ERR_UNKNOWN     (-99)

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_EXIT " " __VA_ARGS__)

/* k should vary the fastest to reduce cache miss. */
#define IDX(Ng,i,j,k)      (((size_t) (i) * (Ng) + (j)) * (Ng) + (k))

#endif

