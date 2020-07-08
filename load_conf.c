/*******************************************************************************
* load_conf.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "libcfg.h"
#include "read_file.h"
#include "genr_mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
/* Print the warning and error messages. */
#define P_CFG_WRN(cfg)  cfg_pwarn(cfg, stderr, FMT_WARN);
#define P_CFG_ERR(cfg)  {                                       \
  cfg_perror(cfg, stderr, FMT_ERR);                             \
  cfg_destroy(cfg);                                             \
  return NULL;                                                  \
}


/*============================================================================*\
                    Functions called via command line flags
\*============================================================================*/

/******************************************************************************
Function `usage`:
  Print the usage of command line options.
******************************************************************************/
static void usage(void *args) {
  printf("Usage: POWSPEC [OPTION]\n\
Compute the auto and cross power spectra from catalogs.\n\
  -h, --help\n\
        Display this message and exit\n\
  -t, --template\n\
        Print a template configuration file to the standard output and exit\n\
  -c, --conf            " FMT_KEY(CONFIG_FILE) "     String\n\
        Specify the configuration file (default: `%s')\n\
  -d, --data            " FMT_KEY(DATA_CATALOG) "    String array\n\
        Specify the input data catalogs\n\
  -r, --rand            " FMT_KEY(RAND_CATALOG) "    String array\n\
        Specify the input random catalogs\n\
  -f, --data-format     " FMT_KEY(DATA_FORMAT) "     Integer array\n\
  -F, --rand-format     " FMT_KEY(RAND_FORMAT) "     Integer array\n\
        Format of the input catalogs\n\
      --data-skip       " FMT_KEY(DATA_SKIP) "       Long integer array\n\
      --rand-skip       " FMT_KEY(RAND_SKIP) "       Long integer array\n\
        Number of lines to be skipped for the ASCII format input catalogs\n\
      --data-comment    " FMT_KEY(DATA_COMMENT) "    Character array\n\
      --rand-comment    " FMT_KEY(RAND_COMMENT) "    Character array\n\
        Comment symbols for the ASCII format input catalogs\n\
      --data-formatter  " FMT_KEY(DATA_FORMATTER) "  String array\n\
      --rand-formatter  " FMT_KEY(RAND_FORMATTER) "  String array\n\
        Formatter for columns of ASCII format input catalogs\n\
  -p, --data-pos        " FMT_KEY(DATA_POSITION) "   String array\n\
  -P, --rand-pos        " FMT_KEY(RAND_POSITION) "   String array\n\
        Column indicator for the positions of the input samples\n\
      --data-wt-comp    " FMT_KEY(DATA_WT_COMP) "    String arary\n\
      --rand-wt-comp    " FMT_KEY(RAND_WT_COMP) "    String array\n\
        Column indicator for completeness weights of the input samples\n\
      --data-wt-fkp     " FMT_KEY(DATA_WT_FKP) "     String array\n\
      --rand-wt-fkp     " FMT_KEY(RAND_WT_FKP) "     String array\n\
        Column indicator for FKP weights of the input samples\n\
      --data-nz         " FMT_KEY(DATA_NZ) "         String array\n\
      --rand-nz         " FMT_KEY(RAND_NZ) "         String array\n\
        Column indicator for comoving number densities of the input samples\n\
      --data-select     " FMT_KEY(DATA_SELECTION) "  String array\n\
      --rand-select     " FMT_KEY(RAND_SELECTION) "  String array\n\
        Specify the criteria for sample selections\n\
      --data-convert    " FMT_KEY(DATA_CONVERT) "    Boolean array\n\
      --rand-convert    " FMT_KEY(RAND_CONVERT) "    Boolean array\n\
        Indicate whether coordinate conversions are needed\n\
  -m, --omega-m         " FMT_KEY(OMEGA_M) "         Double\n\
      --omega-l         " FMT_KEY(OMEGA_LAMBDA) "    Double\n\
      --eos-w           " FMT_KEY(DE_EOS_W) "        Double\n\
        Cosmological parameters for coordinate conversion\n\
      --cmvdst-err      " FMT_KEY(CMVDST_ERR) "      Double\n\
        Tolerance of the integration for coordinate conversion\n\
      --cmvdst-file     " FMT_KEY(Z_CMVDST_CNVT) "   String\n\
        Specify the table file for redshift to comoving distance conversion\n\
  -s, --sim             " FMT_KEY(CUBIC_SIM) "       Boolean\n\
        Indicate whether the samples are from a simulation box\n\
      --line-of-sight   " FMT_KEY(LINE_OF_SIGHT) "   Double array\n\
        The unit line-of-sight vector\n\
  -B, --box-size        " FMT_KEY(BOX_SIZE) "        Double array\n\
        Side lengths of the box for the samples to be placed in\n\
      --box-pad         " FMT_KEY(BOX_PAD) "         Double array\n\
        Padding fraction for adaptive box size\n\
  -G, --grid-size       " FMT_KEY(GRID_SIZE) "       Integer\n\
        Number of grids per box side for the density field\n\
  -n, --part-assign     " FMT_KEY(PARTICLE_ASSIGN) " Integer\n\
        Specify the particle assignment scheme\n\
  -i, --interlace       " FMT_KEY(GRID_INTERLACE) "  Boolean\n\
        Indicate whether to generate an interlaced grid\n\
  -l, --pole            " FMT_KEY(MULTIPOLE) "       Integer array\n\
        Specify the order of the multipoles to be evaluated\n\
      --log-scale       " FMT_KEY(LOG_SCALE) "       Boolean\n\
        Indicate whether to use logarithm wave number bins\n\
  -k, --kmin            " FMT_KEY(KMIN) "            Double\n\
  -K, --kmax            " FMT_KEY(KMAX) "            Double\n\
        Specify the range of the output wave number bins\n\
  -b, --bin-size        " FMT_KEY(BIN_SIZE) "        Double\n\
        Specify the size of the output wave number bins\n\
  -a, --auto            " FMT_KEY(OUTPUT_AUTO) "     String array\n\
        Specify the output files for auto power spectrum multipoles\n\
  -x, --cross           " FMT_KEY(OUTPUT_CROSS) "    String\n\
        Specify the output file for cross power spectrum multipoles\n\
      --output-header   " FMT_KEY(OUTPUT_HEADER) "   Boolean\n\
        Indicate whether to save extra information as header of the outputs\n\
  -w, --overwrite       " FMT_KEY(OVERWRITE) "       Integer\n\
        Indicate whether to overwrite existing output files\n\
  -v, --verbose         " FMT_KEY(VERBOSE) "         Boolean\n\
        Indicate whether to display detailed standard outputs\n\
Consult the -t option for more information on the parameters.\n\
Github repository: https://github.com/cheng-zhao/powspec.\n\
Licence: GPLv3.\n",
    DEFAULT_CONF_FILE);
  exit(0);
}

/******************************************************************************
Function `conf_template`:
  Print a template configuration file.
******************************************************************************/
static void conf_template(void *args) {
  printf("# Configuration file for powspec (default: `%s').\n\
# Format: keyword = value # comment\n\
#     or: keyword = [element1, element2]\n\
#    see: https://github.com/cheng-zhao/libcfg for details.\n\
# Some of the entries allow expressions, see\n\
#         https://github.com/cheng-zhao/libast for details.\n\
# NOTE that command line options have priority over this file.\n\
# Unnecessary entries can be left unset.\n\
\n\
####################################################\n\
#  Specifications of the data and random catalogs  #\n\
####################################################\n\
\n\
DATA_CATALOG    = \n\
RAND_CATALOG    = \n\
    # Filename of the input data/random catalogs.\n\
    # String or 2-element string array (for cross power spectra).\n\
DATA_FORMAT     = \n\
RAND_FORMAT     = \n\
    # Format of the input catalogs (unset: %d).\n\
    # Integer, same dimension as `DATA_CATALOG`.\n\
    # Allowed values are:\n\
    # * %d: ASCII text file;\n\
    # * %d: FITS table.\n\
DATA_SKIP       = \n\
RAND_SKIP       = \n\
    # Number of lines to be skipped for ASCII format input files (unset: %d).\n\
    # Long integer, same dimension as `DATA_CATALOG`.\n\
    # If only one of catalogs is in ASCII format, set an arbitrary value\n\
    #   for the other catalog.\n\
DATA_COMMENT    = \n\
RAND_COMMENT    = \n\
    # Characters for indicating comment lines for ASCII files (unset: '').\n\
    # Character, same dimension as `DATA_CATALOG`. '' for disabling comments.\n\
DATA_FORMATTER  = \n\
RAND_FORMATTER  = \n\
    # C99-style formatter for reading columns of ASCII format input files.\n\
    # String, same dimension as `DATA_CATALOG` (e.g. \"%%d %%ld %%f %%lf %%s\").\n\
    # If a column is suppressed by *, it is not counted for the column number.\n\
    #   e.g., for \"%%d %%*s %%f\", the float number corresponds to column $2.\n\
    # See https://en.cppreference.com/w/c/io/fscanf for details on the format.\n\
DATA_POSITION   = \n\
RAND_POSITION   = \n\
    # 3-D coordinates, in the order of [x,y,z] or [RA,Dec,redshift].\n\
    # 3- or 6-element string array, depending on the number of input catalogs.\n\
    # They can be column indicator or expressions (e.g. \"${RA}\" or \"$1 %% 100\").\n\
    # Allowed values enclosed by ${}:\n\
    # * string: column name of a FITS file;\n\
    # * long integer: column number of an ASCII file (starting from 1).\n\
DATA_WT_COMP    = \n\
RAND_WT_COMP    = \n\
    # Completeness weights of data/random (unset: %s).\n\
    # Column indicator or expression, same dimension as `DATA_CATALOG`.\n\
DATA_WT_FKP     = \n\
RAND_WT_FKP     = \n\
    # FKP weights of data/random (unset: %s).\n\
    # Column indicator or expression, same dimension as `DATA_CATALOG`.\n\
DATA_NZ         = \n\
RAND_NZ         = \n\
    # Radial number density of data/random.\n\
    # Column indicator or expression, same dimension as `DATA_CATALOG`.\n\
    # `RAND_NZ` is used by default, if it is unset then `DATA_NZ` is used.\n\
DATA_SELECTION  = \n\
RAND_SELECTION  = \n\
    # Selection criteria for the data/random (unset: no selection).\n\
    # Logical expression, same dimension as `DATA_CATALOG` (e.g. \"$3 > 0.5\").\n\
DATA_CONVERT    = \n\
RAND_CONVERT    = \n\
    # Boolean option, same dimension as `DATA_CATALOG` (unset: %c).\n\
    # True (T) for converting the coordinates from [RA,Dec,redshift] to\n\
    #   the comoving [x,y,z], given the fiducial cosmology.\n\
\n\
##################################################\n\
#  Fiducial cosmology for coordinate conversion  #\n\
##################################################\n\
\n\
OMEGA_M         = \n\
    # Density parameter of matter at z = 0.\n\
    # Double-precision number.\n\
OMEGA_LAMBDA    = \n\
    # Density parameter of Lambda at z = 0 (unset: 1 - OMEGA_M).\n\
    # Double-precision number.\n\
DE_EOS_W        = \n\
    # Dark energy equation of state: w (unset: -1).\n\
    # Double-precision number.\n\
CMVDST_ERR      = \n\
    # Error for comoving distance evaluation (unset: %g).\n\
    # Double-precision number.\n\
Z_CMVDST_CNVT   = \n\
    # Filename of a table for redshift to comoving distance conversion.\n\
    # It must be a text file with two columns: (redshift, comoving distance).\n\
    # If this file is set, the cosmological parameters above are omitted.\n\
    # Lines start with '%c' are omitted.\n\
\n\
##################################################\n\
#  Configurations for power spectra evaluation   #\n\
##################################################\n\
\n\
CUBIC_SIM       = \n\
    # Indicate whether the input catalogs are from cubic simulation boxes.\n\
    # Boolean option, true for omitting the weights and random catalogs.\n\
LINE_OF_SIGHT   = \n\
    # Unit line-of-sight vector for cubic simulation boxes (unset: [0,0,1]).\n\
    # 3-element double array.\n\
BOX_SIZE        = \n\
    # Side length of the box that catalogs are placed in.\n\
    # Double-precision number or 3-element double array.\n\
    # Coordinates from simulation boxes should be in [0, `BOX_SIZE`).\n\
    # For non-simulation case the catalog is placed at the centre of the box.\n\
    # It is mandatory if `CUBIC_SIM` is true.\n\
BOX_PAD         = \n\
    # Fraction of the extent of catalogs to be padded with 0 (unset: %g).\n\
    # Double-precision number or 3-element double array.\n\
    # It is for determing the box size automatically, if `BOX_SIZE` is unset.\n\
GRID_SIZE       = \n\
    # Number of grid cells per box side for the density fields and FFT.\n\
    # Integer, preferably the power of 2.\n\
PARTICLE_ASSIGN = \n\
    # Scheme for assigning particles to grids (unset: %d).\n\
    # Integer, allowed values are\n\
    # * %d: Nearest-Grid-Point (NGP);\n\
    # * %d: Could-In-Cell (CIC);\n\
    # * %d: Triangular Shaped Cloud (TSC);\n\
    # * %d: Piecewise Cubic Spline (PCS).\n\
GRID_INTERLACE  = \n\
    # Boolean option, indicate whether to use interlaced grids (unset: %c).\n\
\n\
#############################\n\
#  Settings for the output  #\n\
#############################\n\
\n\
MULTIPOLE       = \n\
    # Legendre multipoles to be evaluated, up to ell = %d.\n\
    # Integer or integer array (e.g. \"[0, 2, 4]\")\n\
LOG_SCALE       = \n\
    # Indicator for linear or logarithm wave number bins (unset: %c).\n\
    # Boolean option, true for logarithm bins, false for linear bins.\n\
KMIN            = \n\
    # Lower boundary of the first wave number bin (unset: %g).\n\
    # Double-precision number.\n\
KMAX            = \n\
    # Upper boundary of the last wave number bin (unset: Nyquist frequency).\n\
    # Double-precision number.\n\
    # It is rounded to the closest bin edge defined by `KMIN` and `BIN_SIZE`,\n\
    #   and below the Nyquist frequency.\n\
BIN_SIZE        = \n\
    # Width of each wave number bin.\n\
    # Base-10 logarithm of the ratio between two bins for logarithm scale.\n\
OUTPUT_AUTO     = \n\
    # Name of the output files for auto power spectra.\n\
    # String or 2-element string array. Unset or \"\" for only cross correlations.\n\
OUTPUT_CROSS    = \n\
    # String, name of the output file for cross power spectrum.\n\
OUTPUT_HEADER   = \n\
    # Boolean option, indicate whether to write extra information (unset: %c)\n\
    # If true, write metadata of the catalogs and meshes, as well as shot noise\n\
    #   and normalisation terms to the header of the output files.\n\
OVERWRITE       = \n\
    # Integer, indicate whether to overwrite existing files (unset: %d).\n\
    # Allowed values are:\n\
    # * 0: quit the program when an output file exist;\n\
    # * positive: force overwriting output files whenever possible;\n\
    # * negative: notify at most this number of times for existing files.\n\
VERBOSE         = \n\
    # Boolean option, indicate whether show detailed outputs (unset: %c).\n",
  DEFAULT_CONF_FILE, DEFAULT_FILE_FORMAT, POWSPEC_FFMT_ASCII,
  POWSPEC_FFMT_FITS, DEFAULT_FILE_SKIP, DEFAULT_WT_COMP, DEFAULT_WT_FKP,
  DEFAULT_CONVERT ? 'T' : 'F', DEFAULT_CMVDST_ERR, POWSPEC_READ_COMMENT,
  DEFAULT_BOX_PAD, DEFAULT_PARTICLE_ASSIGN, POWSPEC_ASSIGN_NGP,
  POWSPEC_ASSIGN_CIC, POWSPEC_ASSIGN_TSC, POWSPEC_ASSIGN_PCS,
  DEFAULT_GRID_INTERLACE ? 'T' : 'F', POWSPEC_MAX_ELL,
  DEFAULT_LOG_SCALE ? 'T' : 'F', DEFAULT_KMIN, DEFAULT_HEADER ? 'T' : 'F',
  DEFAULT_OVERWRITE, DEFAULT_VERBOSE ? 'T' : 'F');
  exit(0);
}


/*============================================================================*\
                      Function for reading configurations
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing configurations.
Return:
  Address of the structure.
******************************************************************************/
static CONF *conf_init(void) {
  CONF *conf = calloc(1, sizeof *conf);
  if (!conf) return NULL;
  conf->fconf = NULL;
  conf->dfname = conf->rfname = conf->dfmtr = conf->rfmtr = NULL;
  conf->dpos = conf->rpos = conf->dwcomp = conf->rwcomp;
  conf->dwfkp = conf->rwfkp = conf->dnz = conf->rnz = NULL;
  conf->dsel = conf->rsel = conf->oauto = NULL;
  conf->dftype = conf->rftype = conf->poles = NULL;
  conf->dskip = conf->rskip = NULL;
  conf->dcmt = conf->rcmt = conf->fcdst = conf->ocross = NULL;
  conf->dcnvt = conf->rcnvt = NULL;
  conf->los = conf->bsize = conf->bpad = NULL;
  return conf;
}

/******************************************************************************
Function `conf_read`:
  Read configurations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  Interface of libcfg.
******************************************************************************/
static cfg_t *conf_read(CONF *conf, const int argc, char *const *argv) {
  if (!conf) {
    P_ERR("the structure for configurations is not initialised.\n");
    return NULL;
  }
  cfg_t *cfg = cfg_init();
  if (!cfg) P_CFG_ERR(cfg);

  /* Functions to be called via command line flags. */
  const int nfunc = 2;
  const cfg_func_t funcs[] = {
    {   'h',        "help",             usage,          NULL},
    {   't',    "template",     conf_template,          NULL}
  };

  /* Configuration parameters. */
  const int npar = 45;
  const cfg_param_t params[] = {
    {'c', "conf"          , "CONFIG_FILE"    , CFG_DTYPE_STR , &conf->fconf},
    {'d', "data"          , "DATA_CATALOG"   , CFG_ARRAY_STR , &conf->dfname},
    {'r', "rand"          , "RAND_CATALOG"   , CFG_ARRAY_STR , &conf->rfname},
    {'f', "data-format"   , "DATA_FORMAT"    , CFG_ARRAY_INT , &conf->dftype},
    {'F', "rand-format"   , "RAND_FORMAT"    , CFG_ARRAY_INT , &conf->rftype},
    { 0 , "data-skip"     , "DATA_SKIP"      , CFG_ARRAY_LONG, &conf->dskip},
    { 0 , "rand-skip"     , "RAND_SKIP"      , CFG_ARRAY_LONG, &conf->dskip},
    { 0 , "data-comment"  , "DATA_COMMENT"   , CFG_ARRAY_CHAR, &conf->dcmt},
    { 0 , "rand-comment"  , "RAND_COMMENT"   , CFG_ARRAY_CHAR, &conf->rcmt},
    { 0 , "data-formatter", "DATA_FORMATTER" , CFG_ARRAY_STR , &conf->dfmtr},
    { 0 , "rand-formatter", "RAND_FORMATTER" , CFG_ARRAY_STR , &conf->rfmtr},
    {'p', "data-pos"      , "DATA_POSITION"  , CFG_ARRAY_STR , &conf->dpos},
    {'P', "rand-pos"      , "RAND_POSITION"  , CFG_ARRAY_STR , &conf->rpos},
    { 0 , "data-wt-comp"  , "DATA_WT_COMP"   , CFG_ARRAY_STR , &conf->dwcomp},
    { 0 , "rand-wt-comp"  , "RAND_WT_COMP"   , CFG_ARRAY_STR , &conf->rwcomp},
    { 0 , "data-wt-fkp"   , "DATA_WT_FKP"    , CFG_ARRAY_STR , &conf->dwfkp},
    { 0 , "rand-wt-fkp"   , "RAND_WT_FKP"    , CFG_ARRAY_STR , &conf->rwfkp},
    { 0 , "data-nz"       , "DATA_NZ"        , CFG_ARRAY_STR , &conf->dnz},
    { 0 , "rand-nz"       , "RAND_NZ"        , CFG_ARRAY_STR , &conf->rnz},
    { 0 , "data-select"   , "DATA_SELECTION" , CFG_ARRAY_STR , &conf->dsel},
    { 0 , "rand-select"   , "RAND_SELECTION" , CFG_ARRAY_STR , &conf->rsel},
    { 0 , "data-convert"  , "DATA_CONVERT"   , CFG_ARRAY_BOOL, &conf->dcnvt},
    { 0 , "rand-convert"  , "RAND_CONVERT"   , CFG_ARRAY_BOOL, &conf->rcnvt},
    {'m', "omega-m"       , "OMEGA_M"        , CFG_DTYPE_DBL , &conf->omega_m},
    { 0 , "omega-l"       , "OMEGA_LAMBDA"   , CFG_DTYPE_DBL , &conf->omega_l},
    { 0 , "eos-w"         , "DE_EOS_W"       , CFG_DTYPE_DBL , &conf->eos_w},
    { 0 , "cmvdst-err"    , "CMVDST_ERR"     , CFG_DTYPE_DBL , &conf->ecdst},
    { 0 , "cmvdst-file"   , "Z_CMVDST_CNVT"  , CFG_DTYPE_STR , &conf->fcdst},
    {'s', "sim"           , "CUBIC_SIM"      , CFG_DTYPE_BOOL, &conf->issim},
    { 0 , "line-of-sight" , "LINE_OF_SIGHT"  , CFG_ARRAY_DBL , &conf->los},
    {'B', "box-size"      , "BOX_SIZE"       , CFG_ARRAY_DBL , &conf->bsize},
    { 0 , "box-pad"       , "BOX_PAD"        , CFG_ARRAY_DBL , &conf->bpad},
    {'G', "grid-size"     , "GRID_SIZE"      , CFG_DTYPE_INT , &conf->gsize},
    {'n', "part-assign"   , "PARTICLE_ASSIGN", CFG_DTYPE_INT , &conf->assign},
    {'i', "interlace"     , "GRID_INTERLACE" , CFG_DTYPE_BOOL, &conf->intlace},
    {'l', "pole"          , "MULTIPOLE"      , CFG_ARRAY_INT , &conf->poles},
    { 0 , "log-scale"     , "LOG_SCALE"      , CFG_DTYPE_BOOL, &conf->logscale},
    {'k', "kmin"          , "KMIN"           , CFG_DTYPE_DBL , &conf->kmin},
    {'K', "kmax"          , "KMAX"           , CFG_DTYPE_DBL , &conf->kmax},
    {'b', "bin-size"      , "BIN_SIZE"       , CFG_DTYPE_DBL , &conf->kbin},
    {'a', "auto"          , "OUTPUT_AUTO"    , CFG_ARRAY_STR , &conf->oauto},
    {'x', "cross"         , "OUTPUT_CROSS"   , CFG_DTYPE_STR , &conf->ocross},
    { 0 , "output-header" , "OUTPUT_HEADER"  , CFG_DTYPE_BOOL, &conf->oheader},
    {'w', "overwrite"     , "OVERWRITE"      , CFG_DTYPE_INT , &conf->ovwrite},
    {'v', "verbose"       , "VERBOSE"        , CFG_DTYPE_BOOL, &conf->verbose}
  };

  /* Register functions and parameters. */
  if (cfg_set_funcs(cfg, funcs, nfunc)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);
  if (cfg_set_params(cfg, params, npar)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read configurations from command line options. */
  int optidx;
  if (cfg_read_opts(cfg, argc, argv, POWSPEC_PRIOR_CMD, &optidx))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read parameters from configuration file. */
  if (!cfg_is_set(cfg, &conf->fconf)) conf->fconf = DEFAULT_CONF_FILE;
  if (cfg_read_file(cfg, conf->fconf, POWSPEC_PRIOR_FILE)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  return cfg;
}


/*============================================================================*\
                      Functions for parameter verification
\*============================================================================*/

/******************************************************************************
Function `check_input`:
  Check whether an input file can be read.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_input(const char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR("the input " FMT_KEY(%s) " is not set.\n", key);
    return POWSPEC_ERR_CFG;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot access " FMT_KEY(%s) ": `%s'.\n", key, fname);
    return POWSPEC_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `check_output`:
  Check whether an output file can be written.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file;
  * `ovwrite`:  option for overwriting exisiting files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_output(char *fname, const char *key, const int ovwrite) {
  if (!fname || *fname == '\0') {
    P_ERR("the output " FMT_KEY(%s) " is not set.\n", key);
    return POWSPEC_ERR_CFG;
  }

  /* Check if the file exists. */
  if (!access(fname, F_OK)) {
    /* not overwriting */
    if (ovwrite == 0) {
      P_ERR("the output " FMT_KEY(%s) " exists: `%s'.\n", key, fname);
      return POWSPEC_ERR_FILE;
    }
    /* force overwriting */
    else if (ovwrite > 0) {
      P_WRN("the output " FMT_KEY(%s) " will be overwritten: `%s'.\n",
          key, fname);
    }
    /* ask for decision */
    else {
      P_WRN("the output " FMT_KEY(%s) " exists: `%s'.\n", key, fname);
      char confirm = 0;
      for (int i = 0; i != ovwrite; i--) {
        fprintf(stderr, "Are you going to overwrite it? (y/n): ");
        if (scanf("%c", &confirm) != 1) continue;
        int c;
        while((c = getchar()) != '\n' && c != EOF) continue;
        if (confirm == 'n') {
          P_ERR("cannot write to the file.\n");
          return POWSPEC_ERR_FILE;
        }
        else if (confirm == 'y') break;
      }
      if (confirm != 'y') {
        P_ERR("too many failed inputs.\n");
        return POWSPEC_ERR_FILE;
      }
    }

    /* Check file permission for overwriting. */
    if (access(fname, W_OK)) {
      P_ERR("cannot write to file `%s'.\n", fname);
      return POWSPEC_ERR_FILE;
    }
  }
  /* Check if the path permission. */
  else {
    char *end;
    if ((end = strrchr(fname, POWSPEC_PATH_SEP)) != NULL) {
      *end = '\0';
      if (access(fname, X_OK)) {
        P_ERR("cannot access the directory `%s'.\n", fname);
        return POWSPEC_ERR_FILE;
      }
      *end = POWSPEC_PATH_SEP;
    }
  }
  return 0;
}

/******************************************************************************
Function `check_file_fmt`:
  Verify format settings for input files.
Arguments:
  * `cfg`:      interface of libcfg;
  * `num`:      expected number of files;
  * `name`:     namespace of the input files ("DATA" or "RAND");
  * `issim`:    indicate whether the files are from cubic simulation boxes;
  * `fname`:    array of filenames;
  * `has_asc`:  indicate whether there is an ASCII file;
  * ...:        format settings for input files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_file_fmt(const cfg_t *cfg, const int num, const char *name,
    const bool issim, char ***fname, bool *has_asc, int **type, long **skip,
    char **cmt, char ***fmtr, char ***pos, char ***wcomp, char ***wfkp,
    char ***nz, char ***sel, bool **cnvt) {
  int i, n;
  /* Check CATALOG. */
  if (!(n = cfg_get_size(cfg, fname))) {
    P_ERR(FMT_KEY(%s_CATALOG) " is not set.\n", name);
    return POWSPEC_ERR_CFG;
  }
  if (n < num) {
    P_ERR("too few elements of " FMT_KEY(%s_CATALOG) " (" FMT_KEY(OUTPUT_AUTO)
        " or " FMT_KEY(OUTPUT_CROSS) " may require more catalogs).\n", name);
    return POWSPEC_ERR_CFG;
  }
  else if (n > num) {
    P_WRN("omitting the following " FMT_KEY(%s_CATALOG) ":\n", name);
    for (i = num; i < n; i++) fprintf(stderr, "%s\n", (*fname)[i]);
  }

  /* Check FORMAT. */
  if (!(n = cfg_get_size(cfg, type))) {
    P_WRN(FMT_KEY(%s_FORMAT) " is not set.\nUse the default value: %d (%s)\n",
        name, DEFAULT_FILE_FORMAT, powspec_ffmt_names[DEFAULT_FILE_FORMAT]);
  }
  else if (n < num) {
    P_ERR("too few elements of " FMT_KEY(%s_FORMAT) ".\n", name);
    return POWSPEC_ERR_CFG;
  }
  else if (n > num) {
    P_WRN("omitting the following " FMT_KEY(%s_FORMAT) ":", name);
    for (i = num; i < n; i++) fprintf(stderr, " %d", (*type)[i]);
    fprintf(stderr, "\n");
  }

  *has_asc = false;
  for (int j = 0; j < num; j++) {
    int t = (*type) ? (*type)[j] : DEFAULT_FILE_FORMAT;
    switch(t) {
      case POWSPEC_FFMT_ASCII:
        *has_asc = true;
        break;
      case POWSPEC_FFMT_FITS:
#ifdef WITH_CFITSIO
        break;
#else
        P_ERR("FITS format is not enabled.\n\
Please re-compile the code with option -DWITH_CFITSIO.\n");
        return POWSPEC_ERR_CFG;
#endif
      default:
        P_ERR("invalid " FMT_KEY(%s_FORMAT) ": %d\n", name, t);
        return POWSPEC_ERR_CFG;
    }
  }

  /* Check format settings for ASCII catalogs. */
  if (*has_asc) {
    /* Check SKIP. */
    if ((n = cfg_get_size(cfg, skip))) {
      if (n < num) {
        P_ERR("too few elements of " FMT_KEY(%s_SKIP) ".\n", name);
        return POWSPEC_ERR_CFG;
      }
      else if (n > num) {
        P_WRN("omitting the following " FMT_KEY(%s_SKIP) ":", name);
        for (i = num; i < n; i++) fprintf(stderr, " %ld", (*skip)[i]);
        fprintf(stderr, "\n");
      }
      for (i = 0; i < num; i++) {
        if ((*skip)[i] < 0) {
          P_ERR(FMT_KEY(%s_SKIP) " cannot be negative.\n", name);
          return POWSPEC_ERR_CFG;
        }
      }
    }
    /* Check COMMENT. */
    if ((n = cfg_get_size(cfg, cmt))) {
      if (n < num) {
        P_ERR("too few elements of " FMT_KEY(%s_COMMENT) ".\n", name);
        return POWSPEC_ERR_CFG;
      }
      else if (n > num) {
        P_WRN("omitting the following " FMT_KEY(%s_COMMENT) ":", name);
        for (i = num; i < n; i++) fprintf(stderr, " '%c'", (*cmt)[i]);
        fprintf(stderr, "\n");
      }
    }
    /* Check FORMATTER. */
    if (!(n = cfg_get_size(cfg, fmtr))) {
      P_ERR(FMT_KEY(%s_FORMATTER) " is not set.\n", name);
      return POWSPEC_ERR_CFG;
    }
    else if (n < num) {
      P_ERR("too few elements of " FMT_KEY(%s_FORMATTER) ".\n", name);
      return POWSPEC_ERR_CFG;
    }
    else if (n > num) {
      P_WRN("omitting the following " FMT_KEY(%s_FORMATTER) ":\n", name);
      for (i = num; i < n; i++) fprintf(stderr, "%s\n", (*fmtr)[i]);
    }
  }

  /* Check common settings for all file types. */
  /* Check SELECTION. */
  if ((n = cfg_get_size(cfg, sel))) {
    if (n < num) {
      P_ERR("too few elements of " FMT_KEY(%s_SELECTION) ".\n", name);
      return POWSPEC_ERR_CFG;
    }
    else if (n > num) {
      P_WRN("omitting the following " FMT_KEY(%s_SELECTION) ":\n", name);
      for (i = num; i < n; i++) fprintf(stderr, "%s\n", (*sel)[i]);
    }
  }

  /* Check WT_COMP. */
  if ((n = cfg_get_size(cfg, wcomp))) {
    if (n < num) {
      P_ERR("too few elements of " FMT_KEY(%s_WT_COMP) ".\n", name);
      return POWSPEC_ERR_CFG;
    }
    else if (n > num) {
      P_WRN("omitting the following " FMT_KEY(%s_WT_COMP) ":\n", name);
      for (i = num; i < n; i++) fprintf(stderr, "%s\n", (*wcomp)[i]);
    }
  }

  if (!issim) {
    /* Check WT_FKP. */
    if ((n = cfg_get_size(cfg, wfkp))) {
      if (n < num) {
        P_ERR("too few elements of " FMT_KEY(%s_WT_FKP) ".\n", name);
        return POWSPEC_ERR_CFG;
      }
      else if (n > num) {
        P_WRN("omitting the following " FMT_KEY(%s_WT_FKP) ":\n", name);
        for (i = num; i < n; i++) fprintf(stderr, "%s\n", (*wfkp)[i]);
      }
    }
    /* Check NZ. */
    if ((n = cfg_get_size(cfg, nz))) {
      if (n < num) {
        P_ERR("too few elements of " FMT_KEY(%s_NZ) ".\n", name);
        return POWSPEC_ERR_CFG;
      }
      else if (n > num) {
        P_WRN("omitting the following " FMT_KEY(%s_NZ) ":\n", name);
        for (i = num; i < n; i++) fprintf(stderr, "%s\n", (*nz)[i]);
      }
    }
    /* Check CONVERT. */
    if ((n = cfg_get_size(cfg, cnvt))) {
      if (n < num) {
        P_ERR("too few elements of " FMT_KEY(%s_CONVERT) ".\n", name);
        return POWSPEC_ERR_CFG;
      }
      else if (n > num) {
        P_WRN("omitting the following " FMT_KEY(%s_CONVERT) ":", name);
        for (i = num; i < n; i++)
          fprintf(stderr, " %c", (*cnvt)[i] ? 'T' : 'F');
        fprintf(stderr, "\n");
      }
    }
  }

  /* Check POSITION. */
  if (!(n = cfg_get_size(cfg, pos))) {
    P_ERR(FMT_KEY(%s_POSITION) " is not set.\n", name);
    return POWSPEC_ERR_CFG;
  }
  else if (n < num * 3) {
    P_ERR("too few elements of " FMT_KEY(%s_POSITION) ".\n", name);
    return POWSPEC_ERR_CFG;
  }
  else if (n > num * 3) {
    P_WRN("omitting the following " FMT_KEY(%s_POSITION) ":\n", name);
    for (i = num * 3; i < n; i++) fprintf(stderr, "%s\n", (*pos)[i]);
  }
  return 0;
}

/******************************************************************************
Function `check_cosmo`:
  Verify cosmological parameters for coordinate conversion.
Arguments:
  * `cfg`:      interface of libcfg;
  * `omega_m`:  Omega_m;
  * `omega_l`:  Omega_l;
  * `omega_k`:  Omega_k;
  * `eos_w`:    w.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_cosmo(const cfg_t *cfg, double *omega_m, double *omega_l,
    double *omega_k, double *eos_w) {
  /* Check OMEGA_M. */
  if (!cfg_is_set(cfg, omega_m)) {
    P_ERR(FMT_KEY(OMEGA_M) " is not set.\n");
    return POWSPEC_ERR_CFG;
  }
  if (*omega_m <= 0 || *omega_m > 1) {
    P_ERR(FMT_KEY(OMEGA_M) " must be > 0 and <= 1.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Check OMEGA_LAMBDA */
  if (!cfg_is_set(cfg, omega_l)) {
    *omega_l = 1 - *omega_m;
    *omega_k = 0;
  }
  else if (*omega_l < 0) {
    P_ERR(FMT_KEY(OMEGA_LAMBDA) " must be >= 0.\n");
    return POWSPEC_ERR_CFG;
  }
  else *omega_k = 1 - *omega_m - *omega_l;

  /* Check DE_EOS_W. */
  if (!cfg_is_set(cfg, eos_w)) *eos_w = -1;
  else if (*eos_w > -1 / (double) 3) {
    P_ERR(FMT_KEY(DE_EOS_W) " must be <= -1/3.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Finally, make sure that H^2 (z) > 0. */
  double w3 = *eos_w * 3;
  double widx = w3 + 1;
  if (*omega_k * pow(*omega_l * (-widx), widx / w3) <=
      *omega_l * w3 * pow(*omega_m, widx / w3)) {
    P_ERR("negative H^2 given the cosmological parameters.\n");
    return POWSPEC_ERR_CFG;
  }
  return 0;
}

/******************************************************************************
Function `check_box`:
  Verify specifications of the box for the catalogs to be placed in.
Arguments:
  * `cfg`:      interface of libcfg;
  * `issim`:    indicate whether the files are from cubic simulation boxes;
  * `bsize`:    box size;
  * `bpad`:     box padding fraction.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_box(const cfg_t *cfg, const bool issim, double **bsize,
    double **bpad) {
  int num;
  if (!(num = cfg_get_size(cfg, bsize))) {
    /* Box size is mandatory for simulation boxes. */
    if (issim) {
      P_ERR(FMT_KEY(BOX_SIZE) " is not set.\n");
      return POWSPEC_ERR_CFG;
    }

    /* Determine box size automatically using BOX_PAD. */
    int npad;
    if (!(npad = cfg_get_size(cfg, bpad))) {
      *bpad = malloc(sizeof(double) * 3);
      if (!(*bpad)) {
        P_ERR("failed to allocate memory for box padding.\n");
        return POWSPEC_ERR_MEMORY;
      }
      (*bpad)[0] = (*bpad)[1] = (*bpad)[2] = DEFAULT_BOX_PAD;
    }
    else if (npad == 1) {
      if (**bpad < 0) {
        P_ERR(FMT_KEY(BOX_PAD) " must be >= 0.\n");
        return POWSPEC_ERR_CFG;
      }
      double *tmp = realloc(*bpad, sizeof(double) * 3);
      if (!tmp) {
        P_ERR("failed to allocate memory for box padding.\n");
        return POWSPEC_ERR_MEMORY;
      }
      *bpad = tmp;
      (*bpad)[1] = (*bpad)[2] = (*bpad)[0];
    }
    else if (npad == 3) {
      if ((*bpad)[0] < 0 || (*bpad)[1] < 0 || (*bpad)[2] < 0) {
        P_ERR(FMT_KEY(BOX_PAD) " must be >= 0.\n");
        return POWSPEC_ERR_CFG;
      }
    }
    else {
      P_ERR("wrong dimension of " FMT_KEY(BOX_PAD));
      return POWSPEC_ERR_CFG;
    }
  }
  else if (num == 1) {
    if (**bsize <= 0) {
      P_ERR(FMT_KEY(BOX_SIZE) " must be > 0.\n");
      return POWSPEC_ERR_CFG;
    }
    double *tmp = realloc(*bsize, sizeof(double) * 3);
    if (!tmp) {
      P_ERR("failed to allocate memory for the box size.\n");
      return POWSPEC_ERR_MEMORY;
    }
    *bsize = tmp;
    (*bsize)[1] = (*bsize)[2] = (*bsize)[0];
  }
  else if (num == 3) {
    if ((*bsize)[0] <= 0 || (*bsize)[1] <= 0 || (*bsize)[2] <= 0) {
      P_ERR(FMT_KEY(BOX_SIZE) " must be > 0.\n");
      return POWSPEC_ERR_CFG;
    }
  }
  else {
    P_ERR("wrong dimension of " FMT_KEY(BOX_SIZE));
    return POWSPEC_ERR_CFG;
  }
  return 0;
}

/******************************************************************************
Function `conf_verify`:
  Verify configuration parameters.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int conf_verify(const cfg_t *cfg, CONF *conf) {
  int i, e, num;
  /* Check CUBIC_SIM first, since it decides whether random is needed. */
  if (!cfg_is_set(cfg, &conf->issim)) {
    P_ERR(FMT_KEY(CUBIC_SIM) " is not set.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Check the outputs for the required number of data catalogs. */
  /* Check OVERWRITE. */
  if (!cfg_is_set(cfg, &conf->ovwrite)) conf->ovwrite = DEFAULT_OVERWRITE;

  /* Check OUTPUT_CROSS. */
  conf->iscross = false;
  conf->ndata = 0;
  if (cfg_is_set(cfg, &conf->ocross) && *conf->ocross) {
    if ((e = check_output(conf->ocross, "OUTPUT_CROSS", conf->ovwrite)))
      return e;
    conf->iscross = true;
    conf->ndata = 2;
  }

  /* Check OUTPUT_AUTO. */
  conf->isauto[0] = conf->isauto[1] = false;
  if ((num = cfg_get_size(cfg, &conf->oauto))) {
    if (!conf->iscross) conf->ndata = (num > 2) ? 2 : num;
    if (num > 2) {
      P_WRN("omitting the following " FMT_KEY(OUTPUT_AUTO) ":\n");
      for (i = 2; i < num; i++) fprintf(stderr, "%s\n", conf->oauto[i]);
      num = 2;
    }
    for (i = 0; i < num; i++) {
      /* OUTPUT_AUTO can be omitted by passing an empty string. */
      if (conf->oauto[i][0]) {
        if ((e = check_output(conf->oauto[i], "OUTPUT_AUTO", conf->ovwrite)))
          return e;
        conf->isauto[i] = true;
      }
    }
  }
  if (!conf->ndata) {
    P_ERR("no output file specified.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Check format settings of DATA_CATALOG. */
  if ((e = check_file_fmt(cfg, conf->ndata, "DATA", conf->issim,
      &conf->dfname, conf->has_asc, &conf->dftype, &conf->dskip,
      &conf->dcmt, &conf->dfmtr, &conf->dpos, &conf->dwcomp, &conf->dwfkp,
      &conf->dnz, &conf->dsel, &conf->dcnvt))) return e;

  for (i = 0; i < conf->ndata; i++)
    if ((e = check_input(conf->dfname[i], "DATA_CATALOG"))) return e;

  if (!conf->issim) {
    /* Check format settings of RAND_CATALOG. */
    if ((e = check_file_fmt(cfg, conf->ndata, "RAND", conf->issim,
        &conf->rfname, conf->has_asc + 1, &conf->rftype, &conf->rskip,
        &conf->rcmt, &conf->rfmtr, &conf->rpos, &conf->rwcomp, &conf->rwfkp,
        &conf->rnz, &conf->rsel, &conf->rcnvt))) return e;

    for (i = 0; i < conf->ndata; i++)
      if ((e = check_input(conf->dfname[i], "RAND_CATALOG"))) return e;

    /* Check DATA_NZ and RAND_NZ. */
    if (!(conf->dnz) && !(conf->rnz)) {
      P_ERR(FMT_KEY(DATA_NZ) " and " FMT_KEY(RAND_NZ) " are both not set.\n");
      return POWSPEC_ERR_CFG;
    }
  }

  /* Check if coordinate coversion is required. */
  conf->cnvt = false;
  if (!conf->issim) {
    if (conf->dcnvt) {
      for (i = 0; i < conf->ndata; i++) {
        if (conf->dcnvt) {
          conf->cnvt = true;
          break;
        }
      }
    }
    else conf->cnvt = DEFAULT_CONVERT;
    if (!conf->cnvt) {
      if (conf->rcnvt) {
        for (i = 0; i < conf->ndata; i++) {
          if (conf->rcnvt) {
            conf->cnvt = true;
            break;
          }
        }
      }
      else conf->cnvt = DEFAULT_CONVERT;
    }
  }
  /* Check the fiducial cosmology. */
  if (conf->cnvt) {
    /* Check Z_CMVDST_CNVT. */
    if (cfg_is_set(cfg, &conf->fcdst)) {
      if ((e = check_input(conf->fcdst, "Z_CMVDST_CNVT"))) return e;
    }
    else {
      if ((e = check_cosmo(cfg, &conf->omega_m, &conf->omega_l, &conf->omega_k,
          &conf->eos_w))) return e;

      /* Check CMVDST_ERR. */
      if (!cfg_is_set(cfg, &conf->ecdst)) conf->ecdst = DEFAULT_CMVDST_ERR;
      else if (conf->ecdst < DOUBLE_EPSILON) {
        P_ERR(FMT_KEY(CMVDST_ERR) " is smaller than the machine epsilon.\n");
        return POWSPEC_ERR_CFG;
      }
    }
  }

  /* Check LINE_OF_SIGHT. */
  if (conf->issim) {
    if ((num = cfg_get_size(cfg, &conf->los))) {
      if (num != 3) {
        P_ERR(FMT_KEY(LINE_OF_SIGHT) " must be a 3-element array.\n");
        return POWSPEC_ERR_CFG;
      }
      double tmp = conf->los[0] * conf->los[0] + conf->los[1] * conf->los[1]
        + conf->los[2] * conf->los[2];
      if (tmp < 1 - DOUBLE_TOL || tmp > 1 + DOUBLE_TOL) {
        P_ERR(FMT_KEY(LINE_OF_SIGHT) " must be a unit vector.\n");
        return POWSPEC_ERR_CFG;
      }
    }
    else {      /* default line of sight: [0,0,1] */
      conf->los = calloc(3, sizeof(double));
      if (!conf->los) {
        P_ERR("failed to alocate memory for " FMT_KEY(LINE_OF_SIGHT) ".\n");
        return POWSPEC_ERR_MEMORY;
      }
      conf->los[2] = 1;
    }
  }

  /* Check BOX_SIZE and BOX_PAD. */
  if ((e = check_box(cfg, conf->issim, &conf->bsize, &conf->bpad))) return e;

  /* Check GRID_SIZE. */
  if (!cfg_is_set(cfg, &conf->gsize)) {
    P_ERR(FMT_KEY(GRID_SIZE) " is not set.\n");
    return POWSPEC_ERR_CFG;
  }
  else if (conf->gsize <= 1) {
    P_ERR(FMT_KEY(GRID_SIZE) " must be > 1.\n");
    return POWSPEC_ERR_CFG;
  }
  else if (conf->gsize > POWSPEC_MAX_GSIZE) {
    P_ERR(FMT_KEY(GRID_SIZE) " must be smaller than the preset limit: %d.\n",
        POWSPEC_MAX_GSIZE);
    return POWSPEC_ERR_CFG;
  }

  /* Check PARTICLE_ASSIGN. */
  if (!cfg_is_set(cfg, &conf->assign)) conf->assign = DEFAULT_PARTICLE_ASSIGN;
  else {
    switch (conf->assign) {
      case POWSPEC_ASSIGN_NGP:
      case POWSPEC_ASSIGN_CIC:
      case POWSPEC_ASSIGN_TSC:
      case POWSPEC_ASSIGN_PCS:
        break;
      default:
        P_ERR("invalid " FMT_KEY(PARTICLE_ASSIGN) ": %d\n", conf->assign);
        return POWSPEC_ERR_CFG;
    }
  }

  /* Check GRID_INTERLACE. */
  if (!cfg_is_set(cfg, &conf->intlace)) conf->intlace = DEFAULT_GRID_INTERLACE;

  /* Check MULTIPOLE. */
  if (!(conf->npole = cfg_get_size(cfg, &conf->poles))) {
    P_ERR(FMT_KEY(MULTIPOLE) " is not set.\n");
    return POWSPEC_ERR_CFG;
  }
  /* Sort multipoles and remove duplicates. */
  if (conf->npole > 1) {
    /* 5-line insertion sort from https://doi.org/10.1145/3812.315108 */
    for (num = 1; num < conf->npole; num++) {
      int tmp = conf->poles[num];
      for (i = num; i > 0 && conf->poles[i - 1] > tmp; i--)
        conf->poles[i] = conf->poles[i - 1];
      conf->poles[i] = tmp;
    }
    /* Remove duplicates from the sorted array. */
    i = 0;
    for (num = 1; num < conf->npole; num++) {
      if (conf->poles[num] != conf->poles[i]) {
        i++;
        conf->poles[i] = conf->poles[num];
      }
    }
    conf->npole = i + 1;
  }
  if (conf->poles[0] < 0 || conf->poles[conf->npole - 1] > POWSPEC_MAX_ELL) {
    P_ERR(FMT_KEY(MULTIPOLE) " must be between 0 and %d.\n", POWSPEC_MAX_ELL);
    return POWSPEC_ERR_CFG;
  }

  /* Check KMIN. */
  if (!cfg_is_set(cfg, &conf->kmin)) conf->kmin = DEFAULT_KMIN;
  else if (conf->kmin < 0) {
    P_ERR(FMT_KEY(KMIN) " must be >= 0.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Check LOG_SCALE. */
  if (!cfg_is_set(cfg, &conf->logscale)) conf->logscale = DEFAULT_LOG_SCALE;
  if (conf->logscale && conf->kmin == 0) {
    P_ERR(FMT_KEY(KMIN) " must be > 0 for log scale.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Check BIN_SIZE. */
  if (!cfg_is_set(cfg, &conf->kbin)) {
    P_ERR(FMT_KEY(BIN_SIZE) " is not set.\n");
    return POWSPEC_ERR_CFG;
  }
  else if (conf->kbin <= 0) {
    P_ERR(FMT_KEY(BIN_SIZE) " must be > 0.\n");
    return POWSPEC_ERR_CFG;
  }

  /* Check KMAX. */
  if (!cfg_is_set(cfg, &conf->kmax)) conf->kmax = POWSPEC_KMAX_UNSET_VAL;
  else {
    if (conf->kmax <= 0) {
      P_ERR(FMT_KEY(KMAX) " must be > 0.\n");
      return POWSPEC_ERR_CFG;
    }
    if (conf->logscale) {
      if (log10(conf->kmax) < log10(conf->kmin) + conf->kbin) {
        P_ERR("log10(" FMT_KEY(KMAX) ") must not be smaller than log10("
            FMT_KEY(KMIN) ") + " FMT_KEY(BIN_SIZE) ".\n");
        return POWSPEC_ERR_CFG;
      }
    }
    else if (conf->kmax < conf->kmin + conf->kbin) {
      P_ERR(FMT_KEY(KMAX) " must not be smaller than " FMT_KEY(KMIN)
          " + " FMT_KEY(BIN_SIZE) ".\n");
      return POWSPEC_ERR_CFG;
    }
  }

  /* Check OUTPUT_HEADER. */
  if (!cfg_is_set(cfg, &conf->oheader)) conf->oheader = DEFAULT_HEADER;

  /* Check VERBOSE. */
  if (!cfg_is_set(cfg, &conf->verbose)) conf->verbose = DEFAULT_VERBOSE;
  return 0;
}


/*============================================================================*\
                      Function for printing configurations
\*============================================================================*/

/******************************************************************************
Function `conf_print`:
  Print configuration parameters.
Arguments:
  * `conf`:     structure for storing configurations.
******************************************************************************/
static void conf_print(const CONF *conf) {
  /* Configuration file */
  printf("\n  CONFIG_FILE     = %s", conf->fconf);

  /* Data catalog. */
  const bool twocat = (conf->ndata == 2) ? true : false;
  printf("\n  DATA_CATALOG    = %s", conf->dfname[0]);
  if (twocat) printf("\n                    %s", conf->dfname[1]);
  int tmp = (conf->dftype) ? conf->dftype[0] : DEFAULT_FILE_FORMAT;
  printf("\n  DATA_FORMAT     = %s", powspec_ffmt_names[tmp]);
  if (twocat) {
    tmp = (conf->dftype) ? conf->dftype[1] : DEFAULT_FILE_FORMAT;
    printf(" , %s", powspec_ffmt_names[tmp]);
  }

  if (conf->has_asc[0]) {
    if (conf->dskip) {
      printf("\n  DATA_SKIP       = %ld", conf->dskip[0]);
      if (twocat) printf(" , %ld", conf->dskip[1]);
    }
    if (conf->dcmt) {
      printf("\n  DATA_COMMENT    = '%c'", conf->dcmt[0]);
      if (twocat) printf(" , '%c'", conf->dcmt[1]);
    }
    printf("\n  DATA_FORMATTER  = %s", conf->dfmtr[0]);
    if (twocat) printf("\n                    %s", conf->dfmtr[1]);
  }
  printf("\n  DATA_POSITION   = %s , %s , %s",
      conf->dpos[0], conf->dpos[1], conf->dpos[2]);
  if (twocat) printf("\n                    %s , %s , %s",
      conf->dpos[3], conf->dpos[4], conf->dpos[5]);

  if (conf->dwcomp) {
    printf("\n  DATA_WT_COMP    = %s", conf->dwcomp[0]);
    if (twocat) printf(" , %s", conf->dwcomp[1]);
  }

  if (!conf->issim) {
    if (conf->dwfkp) {
      printf("\n  DATA_WT_FKP     = %s", conf->dwfkp[0]);
      if (twocat) printf(" , %s", conf->dwfkp[1]);
    }
    if (conf->dnz) {
      printf("\n  DATA_NZ         = %s", conf->dnz[0]);
      if (twocat) printf(" , %s", conf->dnz[1]);
    }
  }

  if (conf->dsel) {
    printf("\n  DATA_SELECTION  = %s", conf->dsel[0]);
    if (twocat) printf("\n                    %s", conf->dsel[1]);
  }
  if (!conf->issim) {
    if (conf->dcnvt) {
      printf("\n  DATA_CONVERT    = %c", conf->dcnvt[0] ? 'T' : 'F');
      if (twocat) printf(" , %c", conf->dcnvt[1] ? 'T' : 'F');
    }
    else {
      printf("\n  DATA_CONVERT    = F");
      if (twocat) printf(" , F");
    }
  }

  /* Random catalog. */
  if (!conf->issim) {
    printf("\n  RAND_CATALOG    = %s", conf->rfname[0]);
    if (twocat) printf("\n                    %s", conf->rfname[1]);
    tmp = (conf->rftype) ? conf->rftype[0] : DEFAULT_FILE_FORMAT;
    printf("\n  RAND_FORMAT     = %s", powspec_ffmt_names[tmp]);
    if (twocat) {
      tmp = (conf->rftype) ? conf->rftype[1] : DEFAULT_FILE_FORMAT;
      printf(" , %s", powspec_ffmt_names[tmp]);
    }
  
    if (conf->has_asc[1]) {
      if (conf->rskip) {
        printf("\n  RAND_SKIP       = %ld", conf->rskip[0]);
        if (twocat) printf(" , %ld", conf->rskip[1]);
      }
      if (conf->rcmt) {
        printf("\n  RAND_COMMENT    = '%c'", conf->rcmt[0]);
        if (twocat) printf(" , '%c'", conf->rcmt[1]);
      }
      printf("\n  RAND_FORMATTER  = %s", conf->rfmtr[0]);
      if (twocat) printf("\n                    %s", conf->rfmtr[1]);
    }
    printf("\n  RAND_POSITION   = %s , %s , %s",
        conf->rpos[0], conf->rpos[1], conf->rpos[2]);
    if (twocat) printf("\n                    %s , %s , %s",
        conf->rpos[3], conf->rpos[4], conf->rpos[5]);
  
    if (conf->rwcomp) {
      printf("\n  RAND_WT_COMP    = %s", conf->rwcomp[0]);
      if (twocat) printf(" , %s", conf->rwcomp[1]);
    }
    if (conf->rwfkp) {
      printf("\n  RAND_WT_FKP     = %s", conf->rwfkp[0]);
      if (twocat) printf(" , %s", conf->rwfkp[1]);
    }

    if (conf->rnz) {
      printf("\n  RAND_NZ         = %s", conf->rnz[0]);
      if (twocat) printf(" , %s", conf->rnz[1]);
    }
  
    if (conf->rsel) {
      printf("\n  RAND_SELECTION  = %s", conf->rsel[0]);
      if (twocat) printf("\n                    %s", conf->rsel[1]);
    }
    if (!conf->issim) {
      if (conf->rcnvt) {
        printf("\n  RAND_CONVERT    = %c", conf->rcnvt[0] ? 'T' : 'F');
        if (twocat) printf(" , %c", conf->rcnvt[1] ? 'T' : 'F');
      }
      else {
        printf("\n  RAND_CONVERT    = F");
        if (twocat) printf(" , F");
      }
    }
  }

  /* Fiducial cosmology. */
  if (conf->cnvt) {
    if (conf->fcdst) printf("\n  Z_CMVDST_CNVT   = %s", conf->fcdst);
    else {
      printf("\n  OMEGA_M         = " OFMT_DBL, conf->omega_m);
      printf("\n  OMEGA_LAMBDA    = " OFMT_DBL, conf->omega_l);
      if (conf->eos_w != -1)
        printf("\n  DE_EOS_W        = " OFMT_DBL, conf->eos_w);
      printf("\n  CMVDST_ERR      = " OFMT_DBL, conf->ecdst);
    }
  }

  /* Power spectra evaluation. */
  printf("\n  CUBIC_SIM       = %c", conf->issim ? 'T' : 'F');
  if (conf->issim) {
    printf("\n  LINE_OF_SIGHT   = " OFMT_DBL " , " OFMT_DBL " , " OFMT_DBL,
        conf->los[0], conf->los[1], conf->los[2]);
  }
  if (conf->bsize) {
    printf("\n  BOX_SIZE        = " OFMT_DBL " , " OFMT_DBL " , " OFMT_DBL,
        conf->bsize[0], conf->bsize[1], conf->bsize[2]);
  }
  else {
    printf("\n  BOX_PAD         = " OFMT_DBL " , " OFMT_DBL " , " OFMT_DBL,
        conf->bpad[0], conf->bpad[1], conf->bpad[2]);
  }
  printf("\n  GRID_SIZE       = %d", conf->gsize);
  printf("\n  PARTICLE_ASSIGN = %s", powspec_assign_names[conf->assign]);
  printf("\n  GRID_INTERLACE  = %c", conf->intlace ? 'T' : 'F');

  /* Output. */
  printf("\n  MULTIPOLE       = %d", conf->poles[0]);
  for (int i = 1; i < conf->npole; i++) printf(" , %d", conf->poles[i]);
  printf("\n  KMIN            = " OFMT_DBL, conf->kmin);
  if (POWSPEC_KMAX_ISSET(conf->kmax))
    printf("\n  KMAX            = " OFMT_DBL, conf->kmax);
  printf("\n  LOG_SCALE       = %c", conf->logscale ? 'T' : 'F');
  printf("\n  BIN_SIZE        = " OFMT_DBL, conf->kbin);
  if (conf->isauto[0] || conf->isauto[1]) {
    printf("\n  OUTPUT_AUTO     = %s", conf->oauto[0]);
    if (twocat) printf("\n                    %s", conf->oauto[1]);
  }
  if (conf->iscross) printf("\n  OUTPUT_CROSS    = %s", conf->ocross);
  printf("\n  OUTPUT_HEADER   = %c", conf->oheader ? 'T' : 'F');
  printf("\n  OVERWRITE       = %d\n", conf->ovwrite);
}


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
CONF *load_conf(const int argc, char *const *argv) {
  CONF *conf = conf_init();
  if (!conf) return NULL;

  cfg_t *cfg = conf_read(conf, argc, argv);
  if (!cfg) {
    conf_destroy(conf);
    return NULL;
  }

  printf("Loading configurations ...");
  fflush(stdout);

  if (conf_verify(cfg, conf)) {
    conf_destroy(conf);
    cfg_destroy(cfg);
    return NULL;
  }

  if (conf->verbose) conf_print(conf);

  if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
  cfg_destroy(cfg);

  /* Compute the wave number limits in logarithm scale after printing. */
  if (conf->logscale) {
    conf->kmin = log10(conf->kmin);
    if (POWSPEC_KMAX_ISSET(conf->kmax)) conf->kmax = log10(conf->kmax);
  }

  printf(FMT_DONE);
  return conf;
}

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf) {
  if (!conf) return;
  if (conf->dfname) {
    free(*conf->dfname); free(conf->dfname);
  }
  if (conf->rfname) {
    free(*conf->rfname); free(conf->rfname);
  }
  if (conf->dftype) free(conf->dftype);
  if (conf->rftype) free(conf->rftype);
  if (conf->dskip) free(conf->dskip);
  if (conf->rskip) free(conf->rskip);
  if (conf->dcmt) free(conf->dcmt);
  if (conf->rcmt) free(conf->rcmt);
  if (conf->dfmtr) {
    free(*conf->dfmtr); free(conf->dfmtr);
  }
  if (conf->rfmtr) {
    free(*conf->rfmtr); free(conf->rfmtr);
  }
  if (conf->dpos) {
    free(*conf->dpos); free(conf->dpos);
  }
  if (conf->rpos) {
    free(*conf->rpos); free(conf->rpos);
  }
  if (conf->dwcomp) {
    free(*conf->dwcomp); free(conf->dwcomp);
  }
  if (conf->rwcomp) {
    free(*conf->rwcomp); free(conf->rwcomp);
  }
  if (conf->dwfkp) {
    free(*conf->dwfkp); free(conf->dwfkp);
  }
  if (conf->rwfkp) {
    free(*conf->rwfkp); free(conf->rwfkp);
  }
  if (conf->dnz) {
    free(*conf->dnz); free(conf->dnz);
  }
  if (conf->rnz) {
    free(*conf->rnz); free(conf->rnz);
  }
  if (conf->dsel) {
    free(*conf->dsel); free(conf->dsel);
  }
  if (conf->rsel) {
    free(*conf->rsel); free(conf->rsel);
  }
  if (conf->dcnvt) free(conf->dcnvt);
  if (conf->rcnvt) free(conf->rcnvt);
  if (conf->fcdst) free(conf->fcdst);
  if (conf->los) free(conf->los);
  if (conf->bsize) free(conf->bsize);
  if (conf->bpad) free(conf->bpad);
  if (conf->poles) free(conf->poles);
  if (conf->oauto) {
    free(*conf->oauto); free(conf->oauto);
  }
  if (conf->ocross) free(conf->ocross);
  free(conf);
}

