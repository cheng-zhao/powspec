/*******************************************************************************
* read_ascii.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include "ascii_fmtr.h"
#include "libast.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                  Definitions for libast compatible data types
\*============================================================================*/
#define AST_DTYPE_NULL          0               /* not supported by libast */
#define AST_DTYPE_STR_SPACE     (-AST_DTYPE_STRING)  /* string with spaces */
#define AST_DTYPE(x)            ((x < 0) ? -(unsigned)(x) : x)   /* abs(x) */

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
#define P_AST_ERR(ast)  ast_perror(ast, stderr, FMT_ERR);
#define CLEAN_PTR  {                                                    \
  ast_destroy(ast_pos[0]); ast_destroy(ast_pos[1]);                     \
  ast_destroy(ast_pos[2]); ast_destroy(ast_sel);                        \
  ast_destroy(ast_wc); ast_destroy(ast_wfkp); ast_destroy(ast_nz);      \
  ascii_arg_destroy(arg, nc); free(col);                                \
  if (chunk) free(chunk); if (fp) fclose(fp); if (dat) free(dat);       \
}

#ifdef OMP
#define POWSPEC_QUIT(x) {                                               \
  printf(FMT_FAIL);                                                     \
  P_EXT("failed to read the ASCII file\n");                             \
  exit(x);                                                              \
}
#endif

/*============================================================================*\
                        Data structure for ASCII columns
\*============================================================================*/

/* Structure for recording libast compatible data types. */
typedef struct {
  int dtype;            /* AST_DTYPE_NULL or libast compatible data type */
  union {               /* variable for storing the value */
    int ival; long lval; float fval; double dval;
    struct ast_string_struct_t { int len; const char *str; } sval;
  } v;
} asc_col_t;


/*============================================================================*\
                      Functions for parsing ASCII columns
\*============================================================================*/

/******************************************************************************
Function `ascii_col_init`:
  Initialise columns according to the arguments parsed from the formatter.
Arguments:
  * `arg`:      arguments parsed from the formatter;
  * `num`:      number of parsed arguments;
  * `rnum`:     number of arguments that are not suppressed.
Return:
  Pointer to the structure array for ASCII columns on success; NULL on error.
******************************************************************************/
static asc_col_t *ascii_col_init(asc_arg_t *arg, const int num,
    const int rnum) {
  if (rnum <= 0) return NULL;
  asc_col_t *col = malloc(rnum * sizeof(asc_col_t));
  if (!col) return NULL;

  int j = -1;
  for (int i = 0; i < num; i++) {
    if (arg[i].dtype == ASCII_DTYPE_SKIP) continue;
    if (++j >= rnum) {
      P_ERR("unknown error for identifying ASCII columns\n");
      free(col);
      return NULL;
    }
    /* Convert the `fscanf` type to libast type. */
    switch (arg[i].dtype) {
      case ASCII_DTYPE_INT: col[j].dtype = AST_DTYPE_INT; break;
      case ASCII_DTYPE_LONG: col[j].dtype = AST_DTYPE_LONG; break;
      case ASCII_DTYPE_FLT: col[j].dtype = AST_DTYPE_FLOAT; break;
      case ASCII_DTYPE_DBL: col[j].dtype = AST_DTYPE_DOUBLE; break;
      case ASCII_DTYPE_STR: col[j].dtype = AST_DTYPE_STRING; break;
      case ASCII_DTYPE_CHAR: col[j].dtype = AST_DTYPE_STRING; break;
      default:          /* data types that are not supported by libast */
        col[j].dtype = AST_DTYPE_NULL;
        P_WRN("unsupported formatter `%s', "
            "the corresponding column is not used\n", arg[i].fmtr);
        break;
    }

    /* Add a suppressing symbol '*' to the formatter string. */
    if (col[j].dtype == AST_DTYPE_NULL || col[j].dtype == AST_DTYPE_STRING) {
      int len = strlen(arg[i].fmtr) + 1;        /* fmtr is null terminated */
      char *tmp = realloc(arg[i].fmtr, (len + 1) * sizeof(char));
      if (!tmp) {
        free(col);
        return NULL;
      }
      arg[i].fmtr = tmp;

      /* Parse the formatter string to find the right place for '*'. */
      char *fmt = tmp;
      while (*fmt) {
        char c = *fmt++;
        if (c != '%') continue;
        c = *fmt++;
        if (c == '%') continue;
        /* Check if the formatter parses whitespaces. */
        if (col[j].dtype == AST_DTYPE_STRING && fmt - tmp > 2 && !isspace(*tmp))
          col[j].dtype = AST_DTYPE_STR_SPACE;

        /* Now add '*'. */
        memmove(fmt, fmt - 1, tmp + len - fmt + 1);
        *(fmt - 1) = '*';

        /* The formatter can still parse whitespaces with '[]'. */
        if (col[j].dtype == AST_DTYPE_STRING) {
          while (c != '[' && c != '\0') c = *fmt++;
          if (c == '[') {
            if (*fmt == '^') ++fmt;
            if (*fmt == ']') ++fmt;
            while ((c = *fmt++) != '\0' && c != ']' && c != ' ');
            if (c == ' ') col[j].dtype = AST_DTYPE_STR_SPACE;
          }
        }
        break;
      }
    }
  }
  return col;
}

/******************************************************************************
Function `ascii_read_line`:
  Read a string line to column variables.
Arguments:
  * `line`:     the line to be read;
  * `arg`:      arguments parsed from the formatter;
  * `num`:      number of columns to be read;
  * `col`:      variables and their data types for each column;
  * `end`:      pointer to the first character that is not interpreted.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int ascii_read_line(const char *line, const asc_arg_t *arg,
    const int num, asc_col_t *col, const char **end) {
  int j = 0;
  for (int i = 0; i < num; i++) {
    int n = 0;
    /* User-suppressed columns. */
    if (arg[i].dtype == ASCII_DTYPE_SKIP) {
      if (sscanf(line, arg[i].fmtr, &n) != 0 || n == 0) {
        *end = line;
        return POWSPEC_ERR_ASCII;
      }
    }
    else {
      switch (col[j].dtype) {
        case AST_DTYPE_NULL:
          /* Suppressed columns due to compatibility with libast. */
          if (sscanf(line, arg[i].fmtr, &n) != 0 || n == 0) {
            *end = line;
            return POWSPEC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_INT:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.ival), &n) != 1 || n == 0) {
            *end = line;
            return POWSPEC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_LONG:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.lval), &n) != 1 || n == 0) {
            *end = line;
            return POWSPEC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_FLOAT:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.fval), &n) != 1 || n == 0) {
            *end = line;
            return POWSPEC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_DOUBLE:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.dval), &n) != 1 || n == 0) {
            *end = line;
            return POWSPEC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_STRING:
          while (isspace(*line)) line++;        /* omit whitespaces */
        case AST_DTYPE_STR_SPACE:
          /* Do not actually save the string, but compute the position. */
          if (sscanf(line, arg[i].fmtr, &n) != 0 || n == 0) {
            *end = line;
            return POWSPEC_ERR_ASCII;
          }
          col[j].v.sval.str = line;
          col[j].v.sval.len = n;
          break;
        default:
          P_ERR("unknown data type for column: ${%d}\n", j + 1);
          *end = line;
          return POWSPEC_ERR_UNKNOWN;
      }
      j++;
    }
    line += n;
  }
  return 0;
}

/******************************************************************************
Function `ascii_read_sel`:
  Read the selection criteria based on the AST and ASCII columns.
Arguments:
  * `ast`:      abstract syntax tree for the selection criteria;
  * `col`:      columns to be parsed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int ascii_read_sel(ast_t *ast, const asc_col_t *col, bool *res) {
  for (int i = 0; i < ast->nvar; i++) {
    long idx = ast->vidx[i];
    const asc_col_t *c = col + (idx - 1);
    switch (AST_DTYPE(c->dtype)) {
      case AST_DTYPE_INT:
        if (ast_set_var(ast, idx, &c->v.ival, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_LONG:
        if (ast_set_var(ast, idx, &c->v.lval, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_FLOAT:
        if (ast_set_var(ast, idx, &c->v.fval, 0, AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_DOUBLE:
        if (ast_set_var(ast, idx, &c->v.dval, 0, AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_STRING:
        if (ast_set_var(ast, idx, c->v.sval.str, c->v.sval.len,
              AST_DTYPE_STRING)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      default:
        P_ERR("column ${%ld} not appropriate for selection\n", idx);
        return POWSPEC_ERR_AST;
    }
  }

  bool val = false;
  if (ast_eval(ast, &val)) {
    P_AST_ERR(ast);
    return POWSPEC_ERR_AST;
  }
  *res = val;
  return 0;
}

/******************************************************************************
Function `ascii_read_double`:
  Read a double number based on the AST and ASCII columns.
Arguments:
  * `ast`:      abstract syntax tree for the number evaluation;
  * `col`:      columns to be parsed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int ascii_read_double(ast_t *ast, const asc_col_t *col,
    double *res) {
  for (int i = 0; i < ast->nvar; i++) {
    long idx = ast->vidx[i];
    const asc_col_t *c = col + (idx - 1);
    switch (AST_DTYPE(c->dtype)) {
      case AST_DTYPE_INT:
        if (ast_set_var(ast, idx, &c->v.ival, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_LONG:
        if (ast_set_var(ast, idx, &c->v.lval, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_FLOAT:
        if (ast_set_var(ast, idx, &c->v.fval, 0, AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      case AST_DTYPE_DOUBLE:
        if (ast_set_var(ast, idx, &c->v.dval, 0, AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast);
          return POWSPEC_ERR_AST;
        }
        break;
      default:
        P_ERR("column ${%ld} not appropriate for the numerical evaluation\n",
            idx);
        return POWSPEC_ERR_AST;
    }
  }

  double val = 0;
  if (ast_eval(ast, &val)) {
    P_AST_ERR(ast);
    return POWSPEC_ERR_AST;
  }
  *res = val;
  return 0;
}


/*============================================================================*\
                      Functions for reading file by chunks
\*============================================================================*/

/******************************************************************************
Function `chunk_resize`:
  Enlarge the size of a chunk.
Arguments:
  * `chunk`:    address of the chunk;
  * `size`:     size of the chunk.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int chunk_resize(char **chunk, size_t *size) {
  /* Assume the arguments are not NULL. */
  size_t num;
  if (!(*chunk)) num = POWSPEC_FILE_CHUNK;
  else {
    if (POWSPEC_MAX_CHUNK / 2 < *size) return POWSPEC_ERR_FILE;
    num = *size << 1;
  }

  char *tmp = realloc(*chunk, num * sizeof(char));
  if (!tmp) return POWSPEC_ERR_MEMORY;

  *chunk = tmp;
  *size = num;
  return 0;
}

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
    const int verb) {
  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    return POWSPEC_ERR_FILE;
  }

  if (verb) printf("\n  Reading samples from file: %s\n", fname);

  /* Prepare for the chunk. */
  char *chunk = NULL;
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the file by chunk\n");
    fclose(fp);
    return POWSPEC_ERR_MEMORY;
  }

  /* Allocate memory for the data. */
  size_t max = POWSPEC_DATA_INIT_NUM;
  double *nx = malloc(max * sizeof(double));
  double *ny = malloc(max * sizeof(double));
  if (!nx || !ny) {
    P_ERR("failed to allocate memory for the samples\n");
    fclose(fp); free(chunk);
    if (nx) free(nx);
    if (ny) free(ny);
    return POWSPEC_ERR_MEMORY;
  }

  size_t n, nread, nrest;
  n = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == POWSPEC_READ_COMMENT || *p == '\0') {   /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      if (sscanf(p, "%lf %lf", nx + n, ny + n) != 2) {
        P_ERR("failed to read line: %s\n", p);
        fclose(fp); free(chunk); free(nx); free(ny);
        return POWSPEC_ERR_FILE;
      }

      /* Enlarge the memory for the data if necessary. */
      if (++n >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many samples in the file: `%s'\n", fname);
          fclose(fp); free(chunk); free(nx); free(ny);
          return POWSPEC_ERR_FILE;
        }
        max <<= 1;
        double *tmp = realloc(nx, sizeof(double) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the samples\n");
          fclose(fp); free(chunk); free(nx); free(ny);
          return POWSPEC_ERR_MEMORY;
        }
        nx = tmp;
        tmp = realloc(ny, sizeof(double) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the samples\n");
          fclose(fp); free(chunk); free(nx); free(ny);
          return POWSPEC_ERR_MEMORY;
        }
        ny = tmp;
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        fclose(fp); free(chunk); free(nx); free(ny);
        return POWSPEC_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    fclose(fp); free(chunk); free(nx); free(ny);
    return POWSPEC_ERR_FILE;
  }

  free(chunk);
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);
  if (verb) printf("  Number of samples: %zu\n", n);

  *x = nx;
  *y = ny;
  *num = n;
  return 0;
}

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
  * `sim`:      indicate whether the catalog is from a simulation box;
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
    const char *nz, const char *sel, const bool sim, DATA **data, size_t *num,
    double *sumw, double *I12, double *I22, const int verb) {
  /* Parse the formatter. */
  asc_arg_t *arg;
  int nc, rnc;
  nc = rnc = 0;
  if (!(arg = parse_ascii_fmtr(fmtr, &nc, &rnc))) return POWSPEC_ERR_ASCII;
  if (!rnc) {
    P_ERR("no column to be read given the formatter: `%s'\n", fmtr);
    return POWSPEC_ERR_ASCII;
  }
  if (rnc < 3) P_WRN("reading coordinates from less than 3 columns\n");

  /* Record libast compatible columns. */
  asc_col_t *col;
  if (!(col = ascii_col_init(arg, nc, rnc))) {
    ascii_arg_destroy(arg, nc);
    return POWSPEC_ERR_ASCII;
  }

  /* Initialise all variables for easy error handling. */
  FILE *fp = NULL;
  char *chunk = NULL;
  DATA *dat = NULL;
  ast_t *ast_pos[3] = {NULL, NULL, NULL};
  ast_t *ast_sel, *ast_wc, *ast_wfkp, *ast_nz;
  ast_sel = ast_wc = ast_wfkp = ast_nz = NULL;

  /* Construct the ASTs for positions and selections. */
  for (int i = 0; i < 3; i++) {
    if (!(ast_pos[i] = ast_init())) {
      P_AST_ERR(ast_pos[i]); CLEAN_PTR; return POWSPEC_ERR_AST;
    }
    if (ast_build(ast_pos[i], pos[i], AST_DTYPE_DOUBLE, true)) {
      P_AST_ERR(ast_pos[i]); CLEAN_PTR; return POWSPEC_ERR_AST;
    }
  }
  if (sel && *sel) {
    if (!(ast_sel = ast_init())) {
      P_AST_ERR(ast_sel); CLEAN_PTR; return POWSPEC_ERR_AST;
    }
    if (ast_build(ast_sel, sel, AST_DTYPE_BOOL, true)) {
      P_AST_ERR(ast_sel); CLEAN_PTR; return POWSPEC_ERR_AST;
    }
  }
  /* Construct the ASTs for weights. */
  if (wc && *wc) {
    if (!(ast_wc = ast_init())) {
      P_AST_ERR(ast_wc); CLEAN_PTR; return POWSPEC_ERR_AST;
    }
    if (ast_build(ast_wc, wc, AST_DTYPE_DOUBLE, true)) {
      P_AST_ERR(ast_wc); CLEAN_PTR; return POWSPEC_ERR_AST;
    }
  }
  if (!sim) {
    if (wfkp && *wfkp) {
      if (!(ast_wfkp = ast_init())) {
        P_AST_ERR(ast_wfkp); CLEAN_PTR; return POWSPEC_ERR_AST;
      }
      if (ast_build(ast_wfkp, wfkp, AST_DTYPE_DOUBLE, true)) {
        P_AST_ERR(ast_wfkp); CLEAN_PTR; return POWSPEC_ERR_AST;
      }
    }
    if (nz && *nz) {
      if (!(ast_nz = ast_init())) {
        P_AST_ERR(ast_nz); CLEAN_PTR; return POWSPEC_ERR_AST;
      }
      if (ast_build(ast_nz, nz, AST_DTYPE_DOUBLE, true)) {
        P_AST_ERR(ast_nz); CLEAN_PTR; return POWSPEC_ERR_AST;
      }
    }
  }

  /* Check number of variables for expressions. */
  int max_col = 0;
  for (int i = 0; i < 3; i++) {
    if (ast_pos[i]->nvar == 0) {
      P_ERR("the expression for position coordinate %d is a constant: `%s'\n",
          i + 1, pos[i]);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (rnc < ast_pos[i]->vidx[ast_pos[i]->nvar - 1]) {
      P_ERR("not enough columns for position coordinate %d: `%s'\n",
          i + 1, pos[i]);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (max_col < ast_pos[i]->vidx[ast_pos[i]->nvar - 1])
      max_col = ast_pos[i]->vidx[ast_pos[i]->nvar - 1];
  }
  if (ast_sel) {
    if (ast_sel->nvar == 0) {
      P_ERR("the expression for data selection is a constant: `%s'\n", sel);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (rnc < ast_sel->vidx[ast_sel->nvar - 1]) {
      P_ERR("not enough columns for data selection: `%s'\n", sel);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (max_col < ast_sel->vidx[ast_sel->nvar - 1])
      max_col = ast_sel->vidx[ast_sel->nvar - 1];
  }
  if (ast_wc) {
    if (ast_wc->nvar && rnc < ast_wc->vidx[ast_wc->nvar - 1]) {
      P_ERR("not enough columns for completeness weight: `%s'\n", wc);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (ast_wc->nvar && max_col < ast_wc->vidx[ast_wc->nvar - 1])
      max_col = ast_wc->vidx[ast_wc->nvar - 1];
  }
  if (ast_wfkp) {
    if (ast_wfkp->nvar && rnc < ast_wfkp->vidx[ast_wfkp->nvar - 1]) {
      P_ERR("not enough columns for FKP weight: `%s'\n", wfkp);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (ast_wfkp->nvar && max_col < ast_wfkp->vidx[ast_wfkp->nvar - 1])
      max_col = ast_wfkp->vidx[ast_wfkp->nvar - 1];
  }
  if (ast_nz) {
    if (ast_nz->nvar && rnc < ast_nz->vidx[ast_nz->nvar - 1]) {
      P_ERR("not enough columns for radial number density: `%s'\n", nz);
      CLEAN_PTR; return POWSPEC_ERR_CFG;
    }
    if (ast_nz->nvar && max_col < ast_nz->vidx[ast_nz->nvar - 1])
      max_col = ast_nz->vidx[ast_nz->nvar - 1];
  }
  /* Remove columns that are not needed. */
  if (max_col < rnc) {
    int i, j;
    for (i = j = 0; i < nc; i++) {
      if (arg[i].dtype != ASCII_DTYPE_SKIP) {
        if (++j == max_col) break;
      }
    }
    for (j = i + 1; j < nc; j++) if(arg[j].fmtr) free(arg[j].fmtr);
    nc = i + 1;
    rnc = max_col;
  }

  /* Prepare for the chunk. */
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the file by chunk\n");
    CLEAN_PTR; return POWSPEC_ERR_MEMORY;
  }

  /* Allocate memory for the data. */
  size_t max = POWSPEC_DATA_INIT_NUM;
  if (!(dat = malloc(max * sizeof(DATA)))) {
    P_ERR("failed to allocate memory for the data\n");
    CLEAN_PTR; return POWSPEC_ERR_MEMORY;
  }

  /* Dynamic allocations for OpenMP threads. */
#ifdef OMP
  const int nomp = omp_get_max_threads();
  /* Construct the ASCII columns for non-master threads. */
  asc_col_t **pcol = NULL;
  if (nomp > 1) {
    if (!(pcol = malloc(sizeof(asc_col_t *) * (nomp - 1)))) {
      P_ERR("failed to allocate memory for thread-private columns\n");
      POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
    }
    for (int j = 0; j < nomp - 1; j++) {
      if (!(pcol[j] = malloc(sizeof(asc_col_t) * rnc))) {
        P_ERR("failed to allocate memory for thread-private columns\n");
        POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
      }
      memcpy(pcol[j], col, sizeof(asc_col_t) * rnc);
    }
  }
  /* Construct the ASTs again for non-master threads. */
  ast_t **ast_ppos, **ast_psel, **ast_pwc, **ast_pwfkp, **ast_pnz;
  ast_ppos = ast_psel = ast_pwc = ast_pwfkp = ast_pnz = NULL;
  if (nomp > 1) {
    ast_ppos = malloc(sizeof(ast_t *) * (nomp - 1) * 3);
    ast_psel = malloc(sizeof(ast_t *) * (nomp - 1));
    ast_pwc = malloc(sizeof(ast_t *) * (nomp - 1));
    ast_pwfkp = malloc(sizeof(ast_t *) * (nomp - 1));
    ast_pnz = malloc(sizeof(ast_t *) * (nomp - 1));
    if (!ast_ppos || !ast_psel || !ast_pwc || !ast_pwfkp || !ast_pnz) {
      P_ERR("failed to allocate memory for thread-private ASTs\n");
      POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
    }
    for (int j = 0; j < nomp - 1; j++) {
      /* Positions. */
      for (int i = 0; i < 3; i++) {
        int k = j * 3 + i;
        if (!(ast_ppos[k] = ast_init())) {
          P_AST_ERR(ast_ppos[k]); POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
        if (ast_build(ast_ppos[k], pos[i], AST_DTYPE_DOUBLE, true)) {
          P_AST_ERR(ast_ppos[k]); POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
      }
      /* Selections. */
      if (sel && *sel) {
        if (!(ast_psel[j] = ast_init())) {
          P_AST_ERR(ast_psel[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
        if (ast_build(ast_psel[j], sel, AST_DTYPE_BOOL, true)) {
          P_AST_ERR(ast_psel[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
      }
      else ast_psel[j] = NULL;
      /* Weights. */
      if (wc && *wc) {
        if (!(ast_pwc[j] = ast_init())) {
          P_AST_ERR(ast_pwc[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
        if (ast_build(ast_pwc[j], wc, AST_DTYPE_DOUBLE, true)) {
          P_AST_ERR(ast_pwc[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
      }
      else ast_pwc[j] = NULL;
      if (!sim) {
        if (wfkp && *wfkp) {
          if (!(ast_pwfkp[j] = ast_init())) {
            P_AST_ERR(ast_pwfkp[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
          if (ast_build(ast_pwfkp[j], wfkp, AST_DTYPE_DOUBLE, true)) {
            P_AST_ERR(ast_pwfkp[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
        }
        else ast_pwfkp[j] = NULL;
        if (nz && *nz) {
          if (!(ast_pnz[j] = ast_init())) {
            P_AST_ERR(ast_pnz[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
          if (ast_build(ast_pnz[j], nz, AST_DTYPE_DOUBLE, true)) {
            P_AST_ERR(ast_pnz[j]); POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
        }
        else ast_pnz[j] = NULL;
      }
      else ast_pwfkp[j] = ast_pnz[j] = NULL;
    }
  }
  /* Construct the private data pool. */
  DATA **pdata = malloc(sizeof(DATA *) * nomp);
  size_t *pndata = calloc(nomp, sizeof(size_t));
  double *psumw = calloc(nomp, sizeof(double));
  double *psumw2 = calloc(nomp, sizeof(double));
  double *psumw2n = calloc(nomp, sizeof(double));
  if (!pdata || !pndata || !psumw || !psumw2 || !psumw2n) {
    P_ERR("failed to allocate memory for the thread-private data\n");
    POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
  }
  if (!(pdata[0] = malloc(sizeof(DATA) * nomp * POWSPEC_DATA_THREAD_NUM))) {
    P_ERR("failed to allocate memory for the thread-private data\n");
    POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
  }
  for (int j = 0; j < nomp; j++)
    pdata[j] = pdata[0] + j * POWSPEC_DATA_THREAD_NUM;
  /* Construct the pool for file lines. */
  size_t nlmax = POWSPEC_DATA_INIT_NUM;
  size_t nl = 0;
  char **lines = malloc(sizeof(char *) * nlmax);
  if (!lines) {
    P_ERR("failed to allocate memory for the thread-private lines\n");
    POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
  }
#endif

  /* Open the file for reading. */
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    CLEAN_PTR; return POWSPEC_ERR_FILE;
  }
  if (verb) printf("  Filename: %s\n", fname);

  size_t n, nline, nread, nrest;
  n = nline = nrest = 0;
  double sumwt, sumw2, sumw2n;
  sumwt = sumw2 = sumw2n = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      /* Skip header lines. */
      if (nline++ < skip) {
        p = endl + 1; continue;
      }
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == comment || *p == '\0') {        /* comment or empty */
        p = endl + 1; continue;
      }

#ifdef OMP
      /* Append the line to the pool. */
      lines[nl++] = p;
      /* Enlarge the pool if necessary. */
      if (nl >= nlmax) {
        if (SIZE_MAX / 2 < nlmax) {
          P_ERR("too many lines in the file\n");
          POWSPEC_QUIT(POWSPEC_ERR_FILE);
        }
        nlmax <<= 1;
        char **tmp = realloc(lines, sizeof(char *) * nlmax);
        if (!tmp) {
          P_ERR("failed to allocate memory for the thread-private lines\n");
          POWSPEC_QUIT(POWSPEC_ERR_MEMORY);
        }
        lines = tmp;
      }
#else
      /* Parse the line. */
      const char *stop = NULL;
      if (ascii_read_line(p, arg, nc, col, &stop)) {
        P_ERR("failed to read the line with format `%s':\n", fmtr);
        fprintf(stderr, "%s\n", p);
        if (stop > p) for (int k = 0; k < stop - p; k++) fprintf(stderr, " ");
        fprintf(stderr, "^\n");
        CLEAN_PTR; return POWSPEC_ERR_ASCII;
      }

      /* Apply selection. */
      if (ast_sel) {
        bool keep = false;
        if (ascii_read_sel(ast_sel, col, &keep)) {
          CLEAN_PTR; return POWSPEC_ERR_AST;
        }
        if (!keep) {
          p = endl + 1; continue;
        }
      }

      /* Record the coordinates. */
      for (int i = 0; i < 3; i++) {
        if (ascii_read_double(ast_pos[i], col, &(dat[n].x[i]))) {
          CLEAN_PTR; return POWSPEC_ERR_AST;
        }
      }

      /* Compute weights. */
      double twc;
      if (ast_wc) {
        if (ascii_read_double(ast_wc, col, &twc)) {
          CLEAN_PTR; return POWSPEC_ERR_AST;
        }
      }
      else twc = 1;
      if (!sim) {
        double twfkp, tnz;
        if (ast_wfkp) {
          if (ascii_read_double(ast_wfkp, col, &twfkp)) {
            CLEAN_PTR; return POWSPEC_ERR_AST;
          }
        }
        else twfkp = 1;
        if (ast_nz) {
          if (ascii_read_double(ast_nz, col, &tnz)) {
            CLEAN_PTR; return POWSPEC_ERR_AST;
          }
        }
        else tnz = 0;

        dat[n].w = twc * twfkp;
        sumwt += twc;
        sumw2 += dat[n].w * dat[n].w;
        sumw2n += twc * twfkp * twfkp * tnz;
      }
      else sumwt += (dat[n].w = twc);

      /* Enlarge the memory for the data if necessary. */
      if (++n >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many objects in the file: `%s'\n", fname);
          CLEAN_PTR;
          return POWSPEC_ERR_FILE;
        }
        max <<= 1;
        DATA *tmp = realloc(dat, sizeof(DATA) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the data\n");
          CLEAN_PTR;
          return POWSPEC_ERR_MEMORY;
        }
        dat = tmp;
      }
#endif
      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        CLEAN_PTR;
        return POWSPEC_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

#ifdef OMP
    /* Enlarge the memory for the data if necessary. */
    if (nl + n > max) {
      while (nl + n > max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many objects in the file: `%s'.\n", fname);
          CLEAN_PTR;
          return POWSPEC_ERR_FILE;
        }
        max <<= 1;
      }
      DATA *tmp = realloc(dat, sizeof(DATA) * max);
      if (!tmp) {
        P_ERR("failed to allocate memory for the data.\n");
        CLEAN_PTR;
        return POWSPEC_ERR_MEMORY;
      }
      dat = tmp;
    }

#pragma omp parallel num_threads(nomp) \
    firstprivate(ast_pos, ast_sel, ast_wc, ast_wfkp, ast_nz, col)
    {
      /* Redirect pointers to the private pools. */
      const int tid = omp_get_thread_num();
      if (tid > 0) {
        ast_pos[0] = ast_ppos[(tid - 1) * 3];
        ast_pos[1] = ast_ppos[(tid - 1) * 3 + 1];
        ast_pos[2] = ast_ppos[(tid - 1) * 3 + 2];
        ast_sel = ast_psel[tid - 1];
        ast_wc = ast_pwc[tid - 1];
        ast_wfkp = ast_pwfkp[tid - 1];
        ast_nz = ast_pnz[tid - 1];
        col = pcol[tid - 1];
      }
      DATA *pdat = pdata[tid];
      size_t *pnum = pndata + tid;
      /* Process lines in parallel. */
#pragma omp for
      for (size_t ii = 0; ii < nl; ii++) {
        /* Parse a line in the pool. */
        const char *pp = lines[ii];
        const char *stop = NULL;
        if (ascii_read_line(pp, arg, nc, col, &stop)) {
#pragma omp critical
          {
            P_ERR("failed to read the line with format `%s':\n", fmtr);
            fprintf(stderr, "%s\n", pp);
            if (stop > pp) for (int k = 0; k < stop - pp; k++)
              fprintf(stderr, " ");
            fprintf(stderr, "^\n");
          }
          POWSPEC_QUIT(POWSPEC_ERR_ASCII);
        }

        /* Apply selection. */
        if (ast_sel) {
          bool keep = false;
          if (ascii_read_sel(ast_sel, col, &keep)) {
            POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
          if (!keep) continue;
        }

        /* Record coordinates to the private data pool. */
        for (int i = 0; i < 3; i++) {
          if (ascii_read_double(ast_pos[i], col, &(pdat[*pnum].x[i]))) {
            POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
        }

        /* Compute weights. */
        double twc;
        if (ast_wc) {
          if (ascii_read_double(ast_wc, col, &twc)) {
            POWSPEC_QUIT(POWSPEC_ERR_AST);
          }
        }
        else twc = 1;
        if (!sim) {
          double twfkp, tnz;
          if (ast_wfkp) {
            if (ascii_read_double(ast_wfkp, col, &twfkp)) {
              POWSPEC_QUIT(POWSPEC_ERR_AST);
            }
          }
          else twfkp = 1;
          if (ast_nz) {
            if (ascii_read_double(ast_nz, col, &tnz)) {
              POWSPEC_QUIT(POWSPEC_ERR_AST);
            }
          }
          else tnz = 0;

          pdat[*pnum].w = twc * twfkp;
          psumw[tid] += twc;
          psumw2[tid] += pdat[*pnum].w * pdat[*pnum].w;
          psumw2n[tid] += twc * twfkp * twfkp * tnz;
        }
        else psumw[tid] += (pdat[*pnum].w = twc);

        /* Record the private data and clear the pool if necessary. */
        if (++(*pnum) >= POWSPEC_DATA_THREAD_NUM) {
#pragma omp critical
          {
            /* Enlarge the memory for the data if necessary. */
            if (n + POWSPEC_DATA_THREAD_NUM >= max) {
              if (SIZE_MAX / 2 < max) {
                P_ERR("too many objects in the file: `%s'\n", fname);
                POWSPEC_QUIT(POWSPEC_ERR_AST);
              }
              max <<= 1;
              if (max < n + POWSPEC_DATA_THREAD_NUM)
                max = n + POWSPEC_DATA_THREAD_NUM;
              DATA *tmp = realloc(dat, sizeof(DATA) * max);
              if (!tmp) {
                P_ERR("failed to allocate memory for the data\n");
                POWSPEC_QUIT(POWSPEC_ERR_AST);
              }
              dat = tmp;
            }
            memcpy(dat + n, pdat, sizeof(DATA) * POWSPEC_DATA_THREAD_NUM);
            n += POWSPEC_DATA_THREAD_NUM;
          }
          pndata[tid] = 0;
        }
      }
      nl = 0;
    }
#endif

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

#ifdef OMP
  /* Record the rest of the private data. */
  for (int i = 0; i < nomp; i++) {
    if (pndata[i]) {
      /* Enlarge the memory for the data if necessary. */
      if (n + pndata[i] >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many objects in the file: `%s'\n", fname);
          POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
        max <<= 1;
        if (max < n + pndata[i]) max = n + pndata[i];
        DATA *tmp = realloc(dat, sizeof(DATA) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the data\n");
          POWSPEC_QUIT(POWSPEC_ERR_AST);
        }
        dat = tmp;
      }
      memcpy(dat + n, pdata[i], sizeof(DATA) * pndata[i]);
      n += pndata[i];
    }
    /* Sum the weights. */
    sumwt += psumw[i];
    sumw2 += psumw2[i];
    sumw2n += psumw2n[i];
  }
  /* Release memory for thread-private data structures. */
  for (int i = 0; i < nomp - 1; i++) {
    free(pcol[i]);
    if (ast_psel[i]) ast_destroy(ast_psel[i]);
    if (ast_pwc[i]) ast_destroy(ast_pwc[i]);
    if (ast_pwfkp[i]) ast_destroy(ast_pwfkp[i]);
    if (ast_pnz[i]) ast_destroy(ast_pnz[i]);
  }
  free(pdata[0]); free(pdata); free(pndata);
  for (int i = 0; i < (nomp - 1) * 3; i++) ast_destroy(ast_ppos[i]);
  free(pcol); free(ast_ppos); free(ast_psel);
  free(ast_pwc); free(ast_pwfkp); free(ast_pnz);
  free(psumw); free(psumw2); free(psumw2n);
  free(lines);
#endif

  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    CLEAN_PTR; return POWSPEC_ERR_FILE;
  }

  if (verb) {
    printf("  Number of lines processed in total: %zu\n"
        "  Number of recorded objects: %zu\n", nline, n);
  }

  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);
  ast_destroy(ast_pos[0]); ast_destroy(ast_pos[1]); ast_destroy(ast_pos[2]);
  ast_destroy(ast_sel); ast_destroy(ast_wc); ast_destroy(ast_wfkp);
  ast_destroy(ast_nz);
  ascii_arg_destroy(arg, nc); free(col); free(chunk);

  /* Release pre-allocated memory that is not necessary any more. */
  DATA *tmp = realloc(dat, n * sizeof(DATA));
  if (tmp) *data = tmp;
  else *data = dat;

  *num = n;
  *sumw = sumwt;
  *I12 = sumw2;
  *I22 = sumw2n;
  return 0;
}
