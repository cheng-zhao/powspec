/*******************************************************************************
* read_cata.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include "read_cata.h"
#include <stdlib.h>

/*============================================================================*\
                  Function for initialising the data structure
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing information of the catalogs.
Arguments:
  * `ndata`:    number of data catalogs.
Return:
  Address of the structure.
******************************************************************************/
static CATA *cata_init(const int ndata) {
  if (ndata < 0) return NULL;
  CATA *cat = malloc(sizeof *cat);
  if (!cat) return NULL;

  cat->data = cat->rand = NULL;
  cat->ndata = cat->nrand = NULL;
  cat->wdata = cat->wrand = NULL;
  cat->alpha = cat->shot = cat->norm = NULL;

  cat->num = ndata;
  if (!(cat->data = malloc(ndata * sizeof(DATA)))) {
    free(cat); return NULL;
  }
  if (!(cat->rand = malloc(ndata * sizeof(DATA)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->ndata = calloc(ndata, sizeof(size_t)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->nrand = calloc(ndata, sizeof(size_t)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->wdata = calloc(ndata, sizeof(double)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->wrand = calloc(ndata, sizeof(double)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->alpha = calloc(ndata, sizeof(double)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->shot = calloc(ndata, sizeof(double)))) {
    cata_destroy(cat); return NULL;
  }
  if (!(cat->norm = calloc(ndata, sizeof(double)))) {
    cata_destroy(cat); return NULL;
  }

  for (int i = 0; i < ndata; i++) cat->data[i] = cat->rand[i] = NULL;
  return cat;
}


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
CATA *read_cata(const CONF *conf) {
  printf("Reading catalogs ...");
  if (!conf) {
    P_ERR("configuration parameters not loaded.\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);


  CATA *cat = cata_init(conf->ndata);
  if (!cat) return NULL;

  /* Read all the catalogs. */
  for (int i = 0; i < conf->ndata; i++) {
    double sumw2[2];
    double sumw2n[2];
    sumw2[0] = sumw2[1] = sumw2n[0] = sumw2n[1] = 0;

    /* Read data catalog. */
    if (conf->verbose) {
      if (conf->ndata == 1) printf("  < data catalog >\n");
      else printf(" < data catalog %d >\n", i + 1);
    }

    int ftype = (conf->dftype) ? conf->dftype[i] : DEFAULT_FILE_FORMAT;
    if (ftype == POWSPEC_FFMT_ASCII) {
      long skip = (conf->dskip) ? conf->dskip[i] : DEFAULT_FILE_SKIP;
      char cmt = (conf->dcmt) ? conf->dcmt[i] : DEFAULT_FILE_CMT;
      char *sel = (conf->dsel) ? conf->dsel[i] : NULL;
      char *wc = (conf->dwcomp) ? conf->dwcomp[i] : NULL;
      char *wfkp = (!conf->issim && conf->dwfkp) ? conf->dwfkp[i] : NULL;
      char *nz = (!conf->issim && conf->dnz) ? conf->dnz[i] : NULL;
      if (read_ascii_data(conf->dfname[i], skip, cmt, conf->dfmtr[i],
          conf->dpos + i * 3, wc, wfkp, nz, sel, conf->issim, cat->data + i,
          cat->ndata + i, cat->wdata + i, sumw2, sumw2n, conf->verbose)) {
        cata_destroy(cat);
        return NULL;
      }
    }
#ifdef WITH_CFITSIO
    else if (ftype == POWSPEC_FFMT_FITS) {
      P_ERR("FITS format not implemented yet.\n");      /* TODO */
      cata_destroy(cat);
      return NULL;
    }
#endif
    /* Read random catalog if necessary. */
    if (!conf->issim) {
      if (conf->verbose) {
        if (conf->ndata == 1) printf("\n  < random catalog >\n");
        else printf("\n < random catalog %d >\n", i + 1);
      }

      ftype = (conf->rftype) ? conf->rftype[i] : DEFAULT_FILE_FORMAT;
      if (ftype == POWSPEC_FFMT_ASCII) {
        long skip = (conf->rskip) ? conf->rskip[i] : DEFAULT_FILE_SKIP;
        char cmt = (conf->rcmt) ? conf->rcmt[i] : DEFAULT_FILE_CMT;
        char *sel = (conf->rsel) ? conf->rsel[i] : NULL;
        char *wc = (conf->rwcomp) ? conf->rwcomp[i] : NULL;
        char *wfkp = (conf->rwfkp) ? conf->rwfkp[i] : NULL;
        char *nz = (conf->rnz) ? conf->rnz[i] : NULL;
        if (read_ascii_data(conf->rfname[i], skip, cmt, conf->rfmtr[i],
            conf->rpos + i * 3, wc, wfkp, nz, sel, conf->issim, cat->rand + i,
            cat->nrand + i, cat->wrand + i, sumw2 + 1, sumw2n + 1,
            conf->verbose)) {
          cata_destroy(cat);
          return NULL;
        }
      }
#ifdef WITH_CFITSIO
      else if (ftype == POWSPEC_FFMT_FITS) {
        P_ERR("FITS format not implemented yet.\n");      /* TODO */
        cata_destroy(cat);
        return NULL;
      }
#endif

      /* Compute the weights if necessary. */
      if (cat->wdata[i] == 0 || sumw2[0] == 0) {
        P_ERR("invalid completeness or FKP weights in the data catalog.\n");
        cata_destroy(cat); return NULL;
      }
      if (cat->wrand[i] == 0 || sumw2[1] == 0) {
        P_ERR("invalid completeness or FKP weights in the random catalog.\n");
        cata_destroy(cat); return NULL;
      }
      cat->alpha[i] = cat->wdata[i] / cat->wrand[i];
      /* Shot noise from both data and random */
      cat->shot[i] = sumw2[0] + cat->alpha[i] * cat->alpha[i] * sumw2[1];
      /* Normalise using the random. */
      if (sumw2n[0] == 0) cat->norm[i] = cat->alpha[i] * sumw2n[1];
      /* Normalise using the data. */
      else if (sumw2n[1] == 0) cat->norm[i] = sumw2n[0];
      /* Check consistency (TODO) and Normalise using the random. */
      else {
        cat->norm[i] = cat->alpha[i] * sumw2n[1];
      }
    }
  }

  printf(FMT_DONE);
  return cat;
}

/******************************************************************************
Function `cata_destroy`:
  Release memory allocated for the catalogs.
Arguments:
  * `cat`:      the structure for storing information from the catalogs.
******************************************************************************/
void cata_destroy(CATA *cat) {
  if (!cat) return;
  for (int i = 0; i < cat->num; i++) {
    if (cat->data && cat->data[i]) free(cat->data[i]);
    if (cat->rand && cat->rand[i]) free(cat->rand[i]);
  }
  if (cat->data) free(cat->data);
  if (cat->rand) free(cat->rand);
  if (cat->ndata) free(cat->ndata);
  if (cat->nrand) free(cat->nrand);
  if (cat->wdata) free(cat->wdata);
  if (cat->wrand) free(cat->wrand);
  if (cat->alpha) free(cat->alpha);
  if (cat->shot) free(cat->shot);
  if (cat->norm) free(cat->norm);
  free(cat);
}

