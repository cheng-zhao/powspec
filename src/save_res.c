/*******************************************************************************
* save_res.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "read_cata.h"
#include "genr_mesh.h"
#include "multipole.h"
#include "write_file.h"

/*============================================================================*\
                          Interface for saving results
\*============================================================================*/

/******************************************************************************
Function `save_res`:
  Write the power spectra to files.
Arguments:
  * `conf`:     structure for storing all configurations;
  * `cat`:      structure for storing information of the catalogs;
  * `mesh`:     structure for storing information of the meshes;
  * `pk`:       structure for storing the power spectra.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_res(const CONF *conf, const CATA *cat, const MESH *mesh,
    const PK *pk) {
  printf("Saving outputs ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  OFILE *ofile = output_init();
  if (!ofile) return POWSPEC_ERR_FILE;

  for (int n = 0; n < cat->num; n++) {
    if (!conf->isauto[n]) continue;

    /* Save the auto power spectra. */
    if (output_newfile(ofile, conf->oauto[n])) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }

    /* Write header. */
    if (conf->oheader) {
      /* Meta data of the catalogs. */
      if (output_writeline(ofile, "%c Data catalog: %zu objects, total weight: "
          OFMT_DBL "\n", POWSPEC_SAVE_COMMENT, cat->ndata[n],
          cat->wdata[n])) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }

      if (!conf->issim) {
        if (output_writeline(ofile, "%c Random catalog: %zu objects, "
            "total weight: " OFMT_DBL "\n", POWSPEC_SAVE_COMMENT,
            cat->nrand[n], cat->wrand[n])) {
          output_destroy(ofile); return POWSPEC_ERR_FILE;
        }
      }

      /* Meta data of the meshes. */
      if (output_writeline(ofile, "%c Box size: [%lg, %lg, %lg]\n%c "
          "Box boundaries: [[%lg,%lg], [%lg,%lg], [%lg,%lg]]\n",
          POWSPEC_SAVE_COMMENT, mesh->bsize[0], mesh->bsize[1], mesh->bsize[2],
          POWSPEC_SAVE_COMMENT, mesh->min[0], mesh->min[0] + mesh->bsize[0],
          mesh->min[1], mesh->min[1] + mesh->bsize[1],
          mesh->min[2], mesh->min[2] + mesh->bsize[2])) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }

      if (output_writeline(ofile, "%c Grid size: %d , assignment scheme: %s , "
          "grid interlacing: %s\n", POWSPEC_SAVE_COMMENT, mesh->Ng,
          powspec_assign_names[mesh->assign],
          mesh->intlace ? "enabled" : "disabled")) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }

      /* Meta data of the power spectra evaluation. */
      if (output_writeline(ofile, "%c Shot noise: " OFMT_DBL " , normalisation: "
          OFMT_DBL " \n", POWSPEC_SAVE_COMMENT,
          conf->issim ? cat->shot[n] : cat->shot[n] / cat->norm[n],
          cat->norm[n])) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }
    }

    /* Write column indicators. */
    if (output_writeline(ofile, "%c kcen(1) kmin(2) kmax(3) kavg(4) nmod(5)",
        POWSPEC_SAVE_COMMENT)) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }
    for (int l = 0; l < pk->nl; l++) {
      if (output_writeline(ofile, " P_%d(%d)", pk->poles[l], l + 6)) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }
    }
    if (output_writeline(ofile, "\n")) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }

    /* Write columns. */
    for (int i = 0; i < pk->nbin; i++) {
      if (output_writeline(ofile, OFMT_DBL " " OFMT_DBL " " OFMT_DBL " "
          OFMT_DBL " %zu", pk->k[i], pk->kedge[i], pk->kedge[i + 1], pk->km[i],
          pk->cnt[i])) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }
      for (int l = 0; l < pk->nl; l++) {
        if (output_writeline(ofile, " " OFMT_DBL, pk->pl[n][l][i])) {
          output_destroy(ofile); return POWSPEC_ERR_FILE;
        }
      }
      if (output_writeline(ofile, "\n")) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }
    }

    if (conf->verbose) {
      if (cat->num == 2) {
        printf("  Auto power spectra for catalog %d saved to file: `%s'\n",
          n + 1, conf->oauto[n]);
      }
      else
        printf("  Auto power spectra saved to file: `%s'\n", conf->oauto[n]);
    }
  }

  if (!conf->iscross) {
    output_destroy(ofile);
    printf(FMT_DONE);
    return 0;
  }

  /* Save the cross power spectra. */
  if (output_newfile(ofile, conf->ocross)) {
    output_destroy(ofile); return POWSPEC_ERR_FILE;
  }

  /* Write header. */
  if (conf->oheader) {
    /* Meta data of the catalogs. */
    if (output_writeline(ofile, "%c Data catalog 1: %zu objects, total weight: "
        OFMT_DBL "\n%c Data catalog 2: %zu objects, total weight: " OFMT_DBL
        "\n", POWSPEC_SAVE_COMMENT, cat->ndata[0], cat->wdata[0],
        POWSPEC_SAVE_COMMENT, cat->ndata[1], cat->wdata[1])) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }

    if (!conf->issim) {
      if (output_writeline(ofile, "%c Random catalog 1: %zu objects, "
          "total weight: " OFMT_DBL "\n%c Random catalog 2: %zu objects, "
          "total weight: " OFMT_DBL "\n", POWSPEC_SAVE_COMMENT, cat->nrand[0],
          cat->wrand[0], POWSPEC_SAVE_COMMENT, cat->nrand[1], cat->wrand[1])) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }
    }

    /* Meta data of the meshes. */
    if (output_writeline(ofile, "%c Box size: [%lg, %lg, %lg]\n%c "
        "Box boundaries: [[%lg,%lg], [%lg,%lg], [%lg,%lg]]\n",
        POWSPEC_SAVE_COMMENT, mesh->bsize[0], mesh->bsize[1], mesh->bsize[2],
        POWSPEC_SAVE_COMMENT, mesh->min[0], mesh->min[0] + mesh->bsize[0],
        mesh->min[1], mesh->min[1] + mesh->bsize[1],
        mesh->min[2], mesh->min[2] + mesh->bsize[2])) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }

    if (output_writeline(ofile, "%c Grid size: %d , assignment scheme: %s , "
        "grid interlacing: %s\n", POWSPEC_SAVE_COMMENT, mesh->Ng,
        powspec_assign_names[mesh->assign],
        mesh->intlace ? "enabled" : "disabled")) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }

    /* Meta data of the power spectra evaluation. */
    if (output_writeline(ofile, "%c Normalisation factor 1: " OFMT_DBL
        " , Normalisation factor 2: " OFMT_DBL "\n", POWSPEC_SAVE_COMMENT,
        cat->norm[0], cat->norm[1])) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }
  }

  /* Write column indicators. */
  if (output_writeline(ofile, "%c kcen(1) kmin(2) kmax(3) kavg(4) nmod(5)",
      POWSPEC_SAVE_COMMENT)) {
    output_destroy(ofile); return POWSPEC_ERR_FILE;
  }
  for (int l = 0; l < pk->nl; l++) {
    if (output_writeline(ofile, " P_%d(%d)", pk->poles[l], l + 6)) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }
  }
  if (output_writeline(ofile, "\n")) {
    output_destroy(ofile); return POWSPEC_ERR_FILE;
  }

  /* Write columns. */
  for (int i = 0; i < pk->nbin; i++) {
    if (output_writeline(ofile, OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL
        " %zu", pk->k[i], pk->kedge[i], pk->kedge[i + 1], pk->km[i],
        pk->cnt[i])) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }
    for (int l = 0; l < pk->nl; l++) {
      if (output_writeline(ofile, " " OFMT_DBL, pk->xpl[l][i])) {
        output_destroy(ofile); return POWSPEC_ERR_FILE;
      }
    }
    if (output_writeline(ofile, "\n")) {
      output_destroy(ofile); return POWSPEC_ERR_FILE;
    }
  }

  output_destroy(ofile);

  if (conf->verbose)
    printf("  Cross power spectra saved to file: `%s'\n", conf->ocross);

  printf(FMT_DONE);
  return 0;
}

