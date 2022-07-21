/*******************************************************************************
* powspec.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "read_cata.h"
#include "cnvt_coord.h"
#include "genr_mesh.h"
#include "multipole.h"
#include "save_res.h"
#include <stdio.h>
#include <stdlib.h>

PK *compute_pk(CATA *cata, int *nkbin, int argc, char *argv[]) {
  printf("The following arguments were passed to main():\n");
  printf("argnum \t value \n");
  for (int i = 0; i<argc; i++) printf("%d \t %s \n", i, argv[i]);
  printf("\n");

  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return NULL;
  }

  // CATA *cata;
  // if (!(cata = read_cata(conf))) {
  //   printf(FMT_FAIL);
  //   P_EXT("failed to read the catalogs.\n");
  //   conf_destroy(conf);
  //   return POWSPEC_ERR_CATA;
  // }

  if (cnvt_coord(conf, cata)) {
    printf(FMT_FAIL);
    P_EXT("failed to convert coordinates\n");
    conf_destroy(conf); cata_destroy(cata);
    return NULL;
  }

  MESH *mesh;
  if (!(mesh = genr_mesh(conf, cata))) {
    printf(FMT_FAIL);
    P_EXT("failed to generate the density fields\n");
    conf_destroy(conf); cata_destroy(cata);
    return NULL;
  }

  PK *pk;
  if (!(pk = powspec(conf, cata, mesh))) {
    printf(FMT_FAIL);
    P_EXT("failed to compute the power spectra\n");
    conf_destroy(conf); cata_destroy(cata); mesh_destroy(mesh);
    return NULL;
  }

  // if (save_res(conf, cata, mesh, pk)) {
  //   printf(FMT_FAIL);
  //   P_EXT("failed to write the output to file.\n");
  //   conf_destroy(conf);
  //   mesh_destroy(mesh); powspec_destroy(pk);
  //   return POWSPEC_ERR_SAVE;
  // }

  conf_destroy(conf);
  mesh_destroy(mesh);
  
  return pk;
  
}
 // *nkbin = pk->nbin;
 // int nbin = pk->nbin;//

 // double *pk_array = calloc(4 * nbin, sizeof(size_t));
 // 
 // for (int i = 0; i < nbin; i++) {
 //     pk_array[i] = pk->k[i];
 //     pk_array[nbin + i] = pk->pl[0][0][i];
 //     pk_array[2 * nbin + i] = pk->pl[0][1][i];
 //     pk_array[3 * nbin + i] = pk->pl[0][2][i];     
 // }//

 // powspec_destroy(pk);
 // return pk_array;
//}