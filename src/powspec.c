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

PK *compute_pk(CATA *cata, bool save_out, bool has_randoms, int argc, char *argv[]) {
  //printf("The following arguments were passed to main():\n");
  //printf("argnum \t value \n");
  //for (int i = 0; i<argc; i++) printf("%d \t %s \n", i, argv[i]);
  //printf("\n");

  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return NULL;
  }

  conf->ndata = cata->num; // Trigger XPk computation even if not saved
  if (cata->rand != NULL){
    conf->has_randoms = has_randoms;
  }
 
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

  //for (int i = 0; i < pk->nbin;i++){
  //  for (int k = 0; k < cata->num; k++){
  //    printf("%lf ", pk->pl[k][1][i]);
  //  }
  //  printf("\n");
  //  
  //}


   if (save_out && save_res(conf, cata, mesh, pk)) {
     printf(FMT_FAIL);
     P_EXT("failed to write the output to file.\n");
     conf_destroy(conf);
     mesh_destroy(mesh); powspec_destroy(pk);
     return NULL;

   }

  conf_destroy(conf);
  mesh_destroy(mesh);
  
  return pk;
  
}



PK *compute_pk_mesh(double *raw_mesh, bool save_out, int argc, char *argv[]) {
  //printf("The following arguments were passed to main():\n");
  //printf("argnum \t value \n");
  //for (int i = 0; i<argc; i++) printf("%d \t %s \n", i, argv[i]);
  //printf("\n");

  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return NULL;
  }
  conf->has_randoms = false;
  conf->ndata = 1; // Trigger XPk computation even if not saved
  
  if (!conf->issim){
    printf(FMT_FAIL);
    P_EXT("Power spectrum from mesh only available for simulations.\n");
    return NULL;
  }
  conf-> intlace = false;
  printf("Disabling interlacing for mesh calculations.\n");
  if (conf->intlace){
    printf(FMT_FAIL);
    P_EXT("Power spectrum from mesh does not support interlacing.\n");
    return NULL;
  }
 
  CATA *cata = cata_init(1); //store metadata such as shot noise
  cata->num = 1;
  if (!(cata->shot = calloc(1, sizeof(double)))) {
    cata_destroy(cata); return NULL;
  }
  if (!(cata->norm = calloc(1, sizeof(double)))) {
    cata_destroy(cata); return NULL;
  }

  
  MESH *mesh;
  if (!(mesh = genr_mesh_from_mesh(conf, raw_mesh, cata))) {
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



   if (save_out && save_res(conf, cata, mesh, pk)) {
     printf(FMT_FAIL);
     P_EXT("failed to write the output to file.\n");
     conf_destroy(conf);
     mesh_destroy(mesh); powspec_destroy(pk);
     return NULL;

   }

  conf_destroy(conf);
  mesh_destroy(mesh);
  cata_destroy(cata);
  return pk;
  
}
 