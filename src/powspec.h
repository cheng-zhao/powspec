#ifndef _POWSPEC_H_
#define _POWSPEC_H_


#include "read_cata.h"
#include "multipole.h"

void free_pk_array(double *pk_array);
double *compute_pk(CATA *cata, int *nkbin, int argc, char *argv[]);

#endif