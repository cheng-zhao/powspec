/*******************************************************************************
* read_file.c: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "read_file.h"
#include <stdlib.h>

const char *powspec_ffmt_names[] = {
  "ASCII",
  "FITS"
};

