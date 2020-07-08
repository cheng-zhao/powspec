/*******************************************************************************
* save_res.h: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

*******************************************************************************/

#ifndef _SAVE_RES_H_
#define _SAVE_RES_H_

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
int save_res(const CONF *conf, const CATA *cat, const MESH *mesh, const PK *pk);

#endif
