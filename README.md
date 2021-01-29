# powspec

![Codacy grade](https://img.shields.io/codacy/grade/5a6835fc4cee49d4993380ad1a45ad69.svg)

## Table of Contents

-   [Introduction](#introduction)
-   [Algorithm](#algorithm)
    -   [Simulation box](#simulation-box)
        -   [Particle assignment](#particle-assignment)
        -   [Grid interlacing](#grid-interlacing)
        -   [Multipole evaluation](#multipole-evaluation)
    -   [Survey data](#survey-data)
        -   [Coordinate conversion](#coordinate-conversion)
        -   [Box configuration](#box-configuration)
        -   [Weighted density field](#weighted-density-field)
        -   [Multipole estimator](#multipole-estimator)
-   [Compilation](#compilation)
-   [Configuration parameters and command line options](#configuration-parameters-and-command-line-options)
    -   [Specifications of the input catalogues](#specifications-of-the-input-catalogues)
    -   [Fiducial cosmology for coordinate conversion](#fiducial-cosmology-for-coordinate-conversion)
    -   [Configurations for power spectra evaluation](#configurations-for-power-spectra-evaluation)
    -   [Settings for the output](#settings-for-the-output)
-   [Acknowledgements](#acknowledgements)
-   [References](#references)

## Introduction

This is a C program for computing auto and cross power spectrum multipoles *P*<sub>&ell;</sub>&thinsp;(*k*) given the 3-D positions of tracers. It is applicable for both simulations in periodic boxes, and survey-like data with arbitrary geometry and selection effects.

It is designed in the hope of being (both time and space) efficient, portable, and user-friendly. To this end, various operations are provided for pre-processing the data on-the-fly and in parallel, with a least number of external dependencies. Moreover, little programming knowledge is required for the usage of the code, and the user is only ask for a form-like configuration file.

This program is compliant with the ISO C99 and IEEE POSIX.1-2008 standards. Parallelisation can be enabled with [OpenMP](https://www.openmp.org). Thus it is compatible with most of the modern C compilers and operating systems. It is written by Cheng Zhao (&#36213;&#25104;). And as a whole, it is released under the [GPLv3 license](LICENSE.txt), due to its reliance on the [FFTW library](http://www.fftw.org), though many of the source files are distributed under the [MIT license](LICENSE_MIT.txt), as indicated in their header lines.

If you use this tool in research that results in publications, please cite the following paper:

> Zhao et al. 2020, [arXiv:2007.08997](https://ui.adsabs.harvard.edu/abs/2020arXiv200708997Z/abstract)

<sub>[\[TOC\]](#table-of-contents)</sub>

## Algorithm

### Simulation box

For simulation boxes, the data is assumed to be in a regular cuboid, with periodic boundary conditions on all directions. The side length of the box (*L*<sub>box</sub>) has then to be supplied by the user, and the 3-D coordinates of the input data (*x*<sub>0</sub>, *x*<sub>1</sub>, *x*<sub>2</sub>) must satisfy

0 &le; *x*<sub>*i*</sub> &lt; *L*<sub>box</sub> , &emsp;*i* = 0, 1, 2.

Otherwise, the coordinates of the input catalogue have to be pre-processed, see [Specifications of the input catalogues](#specifications-of-the-input-catalogues) for details.

#### Particle assignment

For a cubic simulation box, tracers are distributed to a regular 3-D mesh for Fast Fourier Transform (FFT), with optionally the following grid assignment schemes <sup>[\[1\]](#ref1)</sup>:

-   Nearest-Grid-Point (NGP):
    -   *w* = 1,&emsp;if *d* &lt; 1/2;
    -   *w* = 0,&emsp;otherwise;
-   Could-In-Cell (CIC):
    -   *w* = 1 &minus; *d*,&emsp;if *d* &lt; 1;
    -   *w* = 0,&emsp;otherwise;
-   Triangular Shaped Cloud (TSC):
    -   *w* = 3/4 &minus; *d*<sup>2</sup>,&emsp;if |*d*| &lt; 1/2;
    -   *w* = (3/2 &minus; *d*)<sup>2</sup>/2,&emsp;if 1/2 &le; *d* &lt; 1;
    -   *w* = 0,&emsp;otherwise;
-   Piecewise Cubic Spline (PCS):
    -   *w* = (4 &minus; 6&thinsp;*d*<sup>2</sup> + 3&thinsp;*d*<sup>3</sup>)/6,&emsp;if |*d*| &lt; 1;
    -   *w* = (2 &minus; *d*)<sup>3</sup>/6,&emsp;if 1 &le; |*d*| &lt; 2;
    -   *w* = 0,&emsp;otherwise;

where *w* denotes the number assigned to a grid point for each dimension, given a particle with distance component (*d* &middot; *H*) to the grid point, with *H* being the side length of one cell.

The window effect induced by the particle assignment scheme is then corrected following [\[2\]](#ref2).

<sub>[\[TOC\]](#table-of-contents)</sub>

#### Grid interlacing

To reduce the "alias" effects of particle assignments, the grid interlacing technique<sup>[\[1\]](#ref1)</sup> is implemented. It generates an addition real-space density field with a shift of half the cell size (*H*/2) on all directions. The two densities fields are combined in Fourier space, and the combined field is then used for power spectrum evaluation.

Since the density field is real, we use the real-to-complex routines of the FFTW library for Fourier transforms, to increase both time and space efficiencies.

<sub>[\[TOC\]](#table-of-contents)</sub>

#### Multipole evaluation

The power spectrum is evaluated via mode counts of the Fourier space density field *&delta;*&thinsp;(***k***). And the multipoles are given by

*P*<sub>&ell;</sub>&thinsp;(*k*) = &langle; |&thinsp;*&delta;*&thinsp;(***k***) *&delta;*&thinsp;(&minus;***k***)&thinsp;| *L*<sub>&ell;</sub>&thinsp;(*&mu;*) &rangle; .

Here, the average &langle;&bull;&rangle; is taken over all grid points in Fourier space with *k*<sub>low</sub> &le; |***k***| &lt; *k*<sub>high</sub>, where *k*<sub>low</sub> and *k*<sub>high</sub> indicate the lower and upper limit of the *k* bin, respectively. And *L*<sub>&ell;</sub> indicates the Legendre polynomial, with

*&mu;* = (***k*** &middot; ***l***)/|***k***| ,

where the unit line-of-sight vector ***l*** is supplied by the user.

<sub>[\[TOC\]](#table-of-contents)</sub>

### Survey-like data

Due to the non-trivial geometry and completeness of survey-like data, random catalogues that encode the same selection functions as the data catalogues are required. And typically various weights are supplied for corrections of the incompleteness of the data. In addition, observational coordinates &mdash; Right Ascension (RA), Declination (Dec), and redshift (*z*) &mdash; are usually provided, and has to be converted to the comoving coordinates *x*<sub>0</sub>, *x*<sub>1</sub>, and *x*<sub>2</sub>.

#### Coordinate conversion

The key of the coordinate conversion process is the relationship between redshift *z* and radial comoving distance *r*. Once this relationship is evaluated, the comoving coordinates are then simply given by

*x*<sub>0</sub> = *r* cos&thinsp;(Dec) cos&thinsp;(RA) ,

*x*<sub>1</sub> = *r* cos&thinsp;(Dec) sin&thinsp;(RA) ,

*x*<sub>2</sub> = *r* sin&thinsp;(Dec) .

This program implements the conversion from *z* to *r* within the framework of a *w*CDM cosmology:

*r* = &int;<sub><sub>0</sub></sub><sup><sup>*z*</sup></sup> *c* d&thinsp;*z*&prime;/\[*H*<sub>0</sub> *E*&thinsp;(*z*&prime;)\] ,

where *H*<sub>0</sub> indicates the Hubble parameter at present (*z* = 0), and *E*&thinsp;(*z*) denotes the reduced Hubble parameter:

*E*<sup>2</sup>&thinsp;(*z*) = &Omega;<sub>m</sub>&thinsp;(1+*z*)<sup>3</sup> + &Omega;<sub>*k*</sub>&thinsp;(1+*z*)<sup>2</sup> + &Omega;<sub>&Lambda;</sub>&thinsp;(1+*z*)<sup>3(1+*w*)</sup> .

Here, &Omega;<sub>m</sub>, &Omega;<sub>*k*</sub>, and &Omega;<sub>&Lambda;</sub> indicate the density parameter for matter, curvature, and dark energy at *z* = 0, respectively. And *w* denotes the dark energy equation of state. In particular, only &Omega;<sub>m</sub>, &Omega;<sub>&Lambda;</sub>, and *w* are supplied by the user, as

&Omega;<sub>*k*</sub> = 1 &minus; &Omega;<sub>m</sub> &minus; &Omega;<sub>&Lambda;</sub> .

Furthermore, to ensure *E*<sup>2</sup>&thinsp;(*z*) &gt; 0, the following condition has to be satisfied:

3*w*&Omega;<sub>m</sub><sup>(3*w*+1)/(3*w*)</sup>&Omega;<sub>&Lambda;</sub> &lt; &Omega;<sub>*k*</sub>\[&minus;(3*w*+1)&Omega;<sub>&Lambda;</sub>\]<sup>(3*w*+1)/(3*w*)</sup> .

The integration for the *z* to *r* conversion is performed using the [*Legendre-Gauss* Quadrature](https://mathworld.wolfram.com/Legendre-GaussQuadrature.html) algorithm. The program uniformly samples 128 (customisable in [define.h](src/define.h#L94)) redshift values in the redshift range of the input catalogues, and checks the maximum order needed for the desired integration precision. This order is then used for the coordinate conversion of all the input objects.

For non-*w*CDM cosmologies, the program itself cannot convert *z* to *r*, but it is able to interpolate a table of (*z*, *r*) pairs for coordinate conversions, where the unit of the radial comoving distance has to be Mpc/*h*. And the interpolation is performed using a cubic-spline algorithm<sup>[\[3\]](#ref3)</sup>.

<sub>[\[TOC\]](#table-of-contents)</sub>

#### Box configuration

To sample the tracer density field on grids, a cuboid that is large enough to include all the data needs to be specified. The side lengths of the cuboid can either be supplied by the user, or determined by the program automatically. For user-specified side lengths, the following condition has to be fulfilled:

*x*<sub>*i*, max</sub> &minus; *x*<sub>*i*, min</sub> &ge; *L*<sub>box, *i*</sub> , &emsp; *i* = 0, 1, 2.

If the box size is not set, adaptive side lengths are evaluated as

*L*<sub>box, *i*</sub> = (*x*<sub>*i*, max</sub> &minus; *x*<sub>*i*, min</sub>) &middot; (1 + *f*<sub>pad, *i*</sub>) ,

where *f*<sub>pad, *i*</sub> indicates the user-supplied padding fraction of the box.

Once *L*<sub>box</sub> is decided, the input catalogues are placed at the centre of the box, and then the [particle assignment](#particle-assignment) scheme and corresponding corrections are applied in the same way as for simulation boxes, to generate the density fields.

<sub>[\[TOC\]](#table-of-contents)</sub>

#### Weighted density field

Given the data and random catalogues, the weighted density field is given by

*F*&thinsp;(***r***) = *w*<sub>FKP</sub>&thinsp;(***r***) &middot; \[*n*<sub>d</sub>&thinsp;(***r***) &minus; *&alpha;*&thinsp;*n*<sub>r</sub>&thinsp;(***r***)\] ,

where *w*<sub>FKP</sub> is the so-called FKP weight, for reducing the variance of the power spectra<sup>[\[4\]](#ref4)</sup>. *n*<sub>d</sub>&thinsp;(***r***) and *n*<sub>r</sub>&thinsp;(***r***) are the density field for the data and random catalogues respectively. And *&alpha;* is the weighted ratio between the data and random samples:

*&alpha;* = &sum;<sub>d</sub> *w*<sub>c, d</sub> / &sum;<sub>r</sub> *w*<sub>c, r</sub> .

Here, *w*<sub>c, d</sub> and *w*<sub>c, r</sub> indicate the user-supplied completeness weights for the samples, and the sums are taken over the corresponding catalogues.

The weighted density field is then Fourier transformed using FFTW. Moreover, if [grid interlacing](#grid-interlacing) is enabled, the combined Fourier space density field is further converted back to configuration space, using the complex-to-real routines, for power spectrum multipole estimations.

<sub>[\[TOC\]](#table-of-contents)</sub>

#### Multipole estimator

The power spectrum multipoles are estimated by<sup>[\[5\]](#ref5)[\[6\]](#ref6)</sup>

*P*<sub>&ell;</sub>&thinsp;(*k*) = \[&thinsp;(2&ell;+1) &langle;*F*<sub>0</sub>&thinsp;(***k***) *F*<sub>&ell;</sub>&thinsp;(***k***)&rangle; &minus; (1+*&alpha;*)&thinsp;*I*<sub>12</sub>&thinsp;\] / *I*<sub>22</sub> ,

where

*I*<sub>*ab*</sub> = &int;&thinsp;d<sup>3</sup>&thinsp;*r* *&ntilde;*<sup>*a*</sup>&thinsp;(***r***) *w*<sup>*b*</sup><sub>FKP</sub>&thinsp;(***r***) ,

*F*<sub>&ell;</sub>&thinsp;(***k***) = &int;&thinsp;d<sup>3</sup>&thinsp;*r* *F*&thinsp;(***r***) exp&thinsp;(*i*&thinsp;***k***&thinsp;&middot;&thinsp;***r***) *L*<sub>&ell;</sub>&thinsp;\[&thinsp;***k***&thinsp;&middot;&thinsp;***r***/(|***k***||***r***|)&thinsp;\] .

And *&ntilde;* denotes the comoving number density of the tracers.

Furthermore, the Legendre polynomials *L*<sub>&ell;</sub> can be decomposed into products of real-form spherical harmonics *Y*<sub>&ell;*m*</sub>&thinsp;, the Fourier space density field multipoles are then<sup>[\[7\]](#ref7)</sup>

*F*<sub>&ell;</sub>&thinsp;(***k***) = 4&pi;/(2&ell;+1) &sum;<sub>*m*</sub> *Y*<sub>&ell;*m*</sub>&thinsp;(***k***/|***k***|) &int;&thinsp;d<sup>3</sup>&thinsp;*r* *F*&thinsp;(***r***) *Y*<sup>*</sup><sub>&ell;*m*</sub>&thinsp;(***r***/|***r***|) exp&thinsp;(*i*&thinsp;***k***&thinsp;&middot;&thinsp;***r***) .

Thus, for each of the multipole, only (2&ell;+1) FFTs are required.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Compilation

The dependencies of this program are listed below:

| Library                                         | Mandatory | Compilation flags<sup>[*](#tab1)</sup>                                           |
|:-----------------------------------------------:|:---------:|:---------------------------------------------------------------------------------|
| [FFTW](http://www.fftw.org)                     | &check;   | `-lfftw3`                                                                        |
| [OpenMP](https://www.openmp.org)                | &cross;   | `-DOMP -fopenmp -lfftw3_omp` (gcc, clang)<br />`-DOMP -openmp -lfftw3_omp` (icc) |
| [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio) | &cross;   | `-DWITH_CFITSIO -lcfitsio`                                                       |

<span id="tab1">*: paths to the header and library files are not included (e.g., `-I/path/to/include -L/path/to/lib`).</span>

The density fields are stored as double-precision floating point numbers by default, with FFTs performed in double precision as well. If the memory consumption is an issue, the program can also be compiled with the `-DSINGLE_PREC` flag, to enabled single-precision density fields and FFTs. Note that in this case FFTW has to be compiled with single precision too (see the [FFTW documentation](http://www.fftw.org/fftw3_doc/Installation-and-Customization.html) for details).

Once the mandatory dependencies are installed, and the corresponding compilation flags are set in the [`Makefile`](Makefile), this program can be compiled with a C compiler that supports the C99 standard, by simply the command

```console
make
```

By default, an executable `POWSPEC` is created on success.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Configuration parameters and command line options

Configuration files and command line options for this program are parsed using the [libcfg](https://github.com/cheng-zhao/libcfg) library. The default name of the configuration file is [`powspec.conf`](powspec.conf), which is customisable in [`define.h`](src/define.h#L49). Custom configuration files can also be passed to the program using the `-c` or `--conf` command line options.

A list of the supported command line options can be displaced using the `-h` or `--help` flags, and a template configuration file is printed with the `-t` or `--template` flags. Please consult [libcfg](https://github.com/cheng-zhao/libcfg) for the formats of the configuration files and command line options.

### Specifications of the input catalogues

#### `DATA_CATALOG`, `RAND_CATALOG`

Filename of the input data and random catalogues. They can be either strings or 2-element string arrays. If the cross power spectrum is required, then the data and random catalogues must be string arrays. For simulation boxes, `RAND_CATALOG` is not mandatory.

*Examples*

```nginx
DATA_CATALOG = input_galaxy_catalog.txt
DATA_CATALOG = [input_galaxy_catalog.txt]
DATA_CATALOG = [ galaxy_catalog_1.txt, galaxy_catalog_2.txt ]
DATA_CATALOG = [ galaxy_catalog_1.txt, \
                 galaxy_catalog_2.txt ]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_FORMAT`, `RAND_FORMAT`

Format of the input catalogs. They must be integers, or 2-element integer arrays, depending on the dimension of `DATA_CATALOG` and `RAND_CATALOG` respectively. The allowed values are:

-   `0`: for ASCII format text files (default);
-   `1`: for FITS tables.

In particular, FITS tables are supported via the CFITSIO library, so the compilation flag `-DWITH_CFITSIO` has to be enabled for reading FITS files.

*Examples*

```nginx
DATA_FORMAT = 0
DATA_FORMAT = [0,1]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_SKIP`, `RAND_SKIP`

Number of lines to be skipped for ASCII format input files. They can be non-negative long integers or long integer arrays, depending on the dimension of `DATA_CATALOG` and `RAND_CATALOG` respectively.

*Examples*

```nginx
DATA_SKIP = 10
DATA_SKIP = [0,5]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_COMMENT`, `RAND_COMMENT`

Indicator of comment lines for ASCII format input files. They can be characters or 2-element character arrays, depending on the dimension of `DATA_CATALOG` and `RAND_CATALOG` respectively. If the first non-whitespace character of a line is the specified character, then the whole line of the input catalogue is omitted. If empty characters (`''`) are supplied, then no line is treated as comments.

*Examples*

```nginx
DATA_COMMENT = '#'
DATA_COMMENT = [ '', '!' ]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_FORMATTER`, `RAND_FORMATTER`

C99-style formatter specifying the format of columns of ASCII format input files. They must be strings or 2-element string arrays, depending on the dimension of `DATA_CATALOG` and `RAND_CATALOG` respectively. Note however that only formatters with the following argument types are supported (see [cppreference.com](https://en.cppreference.com/w/c/io/fscanf) for details):

-   `int *`
-   `long *`
-   `float *`
-   `double *`
-   `char *`

*Examples*

```nginx
DATA_FORMATTER = "%d %ld %f %lf %s"  # for int, long, float, double, and string types
DATA_FORMATTER = "%*d,%10s,%[123]"
        # Column separators are ',';
        # The first column is treated as an integer, but is omitted;
        # The second column is a string with 10 characters;
        # The third column is a string composed of characters '1', '2', and '3'.
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_POSITION`, `RAND_POSITION`

3-D coordinates, in the order of \[*x*<sub>0</sub>, *x*<sub>1</sub>, *x*<sub>2</sub>\] or \[RA, Dec, redshift\], where RA and Dec must be in degrees. They must be 3- or 6-element string arrays. If `DATA_CATALOG` or `RAND_CATALOG` contains only one element, then the corresponding positions must be 3-element arrays. While if there are two elements for `DATA_CATALOG` or `RAND_CATALOG`, there should be 6 elements for the positions.

The strings must be column indicators, or expressions, which are parsed using the [libast](https://github.com/cheng-zhao/libast) library. Columns are indicated by <span><code>${&bull;}</code></span>, where <span><code>&bull;</code></span> must be a number for ASCII format files, and a string for FITS tables. For instance, the 3rd column of an ASCII file can be indicated by `$3`, and the "RA" column of a FITS table can be indicated by `${RA}`. Note that if there are more than one digits for the ASCII column numbers, the number must be enclosed by braces, e.g, `${10}`. And if an ASCII column is omitted via the formatter (e.g. `%*lf`), it is not counted for the column number.

Moreover, expressions are supported for pre-processing the columns, with some basic arithmetic operators, and mathematical functions (see [libast](https://github.com/cheng-zhao/libast) for details).

*Examples*
```nginx
DATA_POSITION = [${RA}*180/3.1415927, ${DEC}, ${Z}, ${RA}, ${DEC}, ${Z}]
DATA_POSITION = [($1+1000)%1000, ($2+1000)%1000, ($3+1000)%1000]
DATA_POSITION = [$1, $2, ($3+$6*(1+0.6)/(100*sqrt(0.31*(1+0.6)**3+0.69))+1000)%1000]
      # The last expression implies real to redshift space conversion
      # with periodic boundary conditions given the box size 1000 Mpc/h,
      # in a FlatLCDM cosmology with Omega_m = 0.31, at redshift 0.6.
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_WT_COMP`, `RAND_WT_COMP`

Completeness weights for the data and randoms (see [Weighted density field](#weighted-density-field)). They can be column indicators or expressions, or the corresponding 2-element string arrays.

*Examples*

```nginx
DATA_WT_COMP = ${WEIGHT_SYSTOT} * ${WEIGHT_NOZ} * ${WEIGHT_CP}
DATA_WT_COMP = [1, $4 * ($5 + $6 - 1)]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_WT_FKP`, `RAND_WT_FKP`

FKP weights for the data and randoms (see [Weighted density field](#weighted-density-field)). They can be column indicators or expressions, or the corresponding 2-element string arrays. They are not used for simulation boxes.

*Examples*

```nginx
DATA_WT_FKP = [$5, $5]
DATA_WT_FKP = 1 / (6000 * ${NZ})
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_NZ`, `RAND_NZ`

Radial comoving number density of the data and randoms. They can be column indicators or expressions, or the corresponding 2-element string arrays. They are used for the normalisation of the power spectra (see [Multipole estimator](#multipole-estimator) for details). By default, `RAND_NZ` is used. And if `RAND_NZ` is unset, then `DATA_NZ` is used for the normalisation. They are not used for simulation boxes.

*Examples*

```nginx
RAND_NZ = 5678000 / 1000**3
RAND_NZ = [${NZ}, ${NZ}]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_SELECTION`, `RAND_SELECTION`

Selection criteria for the input catalogues. They can be column indicators or expressions, or the corresponding 2-element string arrays. Numerical, bitwise, and logical expressions are all supported. Only objects with columns fulfilling the conditions are kept.

*Examples*

```nginx
DATA_SELECTION = isfinite($1) && $3 == "YES" && log($4) < 1.0
DATA_SELECTION = [$5 & 1 != 0, $1 != "no"]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

#### `DATA_CONVERT`, `RAND_CONVERT`

Specify whether coordinate conversion is needed for the data and random catalogues (see [Coordinate conversion](#coordinate-conversion) for details). They must be boolean values or 2-element boolean arrays, depending on the dimension of `DATA_CATALOG` and `RAND_CATALOG` respectively. If the conversion for any of the catalogues is enabled, the parameters specified in the section [Fiducial cosmology for coordinate conversion](#fiducial-cosmology-for-coordinate-conversion) are used.

*Examples*

```nginx
DATA_CONVERT = [T,F]
DATA_CONVERT = [ True, 0 ]  # 0 for false
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### Fiducial cosmology for coordinate conversion

The formulae for the coordinate conversion are detailed in the [Coordinate conversion](#coordinate-conversion) section.

#### `OMEGA_M`

Density parameter of matter at *z* = 0. It must be a double-precision floating point number, in the range (0, 1].

#### `OMEGA_LAMBDA`

Density parameter of dark energy at *z* = 0. It must be a double-precision floating point number, and &ge; 0. By default it is 1 &minus; `OMEGA_M`.

#### `DE_EOS_W`

Dark energy equation of state. It must be a double-precision floating point number, and &le; 1/3. By default it is &minus;1.

#### `CMVDST_ERR`

The tolerance of the integration error. It must be larger than the [machine epsilon](https://en.wikipedia.org/wiki/Machine_epsilon), i.e., around `1e-16`.

#### `Z_CMVDST_CNVT`

Filename of an ASCII table for redshift to radial comoving distance (in Mpc/*h*) conversion. The first two columns of this file have to be redshift and radial comoving distance, respectively. If the columns or units are not appropriate, the file can be passed to the program via command line options and pipe, e.g.

```bash
./POWSPEC --cmvdst-file <(awk '{printf("%lf %lf\n", $3, $4 * 0.676)}' input_cnvt_file.txt)
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### Configurations for power spectra evaluation

#### `CUBIC_SIM`

Indicate whether the input catalogs are from cubic simulation boxes. It must be a boolean value.

*Examples*

```nginx
CUBIC_SUM = T
```

#### `LINE_OF_SIGHT`

A unit vector defining the line of sight for cubic simulation boxes (see [Multipole evaluation](#multipole-evaluation)). By default it is `[0,0,1]`, i.e., the line of sight is along the 3rd Cartesian axis.

#### `BOX_SIZE`

Side length of the box that catalogues are placed in. It can be a double-precision floating-point number, or 3-element double array. In particular, a single double number implies a regular cuboid, or the side length on all directions are identical. Otherwise, box size on different directions are read from the array.

It is mandatory for simulation boxes. And adaptive box sizes are evaluated if it is not supplied for a survey-like data.

*Examples*

```nginx
BOX_SIZE = 1000
BOX_SIZE = [600, 800, 1000]
```

#### `BOX_PAD`

Padding fraction for the adaptive box size. It can be a double-precision floating-point number, or 3-element double array. And it is essentially *f*<sub>pad, *i*</sub> in section [Box configuration](#box-configuration).

#### `GRID_SIZE`

Number of grid cells per box side for the density fields. It must be a positive integer, preferably a power of 2.

#### `PARTICLE_ASSIGN`

Particle assignment scheme (see section [Particle assignment](#particle-assignment)). It must be an integer, and allowed values are:

-   `0`: Nearest-Grid-Point (NGP);
-   `1`: Could-In-Cell (CIC);
-   `2`: Triangular Shaped Cloud (TSC);
-   `3`: Piecewise Cubic Spline (PCS).

#### `GRID_INTERLACE`

A boolean option indicating whether to use interlaced grids (see [Grid interlacing](#grid-interlacing)).

<sub>[\[TOC\]](#table-of-contents)</sub>

### Settings for the output

#### `MULTIPOLE`

Legendre multipoles of the power spectra to be evaluated. It must be an non-negative integer, or integer arrays. The current maximum supported &ell; is `6`.

*Examples*

```nginx
MULTIPOLE = 0
MULTIPOLE = [0,1,2,3,4,5,6]
```

#### `LOG_SCALE`

A boolean value indicating whether to set wave number bins in logarithm scale.

#### `KMIN`

Lower boundary of the first wave number bin. It can be a non-negative double-precision floating-point number. For logarithm scale bins, it must be positive.

#### `KMAX`

Upper boundary of the last wave number bin. It can be a non-negative double-precision floating-point number, and must be larger than `KMIN`. If it is unset, then the Nyquist frequency is used.

#### `BIN_SIZE`

Width of each wave number bin. It must be a double-precision floating-point number. And for logarithm scales, it indicates the base-10 logarithm of the ratio between two adjacent wave number bin edges.

#### `OUTPUT_AUTO`

Name of the output files for auto power spectrum multipoles. It can be a string or 2-element string array. In particular, auto power spectra can be omitted by setting an empty string.

*Examples*
```nginx
OUTPUT_AUTO = output_auto_pk.txt
OUTPUT_AUTO = ["", output_auto_pk2.txt]
        # The auto power spectrum multipoles for the first catalogue are omitted.
```

#### `OUTPUT_CROSS`

String, for the name of the output file for cross power spectrum multipoles.

#### `OUTPUT_HEADER`

A boolean value indicating whether to write extra information as the header of the output files, including number of objects in the input catalogues, configurations of the box and grids, and the shot noise correction and normalisation factors.

#### `OVERWRITE`

An integer value indicating whether to overwrite existing files. Allowed values are

-   `0`: quit the program when an output file exist;
-   postive: force overwriting output files whenever possible;
-   negative: notify at most this number (absolute value) of times, for asking whether overwriting existing files.

#### `VERBOSE`

A boolean value indicating whether to show detailed standard outputs.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Acknowledgements

This program benefits greatly from the open-source [nbodykit](https://nbodykit.readthedocs.io/en/latest/index.html) project<sup>[\[8\]](#ref8)</sup>. And I thank Dr. Chia-Hsun Chuang for helpful discussions on the early-stage development of this program. 

<sub>[\[TOC\]](#table-of-contents)</sub>

## References

<span id="ref1">\[1\]</span> Sefusatti, Crocce, Scoccimarro, Couchman, 2016, [Correcting for the Alias Effect When Measuring the Power Spectrum Using a Fast Fourier Transform](https://doi.org/10.1093/mnras/stw1229), *MNRAS*, 460(4):3624&ndash;3636 ([arXiv:1512.07295](https://arxiv.org/abs/1512.07295))

<span id="ref2">\[2\]</span> Jing, 2005, [Accurate estimators of correlation functions in Fourier space](https://doi.org/10.1086/427087), *ApJ*, 620(2):559&ndash;563 ([arXiv:astro-ph/0409240](https://arxiv.org/abs/astro-ph/0409240))

<span id="ref3">\[3\]</span> Hornbeck, 2020, [Fast Cubic Spline Interpolation](https://arxiv.org/abs/2001.09253), *arXiv e-prints*, 2001.09253 ([source code](https://doi.org/10.5281/zenodo.3611922))

<span id="ref4">\[4\]</span> Feldman, Kaiser, Peacock, 1994, [Power-Spectrum Analysis of Three-dimensional Redshift Surveys](https://doi.org/10.1086/174036), *ApJ*, 426:23 ([arXiv:astro-ph/9304022](https://arxiv.org/abs/astro-ph/9304022))

<span id="ref5">\[5\]</span> Sefusatti, 2005, [Probing fundamental physics with large-scale structure: From galaxy formation to inflation](https://ui.adsabs.harvard.edu/abs/2005PhDT........23S), *PhD Thesis*, New York University)

<span id="ref6">\[6\]</span> Yamamoto, Nakamichi, Kamino, Bassett, Nishioka, 2006, [A Measurement of the Quadrupole Power Spectrum in the Clustering of the 2dF QSO Survey](https://doi.org/10.1093/pasj/58.1.93), *PASJ*, 93:102 ([arXiv:astro-ph/0505115](https://arxiv.org/abs/astro-ph/0505115))

<span id="ref7">\[7\]</span> Hand, Li, Slepian, Seljak, 2017, [An optimal FFT-based anisotropic power spectrum estimator](https://doi.org/10.1088/1475-7516/2017/07/002), *JCAP*, 2017:2 ([arXiv:1704.02357](https://arxiv.org/abs/1704.02357))

<span id="ref8">\[8\]</span> Hand, Feng, Beutler, Li, Modi, Seljak, Slepian, 2018, [nbodykit: An Open-source, Massively Parallel Toolkit for Large-scale Structure](https://doi.org/10.3847/1538-3881/aadae0), *AJ*, 156(4):160 ([arXiv:1712.05834](https://arxiv.org/abs/1712.05834))

<sub>[\[TOC\]](#table-of-contents)</sub>
