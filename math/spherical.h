/*******************************************************************************
* spherical.h: this file is part of the powspec program.

* powspec: C code for auto and cross power spectra evaluation.

* Github repository:
        https://github.com/cheng-zhao/powspec

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef _SPHERICAL_H_
#define _SPHERICAL_H_

/******************************************************************************
  Implementation of the real-form spherical harmonics:
            /  (-1)^m / sqrt(2) * (Ylm + Ylm^*)          ,  m > 0
    YlmR = <   Yl0                                       ,  m = 0
            \  (-1)^m / (i * sqrt(2)) (Yl|m| - Yl|m|^*)  ,  m < 0
******************************************************************************/

/******************************************************************************
Functions `YlmR_l<l>`:
  Compute the real-form spherical harmonics Y(l, m, theta, phi).
Arguments:
  * `m`:        m
  * `cost`:     cos(theta)
  * `sint`:     sin(theta)
  * `cosp`:     cos(phi)
  * `sinp`:     sin(phi)
Return:
  Y(l, m, theta, phi).
******************************************************************************/

inline double YlmR_l1(const int m, const double cost, const double sint,
    const double cosp, const double sinp) {
  switch (m) {
    case -1: return 0x1.f45437857749ap-2 * sint * sinp;
    case 0: return 0x1.f45437857749ap-2 * cost;
    case 1: return 0x1.f45437857749ap-2 * sint * cosp;
    default: return 1;
  }
}

inline double YlmR_l2(const int m, const double cost, const double sint,
    const double cosp, const double sinp) {
  switch (m) {
    case -2: return 0x1.17b14102b0694p+0 * sint * sint * cosp * sinp;
    case -1: return 0x1.17b14102b0694p+0 * cost * sint * sinp;
    case 0: return 0x1.e471027d29b67p-1 * (cost * cost - 0x1.5555555555555p-2);
    case 1: return 0x1.17b14102b0694p+0 * cost * sint * cosp;
    case 2: return 0x1.17b14102b0694p+0 * sint * sint * (cosp * cosp - 0.5);
    default: return 1;
  }
}

inline double YlmR_l3(const int m, const double cost, const double sint,
    const double cosp, const double sinp) {
  switch (m) {
    case -3:
      return 0x1.2e1a3183e613bp+1 * sint * sint * sint
        * sinp * (0.75 - sinp * sinp);
    case -2:
      return 0x1.71ff8e45cad2ep+1 * cost * sint * sint * cosp * sinp;
    case -1:
      return 0x1.2482623faf652p+1 * (cost * cost - 0.2) * sint * sinp;
    case 0:
      return 0x1.ddaa6bb0942dcp+0 * cost * (cost * cost - 0.6);
    case 1:
      return 0x1.2482623faf652p+1 * (cost * cost - 0.2) * sint * cosp;
    case 2:
      return 0x1.71ff8e45cad2ep+1 * cost * sint * sint * (cosp * cosp - 0.5);
    case 3:
      return 0x1.2e1a3183e613bp+1 * sint * sint * sint * cosp
        * (cosp * cosp - 0.75);
    default: return 1;
  }
}

inline double YlmR_l4(const int m, const double cost, const double sint,
    const double cosp, const double sinp) {
  double tmp, tmp2;
  switch (m) {
    case -4:
      tmp = sint * sint;
      return 0x1.406d8aa0d83a4p+2 * tmp * tmp * cosp * (cosp * cosp - 0.5)
        * sinp;
    case -3:
      return 0x1.c5274a45d91d8p+2 * cost * sint * sint * sint * sinp *
        (0.75 - sinp * sinp);
    case -2:
      tmp = cost * cost;
      return 0x1.e471027d29b67p+2 * (tmp - 0.875 * tmp * tmp - 0.125)
        * cosp * sinp;
    case -1:
      return 0x1.2bbb9c2824048p+2 * (cost * cost - 0x1.b6db6db6db6dbp-2)
        * cost * sint * sinp;
    case 0:
      tmp = cost * cost;
      return 0x1.96376c8ddaa91p+1 * (0x1.2aaaaaaaaaaabp+0 * tmp * tmp
          - tmp + 0.1);
    case 1:
      return 0x1.2bbb9c2824048p+2 * (cost * cost - 0x1.b6db6db6db6dbp-2)
        * cost * sint * cosp;
    case 2:
      tmp = cost * cost;
      return 0x1.e471027d29b67p+2 * (tmp - 0.875 * tmp * tmp - 0.125)
        * (cosp * cosp - 0.5);
    case 3:
      return 0x1.c5274a45d91d8p+2 * cost * sint * sint * sint * cosp *
        (cosp * cosp - 0.75);
    case 4:
      tmp = sint * sint;
      tmp2 = cosp * cosp;
      return 0x1.406d8aa0d83a4p+2 * tmp * tmp * (tmp2 * tmp2 - tmp2 + 0.125);
    default: return 1;
  }
}

inline double YlmR_l5(const int m, const double cost, const double sint,
    const double cosp, const double sinp) {
  double t1, t2;
  switch (m) {
    case -5:
      t1 = sint * sint;
      t2 = cosp * cosp;
      return 0x1.50114f179e96cp+3 * t1 * t1 * sint * sinp
        * (t2 * t2 - 0.75 * t2 + 0.0625);
    case -4:
      t1 = sint * sint;
      return 0x1.09af4d7ffa13bp+4 * cost * t1 * t1 * cosp * (cosp * cosp - 0.5)
        * sinp;
    case -3:
      t1 = cost * cost;
      return 0x1.391ccd9010092p+4 * (t1 - 0.9 * t1 * t1 - 0.1) * sint
        * sinp * (0.75 - sinp * sinp);
    case -2:
      t1 = cost * cost;
      return 0x1.32c94e82e889dp+4 * cost * (t1 - 0.75 * t1 * t1 - 0.25)
        * cosp * sinp;
    case -1:
      t1 = cost * cost;
      return 0x1.95d71750038a7p+2 * (1.5 * t1 * t1 - t1 + 0x1.2492492492492p-4)
        * sint * sinp;
    case 0:
      t1 = cost * cost;
      return 0x1.05f7fe2f32cdcp+3 * cost
        * (0.9 * t1 * t1 - t1 + 0x1.b6db6db6db6dbp-3);
    case 1:
      t1 = cost * cost;
      return 0x1.95d71750038a7p+2 * (1.5 * t1 * t1 - t1 + 0x1.2492492492492p-4)
        * sint * cosp;
    case 2:
      t1 = cost * cost;
      return 0x1.32c94e82e889dp+4 * cost * (t1 - 0.75 * t1 * t1 - 0.25)
        * (cosp * cosp - 0.5);
    case 3:
      t1 = cost * cost;
      return 0x1.391ccd9010092p+4 * (t1 - 0.9 * t1 * t1 - 0.1) * sint
        * cosp * (cosp * cosp - 0.75);
    case 4:
      t1 = sint * sint;
      t2 = cosp * cosp;
      return 0x1.09af4d7ffa13bp+4 * cost * t1 * t1 * (t2 * t2 - t2 + 0.125);
    case 5:
      t1 = sint * sint;
      t2 = cosp * cosp;
      return 0x1.50114f179e96cp+3 * t1 * t1 * sint * cosp
        * (t2 * t2 - 1.25 * t2 + 0.3125);
    default: return 1;
  }
}

inline double YlmR_l6(const int m, const double cost, const double sint,
    const double cosp, const double sinp) {
  double t1, t2;
  switch (m) {
    case -6:
      t1 = sint * sint * sint;
      t2 = cosp * cosp;
      return 0x1.5dca4e99e480ep+4 * t1 * t1 * cosp
        * (t2 * t2 - t2 + 0.1875) * sinp;
    case -5:
      t1 = sint * sint;
      t2 = cosp * cosp;
      return 0x1.2eed606fefa75p+5 * cost * t1 * t1 * sint * sinp
        * (t2 * t2 - 0.75 * t2 + 0.0625);
    case -4:
      t1 = cost * cost;
      return 0x1.6336b4652ee84p+5 * (t1 * t1 * t1 - 0x1.0ba2e8ba2e8bap+1
          * t1 * t1 + 0x1.2e8ba2e8ba2e9p+0 * t1 - 0x1.745d1745d1746p-4)
        * cosp * sinp * (cosp * cosp - 0.5);
    case -3:
      t1 = cost * cost;
      return 0x1.9cb3305568dd9p+5 * cost * (t1 - 0x1.9249249249249p-1 * t1 * t1
          - 0x1.b6db6db6db6dbp-3) * sint * sinp * (0.75 - sinp * sinp);
    case -2:
      t1 = cost * cost;
      return 0x1.e66578f6f272ep+4 * (0x1.8ba2e8ba2e8bap+0 * t1 * t1
          - t1 * t1 * t1 - 0x1.26c9b26c9b26dp-1 * t1 + 0x1.f07c1f07c1f08p-6)
        * cosp * sinp;
    case -1:
      t1 = cost * cost;
      return 0x1.339fc3ab0f95dp+4 * cost * (t1 * t1 - 0x1.d1745d1745d17p-1 * t1
          + 0x1.364d9364d9365p-3) * sint * sinp;
    case 0:
      t1 = cost * cost;
      return 0x1.d5e74e9acca4ap+3 * (t1 * (t1 * t1 - 0x1.5d1745d1745d1p+0 * t1
            + 0x1.d1745d1745d17p-2) - 0x1.62a1cd058a873p-6);
    case 1:
      t1 = cost * cost;
      return 0x1.339fc3ab0f95dp+4 * cost * (t1 * t1 - 0x1.d1745d1745d17p-1 * t1
          + 0x1.364d9364d9365p-3) * sint * cosp;
    case 2:
      t1 = cost * cost;
      return 0x1.e66578f6f272ep+4 * (0x1.8ba2e8ba2e8bap+0 * t1 * t1
          - t1 * t1 * t1 - 0x1.26c9b26c9b26dp-1 * t1 + 0x1.f07c1f07c1f08p-6)
        * (cosp * cosp - 0.5);
    case 3:
      t1 = cost * cost;
      return 0x1.9cb3305568dd9p+5 * cost * (t1 - 0x1.9249249249249p-1 * t1 * t1
          - 0x1.b6db6db6db6dbp-3) * sint * cosp * (cosp * cosp - 0.75);
    case 4:
      t1 = cost * cost;
      t2 = cosp * cosp;
      return 0x1.6336b4652ee84p+5 * (t1 * t1 * t1 - 0x1.0ba2e8ba2e8bap+1
          * t1 * t1 + 0x1.2e8ba2e8ba2e9p+0 * t1 - 0x1.745d1745d1746p-4)
        * (t2 * t2 - t2 + 0.125);
    case 5:
      t1 = sint * sint;
      t2 = cosp * cosp;
      return 0x1.2eed606fefa75p+5 * cost * t1 * t1 * sint * cosp
        * (t2 * t2 - 1.25 * t2 + 0.3125);
    case 6:
      t1 = sint * sint * sint;
      t2 = cosp * cosp;
      return 0x1.5dca4e99e480ep+4 * t1 * t1
        * (t2 * t2 * t2 - 1.5 * t2 * t2 + 0.5625 * t2 - 0.03125);
    default: return 1;
  }
}

#endif
