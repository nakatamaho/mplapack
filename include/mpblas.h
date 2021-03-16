/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas.h,v 1.17 2010/08/07 03:15:46 nakatamaho Exp $
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#ifndef _MPBLAS_H_
#define _MPBLAS_H_

#define ___MPLAPACK_DEFAULT_PRECISION___ 512

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
#include <mpblas_mpfr.h>
typedef mplapackint INTEGER;
typedef mpreal REAL;
typedef mpcomplex COMPLEX;
#define Mlsame   Mlsame_mpfr
#define Mxerbla  Mxerbla_mpfr
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
#include <mpblas_gmp.h>
typedef mplapackint INTEGER;
typedef mpc_class COMPLEX;
typedef mpf_class REAL;
#define Mlsame   Mlsame_gmp
#define Mxerbla  Mxerbla_gmp 
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
#include <mpblas_qd.h>
typedef mplapackint INTEGER;
typedef qd_real REAL;
typedef qd_complex COMPLEX;
#define Mlsame   Mlsame_qd
#define Mxerbla  Mxerbla_qd
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
#include <mpblas_dd.h>
typedef mplapackint INTEGER;
typedef dd_real REAL;
typedef dd_complex COMPLEX;
#define Mlsame   Mlsame_dd
#define Mxerbla  Mxerbla_dd
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#include <mpblas_double.h>
typedef mplapackint INTEGER;
typedef double REAL;
typedef std::complex<double> COMPLEX;
#define Mlsame   Mlsame_double
#define Mxerbla  Mxerbla_double
#endif

#if defined ___MPLAPACK_BUILD_WITH_LONGDOUBLE___
#include <mpblas_longdouble.h>
typedef mplapackint INTEGER;
typedef long double REAL;
typedef std::complex<long double> COMPLEX;
#define Mlsame   Mlsame_longdouble
#define Mxerbla  Mxerbla_longdouble
#endif


#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
#include <mpblas__Float128.h>
typedef mplapackint INTEGER;
typedef _Float128 REAL;
typedef std::complex<_Float128> COMPLEX;
#define Mlsame   Mlsame__Float128
#define Mxerbla  Mxerbla__Float128
#endif

#include <algorithm>    // std::max
#include <math.h> // sqrt
#include <cstdlib>
using std::max;
using std::min;
using std::abs;
using std::sqrt;

#endif
