/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mblas.h,v 1.17 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MBLAS_H_
#define _MBLAS_H_

#define ___MPACK_DEFAULT_PRECISION___ 512

#if defined ___MPACK_BUILD_WITH_MPFR___
#include <mblas_mpfr.h>
typedef mpackint INTEGER;
typedef mpreal REAL;
typedef mpcomplex COMPLEX;
#define Mlsame   Mlsame_mpfr
#define Mxerbla  Mxerbla_mpfr
#endif

#if defined ___MPACK_BUILD_WITH_GMP___
#include <mblas_gmp.h>
typedef mpackint INTEGER;
typedef mpc_class COMPLEX;
typedef mpf_class REAL;
#define Mlsame   Mlsame_gmp
#define Mxerbla  Mxerbla_gmp 
#endif

#if defined ___MPACK_BUILD_WITH_QD___
#include <mblas_qd.h>
typedef mpackint INTEGER;
typedef qd_real REAL;
typedef qd_complex COMPLEX;
#define Mlsame   Mlsame_qd
#define Mxerbla  Mxerbla_qd
#endif

#if defined ___MPACK_BUILD_WITH_DD___
#include <mblas_dd.h>
typedef mpackint INTEGER;
typedef dd_real REAL;
typedef dd_complex COMPLEX;
#define Mlsame   Mlsame_dd
#define Mxerbla  Mxerbla_dd
#endif

#if defined ___MPACK_BUILD_WITH_DOUBLE___
#include <mblas_double.h>
typedef mpackint INTEGER;
typedef double REAL;
typedef std::complex<double> COMPLEX;
#define Mlsame   Mlsame_double
#define Mxerbla  Mxerbla_double
#endif

#if defined ___MPACK_BUILD_WITH___FLOAT128___
#include <mblas___float128.h>
typedef mpackint INTEGER;
typedef __float128 REAL;
typedef std::complex<__float128> COMPLEX;
#define Mlsame   Mlsame___float128
#define Mxerbla  Mxerbla___float128
#endif

#include <cstdlib>
using std::max;
using std::min;
using std::abs;
using std::sqrt;

#endif
