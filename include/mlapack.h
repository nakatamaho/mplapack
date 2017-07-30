/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mlapack.h,v 1.28 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MLAPACK_H_
#define _MLAPACK_H_

#if defined ___MPACK_BUILD_WITH_GMP___
#include <mlapack_gmp.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef mpf_class REAL;
typedef mpc_class COMPLEX;
#define Mlsame   Mlsame_gmp
#define Mxerbla  Mxerbla_gmp
#define Rlamch   Rlamch_gmp
#define iMlaenv  iMlaenv_gmp
#endif

#if defined ___MPACK_BUILD_WITH_MPFR___
#include <mlapack_mpfr.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef mpreal REAL;
typedef mpcomplex COMPLEX;
#define Mlsame   Mlsame_mpfr
#define Mxerbla  Mxerbla_mpfr
#define Rlamch   Rlamch_mpfr
#define iMlaenv  iMlaenv_mpfr
#endif

#if defined ___MPACK_BUILD_WITH_QD___
#include <mlapack_qd.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef qd_real REAL;
typedef qd_complex COMPLEX;
#if !defined __MUTILS_CPP__
#define nint __qd_nint
#endif
#define Mlsame   Mlsame_qd
#define Mxerbla  Mxerbla_qd
#define Rlamch   Rlamch_qd
#define iMlaenv  iMlaenv_qd
#endif

#if defined ___MPACK_BUILD_WITH_DD___
#include <mlapack_dd.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef dd_real REAL;
typedef dd_complex COMPLEX;
#if !defined __MUTILS_CPP__
#define nint __dd_nint
#endif
#define Mlsame   Mlsame_dd
#define Mxerbla  Mxerbla_dd
#define Rlamch   Rlamch_dd
#define iMlaenv  iMlaenv_dd
#endif

#if defined ___MPACK_BUILD_WITH_DOUBLE___
#include <mlapack_double.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef double REAL;
typedef std::complex<double> COMPLEX;
#define Mlsame   Mlsame_double
#define Mxerbla  Mxerbla_double
#define Rlamch   Rlamch_double
#define iMlaenv  iMlaenv_double
#endif

#if defined ___MPACK_BUILD_WITH_LONGDOUBLE___
#include <mlapack_longdouble.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef long double REAL;
typedef std::complex<long double> COMPLEX;
#define Mlsame   Mlsame_longdouble
#define Mxerbla  Mxerbla_longdouble
#define Rlamch   Rlamch_longdouble
#define iMlaenv  iMlaenv_longdouble
#endif


#if defined ___MPACK_BUILD_WITH___FLOAT128___
#include <mlapack___float128.h>
typedef mpackint INTEGER;
typedef mpacklogical LOGICAL;
typedef __float128 REAL;
typedef std::complex<__float128> COMPLEX;
#define Mlsame   Mlsame___float128
#define Mxerbla  Mxerbla___float128
#define Rlamch   Rlamch___float128
#define iMlaenv  iMlaenv___float128
#endif

#endif

