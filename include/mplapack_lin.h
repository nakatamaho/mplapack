/*
 * Copyright (c) 2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack.h,v 1.28 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_LIN_H_
#define _MPLAPACK_LIN_H_

#if defined ___MPLAPACK_BUILD_WITH_GMP___
#include <mplapack_lin_gmp.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mpf_class REAL;
typedef mpc_class COMPLEX;
#define Mlsame Mlsame_gmp
#define Mxerbla Mxerbla_gmp
#define Rlamch Rlamch_gmp
#define iMlaenv iMlaenv_gmp
#define iMieeeck iMieeeck_gmp
#define iMparmq iMparmq_gmp
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
#include <mplapack_lin_mpfr.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mpreal REAL;
typedef mpcomplex COMPLEX;
#define Mlsame Mlsame_mpfr
#define Mxerbla Mxerbla_mpfr
#define Rlamch Rlamch_mpfr
#define iMlaenv iMlaenv_mpfr
#define iMieeeck iMieeeck_mpfr
#define iMparmq iMparmq_mpfr
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
#include <mplapack_lin_qd.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef qd_real REAL;
typedef qd_complex COMPLEX;
#define Mlsame Mlsame_qd
#define Mxerbla Mxerbla_qd
#define Rlamch Rlamch_qd
#define iMlaenv iMlaenv_qd
#define iMieeeck iMieeeck_qd
#define iMparmq iMparmq_qd
#if !defined __MUTILS_CPP__
#define nint __qd_nint
#endif
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
#include <mplapack_lin_dd.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef dd_real REAL;
typedef dd_complex COMPLEX;
#define Mlsame Mlsame_dd
#define Mxerbla Mxerbla_dd
#define Rlamch Rlamch_dd
#define iMlaenv iMlaenv_dd
#define iMieeeck iMieeeck_dd
#define iMparmq iMparmq_dd
#if !defined __MUTILS_CPP__
#define nint __dd_nint
#endif
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#include <mplapack_lin_double.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef double REAL;
typedef std::complex<double> COMPLEX;
#define Mlsame Mlsame_double
#define Mxerbla Mxerbla_double
#define Rlamch Rlamch_double
#define iMlaenv iMlaenv_double
#define iMieeeck iMieeeck_double
#define iMparmq iMparmq_double
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
#include <mplapack_lin__Float64x.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef _Float64x REAL;
typedef std::complex<_Float64x> COMPLEX;
#define Mlsame Mlsame__Float64x
#define Mxerbla Mxerbla__Float64x
#define Rlamch Rlamch__Float64x
#define iMlaenv iMlaenv__Float64x
#define iMieeeck iMieeeck__Float64x
#define iMparmq iMparmq__Float64x
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
#include <mplapack_lin__Float128.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef _Float128 REAL;
typedef std::complex<_Float128> COMPLEX;
#define Mlsame Mlsame__Float128
#define Mxerbla Mxerbla__Float128
#define Rlamch Rlamch__Float128
#define iMlaenv iMlaenv__Float128
#define iMieeeck iMieeeck__Float128
#define iMparmq iMparmq__Float128
#endif

#endif
