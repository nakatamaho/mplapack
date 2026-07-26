/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
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

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
#include <mpblas_mpfr.h>
typedef mplapackint INTEGER;
typedef mpfr_class REAL;
typedef mpc_class COMPLEX;
#define Mlsame Mlsame_mpfr
#define Mxerbla Mxerbla_mpfr
#define castINTEGER castINTEGER_mpfr
#define castREAL castREAL_mpfr
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
#include <mpblas_gmp.h>
typedef mplapackint INTEGER;
typedef mpfc_class COMPLEX;
typedef mpf_class REAL;
#define Mlsame Mlsame_gmp
#define Mxerbla Mxerbla_gmp
#define castINTEGER castINTEGER_gmp
#define castREAL castREAL_gmp
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
#include <mpblas_qd.h>
typedef mplapackint INTEGER;
typedef qd_real REAL;
typedef qd_complex COMPLEX;
#define Mlsame Mlsame_qd
#define Mxerbla Mxerbla_qd
#define castINTEGER castINTEGER_qd
#define castREAL castREAL_qd
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
#include <mpblas_dd.h>
typedef mplapackint INTEGER;
typedef dd_real REAL;
typedef dd_complex COMPLEX;
#define Mlsame Mlsame_dd
#define Mxerbla Mxerbla_dd
#define castINTEGER castINTEGER_dd
#define castREAL castREAL_dd
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#include <mpblas_double.h>
typedef mplapackint INTEGER;
typedef double REAL;
typedef std::complex<double> COMPLEX;
#define Mlsame Mlsame_double
#define Mxerbla Mxerbla_double
#define castINTEGER castINTEGER_double
#define castREAL castREAL_double
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
#include <mpblas_binary80.h>
typedef mplapackint INTEGER;
typedef mplapack_binary80_t REAL;
typedef std::complex<mplapack_binary80_t> COMPLEX;
#define Mlsame Mlsame_binary80
#define Mxerbla Mxerbla_binary80
#define castINTEGER castINTEGER_binary80
#define castREAL castREAL_binary80
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#include <mpblas_binary128.h>
typedef mplapackint INTEGER;
typedef mplapack_binary128_t REAL;
typedef std::complex<mplapack_binary128_t> COMPLEX;
#define Mlsame Mlsame_binary128
#define Mxerbla Mxerbla_binary128
#define castINTEGER castINTEGER_binary128
#define castREAL castREAL_binary128
#endif

#endif
