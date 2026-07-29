/*
 * Copyright (c) 2021
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

#ifndef MPLAPACK_EIG_H
#define MPLAPACK_EIG_H

#if defined MPLAPACK_BUILD_WITH_GMP
#include <mplapack_eig_gmp.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mpf_class REAL;
typedef mpfc_class COMPLEX;
#define Mlsame Mlsame_gmp
#define Mxerbla Mxerbla_gmp
#define Rlamch Rlamch_gmp
#define iMlaenv iMlaenv_gmp
#define iMieeeck iMieeeck_gmp
#define iMparmq iMparmq_gmp
#endif

#if defined MPLAPACK_BUILD_WITH_MPFR
#include <mplapack_eig_mpfr.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mpfr_class REAL;
typedef mpc_class COMPLEX;
#define Mlsame Mlsame_mpfr
#define Mxerbla Mxerbla_mpfr
#define Rlamch Rlamch_mpfr
#define iMlaenv iMlaenv_mpfr
#define iMieeeck iMieeeck_mpfr
#define iMparmq iMparmq_mpfr
#endif

#if defined MPLAPACK_BUILD_WITH_QD
#include "mplapack_eig_qd.h"
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
#endif

#if defined MPLAPACK_BUILD_WITH_DD
#include <mplapack_eig_dd.h>
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
#endif

#if defined MPLAPACK_BUILD_WITH_DOUBLE
#include <mplapack_eig_double.h>
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

#if defined MPLAPACK_BUILD_WITH_BINARY80
#include <mplapack_eig_binary80.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mplapack_binary80_t REAL;
typedef std::complex<mplapack_binary80_t> COMPLEX;
#define Mlsame Mlsame_binary80
#define Mxerbla Mxerbla_binary80
#define Rlamch Rlamch_binary80
#define iMlaenv iMlaenv_binary80
#define iMieeeck iMieeeck_binary80
#define iMparmq iMparmq_binary80
#endif

#if defined MPLAPACK_BUILD_WITH_BINARY128
#include <mplapack_eig_binary128.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mplapack_binary128_t REAL;
typedef std::complex<mplapack_binary128_t> COMPLEX;
#define Mlsame Mlsame_binary128
#define Mxerbla Mxerbla_binary128
#define Rlamch Rlamch_binary128
#define iMlaenv iMlaenv_binary128
#define iMieeeck iMieeeck_binary128
#define iMparmq iMparmq_binary128
#endif

#endif
