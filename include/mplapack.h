/*
 * Copyright (c) 2008-2026
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

#ifndef _MPLAPACK_H_
#define _MPLAPACK_H_

#if defined ___MPLAPACK_BUILD_WITH_GMP___
#include <mplapack_gmp.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mpf_class REAL;
typedef mpfc_class COMPLEX;
#define Mlsamen Mlsamen_gmp
#define Mlsame Mlsame_gmp
#define Mxerbla Mxerbla_gmp
#define Rlamch Rlamch_gmp
#define Rroundup_lwork Rroundup_lwork_gmp
#define iMlaenv2stage iMlaenv2stage_gmp
#define iMlaenv iMlaenv_gmp
#define iMlaver iMlaver_gmp
#define iMieeeck iMieeeck_gmp
#define iMparmq iMparmq_gmp
#define iMparam2stage iMparam2stage_gmp
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
#include <mplapack_mpfr.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mpfr_class REAL;
typedef mpc_class COMPLEX;
#define Mlsamen Mlsamen_mpfr
#define Mlsame Mlsame_mpfr
#define Mxerbla Mxerbla_mpfr
#define Rlamch Rlamch_mpfr
#define Rroundup_lwork Rroundup_lwork_mpfr
#define iMlaenv2stage iMlaenv2stage_mpfr
#define iMlaenv iMlaenv_mpfr
#define iMlaver iMlaver_mpfr
#define iMieeeck iMieeeck_mpfr
#define iMparmq iMparmq_mpfr
#define iMparam2stage iMparam2stage_mpfr
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
#include <mplapack_qd.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef qd_real REAL;
typedef qd_complex COMPLEX;
#define Mlsamen Mlsamen_qd
#define Mlsame Mlsame_qd
#define Mxerbla Mxerbla_qd
#define Rlamch Rlamch_qd
#define Rroundup_lwork Rroundup_lwork_qd
#define iMlaenv2stage iMlaenv2stage_qd
#define iMlaenv iMlaenv_qd
#define iMlaver iMlaver_qd
#define iMieeeck iMieeeck_qd
#define iMparmq iMparmq_qd
#define iMparam2stage iMparam2stage_qd
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
#include <mplapack_dd.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef dd_real REAL;
typedef dd_complex COMPLEX;
#define Mlsamen Mlsamen_dd
#define Mlsame Mlsame_dd
#define Mxerbla Mxerbla_dd
#define Rlamch Rlamch_dd
#define Rroundup_lwork Rroundup_lwork_dd
#define iMlaenv2stage iMlaenv2stage_dd
#define iMlaenv iMlaenv_dd
#define iMlaver iMlaver_dd
#define iMieeeck iMieeeck_dd
#define iMparmq iMparmq_dd
#define iMparam2stage iMparam2stage_dd
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#include <mplapack_double.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef double REAL;
typedef std::complex<double> COMPLEX;
#define Mlsame Mlsame_double
#define Mlsamen Mlsamen_double
#define Mxerbla Mxerbla_double
#define Rlamch Rlamch_double
#define Rroundup_lwork Rroundup_lwork_double
#define iMlaenv iMlaenv_double
#define iMlaenv2stage iMlaenv2stage_double
#define iMlaver iMlaver_double
#define iMieeeck iMieeeck_double
#define iMparmq iMparmq_double
#define iMparam2stage iMparam2stage_double
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
#include <mplapack_binary80.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mplapack_binary80_t REAL;
typedef std::complex<mplapack_binary80_t> COMPLEX;
#define Mlsame Mlsame_binary80
#define Mlsamen Mlsamen_binary80
#define Mxerbla Mxerbla_binary80
#define Rlamch Rlamch_binary80
#define Rroundup_lwork Rroundup_lwork_binary80
#define iMlaver iMlaver_binary80
#define iMlaenv iMlaenv_binary80
#define iMlaenv2stage iMlaenv2stage_binary80
#define iMieeeck iMieeeck_binary80
#define iMparmq iMparmq_binary80
#define iMparam2stage iMparam2stage_binary80
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#include <mplapack_binary128.h>
typedef mplapackint INTEGER;
typedef mplapacklogical LOGICAL;
typedef mplapack_binary128_t REAL;
typedef std::complex<mplapack_binary128_t> COMPLEX;
#define Mlsame Mlsame_binary128
#define Mlsamen Mlsamen_binary128
#define Mxerbla Mxerbla_binary128
#define Rlamch Rlamch_binary128
#define Rroundup_lwork Rroundup_lwork_binary128
#define iMlaver iMlaver_binary128
#define iMlaenv iMlaenv_binary128
#define iMlaenv2stage iMlaenv2stage_binary128
#define iMieeeck iMieeeck_binary128
#define iMparmq iMparmq_binary128
#define iMparam2stage iMparam2stage_binary128
#endif

#endif
