/*
 * Copyright (c) 2008-2022
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

#include <sys/types.h>
#if !defined _WIN32
#include <sys/resource.h>
#endif

#define KFLOPS 1e-3
#define MFLOPS 1e-6
#define GFLOPS 1e-9
#define TFLOPS 1e-12

#if defined __APPLE__
#define DYLIB_SUFFIX ".dylib"
#elif defined _WIN32
#define DYLIB_SUFFIX "-0.dll" //XXX
#else
#define DYLIB_SUFFIX ".so"
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
#include <mplapack_benchmark_mpfr.h>
#define MPLAPACK_REF_LIB "libmplapack_mpfr"
#define MPBLAS_REF_LIB "libmpblas_mpfr"
#define ___MPLAPACK_INITIALIZE___ gmp_randinit_default(uniformrandomstate_mpfr);
gmp_randstate_t uniformrandomstate_mpfr;
mpreal randomnumber(mpreal dummy) {
    mpreal mtmp;

    mtmp = urandomb(uniformrandomstate_mpfr);
    mtmp = 2.0 * mtmp - 1.0;

    return mtmp;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
#include <mplapack_benchmark_gmp.h>
#define MPLAPACK_REF_LIB "libmplapack_gmp"
#define MPBLAS_REF_LIB "libmpblas_gmp"
#define ___MPLAPACK_INITIALIZE___ uniformrandomstate_gmp = new gmp_randclass(gmp_randinit_default);
gmp_randclass *uniformrandomstate_gmp;
mpf_class randomnumber(mpf_class dummy) {
    mpf_class mtmp;

    mtmp = uniformrandomstate_gmp->get_f();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
#include <mplapack_benchmark_dd.h>
#define MPLAPACK_REF_LIB "libmplapack_dd"
#define MPBLAS_REF_LIB "libmpblas_dd"
dd_real randomnumber(dd_real dummy) {
    dd_real mtmp;
    mtmp = ddrand();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#define ___MPLAPACK_INITIALIZE___
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
#include <mplapack_benchmark_qd.h>
#define MPLAPACK_REF_LIB "libmplapack_qd"
#define MPBLAS_REF_LIB "libmpblas_qd"
#define ___MPLAPACK_INITIALIZE___
qd_real randomnumber(qd_real dummy) {
    qd_real mtmp;
    mtmp = qdrand(); // uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#include <mplapack_benchmark_double.h>
#define MPLAPACK_REF_LIB "libmplapack_double"
#define MPBLAS_REF_LIB "libmpblas_double"
#define ___MPLAPACK_INITIALIZE___
double randomnumber(double dummy) {
#if defined _WIN32 // XXX
    double mtmp = (double)rand() / (double)RAND_MAX;
#else
    double mtmp = drand48();
#endif
    return mtmp;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
#include <mplapack_benchmark__Float64x.h>
#define MPLAPACK_REF_LIB "libmplapack__Float64x"
#define MPBLAS_REF_LIB "libmpblas__Float64x"
#define ___MPLAPACK_INITIALIZE___
_Float64x randomnumber(_Float64x dummy) {
    _Float64x mtmp;
#if defined _WIN32
    mtmp = ((double)rand() / (double)RAND_MAX);          // uniform random between [0,1] via rand
    mtmp += ((double)rand() / (double)RAND_MAX) * 1e-16; // uniform random between [0,1] via rand
#else
    mtmp = drand48();          // uniform random between [0,1] via lrand48
    mtmp += drand48() * 1e-16; // uniform random between [0,1] via lrand48
#endif
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
#include <mplapack_benchmark__Float128.h>
#define MPLAPACK_REF_LIB "libmplapack__Float128"
#define MPBLAS_REF_LIB "libmpblas__Float128"
#define ___MPLAPACK_INITIALIZE___
_Float128 randomnumber(_Float128 dummy) {
    _Float128 mtmp;
#if defined _WIN32
    mtmp = ((double)rand() / (double)RAND_MAX);          // uniform random between [0,1] via rand
    mtmp += ((double)rand() / (double)RAND_MAX) * 1e-16; // uniform random between [0,1] via rand
    mtmp += ((double)rand() / (double)RAND_MAX) * 1e-32; // uniform random between [0,1] via rand
#else
    mtmp = drand48();          // uniform random between [0,1] via drand48
    mtmp += drand48() * 1e-16; // uniform random between [0,1] via drand48
    mtmp += drand48() * 1e-32; // uniform random between [0,1] via drand48
#endif
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined ___DOUBLE_BENCH___
double randomnumber(double dummy) {
#if defined _WIN32 // XXX
    double mtmp = ((double)rand() / (double)RAND_MAX);
#else
    double mtmp = drand48();
#endif
    return mtmp;
}
#endif
