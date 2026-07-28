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
#define DYLIB_SUFFIX "-0.dll" // XXX
#else
#define DYLIB_SUFFIX ".so"
#endif

using std::max;
using std::min;

#if defined MPLAPACK_BUILD_WITH_MPFR
#include <mplapack_benchmark_mpfr.h>
#define MPLAPACK_REF_LIB "libmplapack_mpfr"
#define MPBLAS_REF_LIB "libmplapack_mpfr"
#define MPLAPACK_INITIALIZE uniformrandomstate_mpfr.seed(0UL);
mpfrxx::mpfr_randclass uniformrandomstate_mpfr(gmp_randinit_default);
mpfr_class randomnumber(mpfr_class dummy) {
    mpfr_class mtmp;

    mtmp = uniformrandomstate_mpfr.get_fr();
    mtmp = 2.0 * mtmp - 1.0;

    return mtmp;
}
#endif

#if defined MPLAPACK_BUILD_WITH_GMP
#include <mplapack_benchmark_gmp.h>
#define MPLAPACK_REF_LIB "libmplapack_gmp"
#define MPBLAS_REF_LIB "libmplapack_gmp"
#define MPLAPACK_INITIALIZE uniformrandomstate_gmp.seed(0UL);
gmpxx::gmp_randclass uniformrandomstate_gmp(gmp_randinit_default);
mpf_class randomnumber(mpf_class dummy) {
    mpf_class mtmp;

    mtmp = uniformrandomstate_gmp.get_f();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined MPLAPACK_BUILD_WITH_DD
#include <mplapack_benchmark_dd.h>
#define MPLAPACK_REF_LIB "libmplapack_dd"
#define MPBLAS_REF_LIB "libmplapack_dd"
dd_real randomnumber(dd_real dummy) {
    dd_real mtmp;
    mtmp = ddrand();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#define MPLAPACK_INITIALIZE
#endif

#if defined MPLAPACK_BUILD_WITH_QD
#include <mplapack_benchmark_qd.h>
#define MPLAPACK_REF_LIB "libmplapack_qd"
#define MPBLAS_REF_LIB "libmplapack_qd"
#define MPLAPACK_INITIALIZE
qd_real randomnumber(qd_real dummy) {
    qd_real mtmp;
    mtmp = qdrand(); // uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}
#endif

#if defined MPLAPACK_BUILD_WITH_DOUBLE
#include <mplapack_benchmark_double.h>
#define MPLAPACK_REF_LIB "libmplapack_double"
#define MPBLAS_REF_LIB "libmplapack_double"
#define MPLAPACK_INITIALIZE
double randomnumber(double dummy) {
#if defined _WIN32 // XXX
    double mtmp = (double)rand() / (double)RAND_MAX;
#else
    double mtmp = drand48();
#endif
    return mtmp;
}
#endif

#if defined MPLAPACK_BUILD_WITH_BINARY80
#include <mplapack_benchmark_binary80.h>
#define MPLAPACK_REF_LIB "libmplapack_binary80"
#define MPBLAS_REF_LIB "libmplapack_binary80"
#define MPLAPACK_INITIALIZE
mplapack_binary80_t randomnumber(mplapack_binary80_t dummy) {
    mplapack_binary80_t mtmp;
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

#if defined MPLAPACK_BUILD_WITH_BINARY128
#include <mplapack_benchmark_binary128.h>
#define MPLAPACK_REF_LIB "libmplapack_binary128"
#define MPBLAS_REF_LIB "libmplapack_binary128"
#define MPLAPACK_INITIALIZE
mplapack_binary128_t randomnumber(mplapack_binary128_t dummy) {
    mplapack_binary128_t mtmp;
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
