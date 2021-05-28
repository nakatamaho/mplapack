/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlamch.cpp,v 1.16 2010/08/07 04:48:32 nakatamaho Exp $
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
#include <mpblas.h>
#include <mplapack.h>

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
#include <float.h>
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
#include <float.h>
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___ && defined ___MPLAPACK_WANT_LIBQUADMATH___
#include <quadmath.h>
#else
#include <float.h>
#endif

#include <stdio.h>

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
REAL RlamchE_mpfr(void) {
    static REAL eps;
    static int called = 0;

    if (called)
        return eps;

    REAL one = 1.0, rtmp;
    mp_prec_t prec;
    prec = one.get_prec();
    rtmp = exp2(prec);
    eps = one / rtmp;
    called = 1;
    return eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchS_mpfr(void) {
    static REAL safemin;
    static int called = 0;

    if (called)
        return safemin;

    REAL one = 1.0, rtmp;
    mp_exp_t emin;
    emin = one.get_emin() + 1;
    safemin = div_2ui(one, -emin);
    return safemin;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchB_mpfr(void) {
    REAL two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchP_mpfr(void) {
    REAL base, eps, prec;

    base = RlamchB_mpfr();
    eps = RlamchE_mpfr();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchN_mpfr(void) {
    REAL r;
    mp_prec_t p = r.get_default_prec();
    return REAL(p);
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchR_mpfr(void) {
    // always rounding in addition on MPFR.
    REAL mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchM_mpfr(void) {
    // Note: in MPFR, no gradual underflow, just happens suddenly.
    REAL rtmp;
    mpreal minexp(rtmp.get_emin());
    return minexp;
}

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchU_mpfr(void) {
    static REAL safemin;
    static int called = 0;

    if (called)
        return safemin;

    REAL one = 1.0, rtmp;
    mp_exp_t emin;
    emin = one.get_emin();
    safemin = div_2si(one, -emin - 1);
    return safemin;
}

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchL_mpfr(void) {
    mpreal dummy;
    mp_exp_t emax;
    emax = dummy.get_emax() - 1;
    return mpreal(emax);
}

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchO_mpfr(void) {
    ///(1-2^(-p))*2^emax
    // in double's case, +1.7976931348623157e+308 = +4.4942328371557898e+307 x 3.999999999999999(16 digits)
    // where +4.4942328371557898e+307 is the reciprocal of 1/safmin.
    // of course +4.4942328371557898e+307 times 4 overflows.

    static REAL overflow;
    static int called = 0;

    if (called)
        return overflow;

    REAL one = 1.0;
    REAL rtmp;
    mp_prec_t prec = one.get_prec();
    mp_exp_t emax = one.get_emax();
    rtmp = exp2(-(long)prec);
    overflow = mul_2si(one, emax - 1);
    overflow = overflow * (one - rtmp) * 2.0;
    return overflow;
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchZ_mpfr(void) {
    REAL mtemp = 0.0;
    return mtemp;
}

REAL Rlamch_mpfr(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE_mpfr();
    if (Mlsame(cmach, "S"))
        return RlamchS_mpfr();
    if (Mlsame(cmach, "B"))
        return RlamchB_mpfr();
    if (Mlsame(cmach, "P"))
        return RlamchP_mpfr();
    if (Mlsame(cmach, "N"))
        return RlamchN_mpfr();
    if (Mlsame(cmach, "R"))
        return RlamchR_mpfr();
    if (Mlsame(cmach, "M"))
        return RlamchM_mpfr();
    if (Mlsame(cmach, "U"))
        return RlamchU_mpfr();
    if (Mlsame(cmach, "L"))
        return RlamchL_mpfr();
    if (Mlsame(cmach, "O"))
        return RlamchO_mpfr();

    Mxerbla("Rlamch", 1);
    return RlamchZ_mpfr();
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
REAL RlamchE_gmp(void) {
    static REAL eps;
    static int called = 0;
    if (called)
        return eps;
    REAL one;
    unsigned long exp2;
    one = 1.0;
    exp2 = mpf_get_prec(one.get_mpf_t());
    mpf_div_2exp(eps.get_mpf_t(), one.get_mpf_t(), exp2);
    called = 1;
    return eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchS_gmp(void) {
    REAL sfmin;
    REAL one = 1.0;
    unsigned long exp2;

    exp2 = (1UL << (mp_bits_per_limb - 8)) - 1; // 6 seems to be the smallest on amd64 but for safty
    mpf_div_2exp(sfmin.get_mpf_t(), one.get_mpf_t(), exp2);
    return sfmin;

    /* following code fragment is to test safe minimum
        REAL largenum;
        for(int p = 60; p>=0; p--) {
          for (int a=16; a<= 5120; a=a+128) {
            sfmin = 0.0;
            mpf_set_default_prec(a);
            exp2 = (1UL << (mp_bits_per_limb-p)) - 1;
            mpf_div_2exp(sfmin.get_mpf_t(), one.get_mpf_t(), exp2);
            largenum  = 1.0 / sfmin;
            gmp_printf("%d, a:%5d, p:%5d sfmin: %16.20Fe largenum: %16.20Fe\n", mp_bits_per_limb, a, p, sfmin.get_mpf_t(), largenum.get_mpf_t());
            if (sfmin < 1.0 ) { printf("sfmin yes\n"); }
            if (largenum > 1.0 ) { printf("largenum yes\n"); }
          }
        }
    */
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchB_gmp(void) {
    REAL two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchP_gmp(void) {
    REAL base, eps, prec;

    base = RlamchB_gmp();
    eps = RlamchE_gmp();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchN_gmp(void) {
    unsigned long int tmp;
    REAL mtmp;
    REAL mtmp2;
    tmp = mpf_get_prec(mtmp.get_mpf_t());
    mtmp2 = tmp;
    return mtmp2;

    /* following is fragment of code to test digits in mantissa
       for (int a=8; a<= 5120; a=a+128) {
         mpf_set_default_prec(a);
         REAL mmtmp;
         tmp = mpf_get_prec(mmtmp.get_mpf_t());
         printf("number of digits in mantissa %d\n", (int)tmp );
       }
    */
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchR_gmp(void) {
    // always rounding in addition on GMP.
    REAL mtmp;

    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchM_gmp(void) {
    unsigned long exp2;
    REAL tmp;
    REAL uflowmin, one = 1.0;
    exp2 = (1UL << (mp_bits_per_limb - 8)) - 1; // 6 seems to be the smallest on amd64 but for safty
    tmp = exp2;
    return -tmp;

    /*
      following code fragment is to test minimum
      exponent before (gradual) underflow...but we just got Bus error
        REAL mtmp1, mtmp2;
        for(int p = 11; p>=0; p--) {
          for (int a=51200; a<= 102400; a=a+128) {
            mpf_set_default_prec(a);
            printf("p %d a %d \n", p, a);
            uflowmin=0.0;
            exp2 = (1UL << (mp_bits_per_limb-p)) - 1;
            mpf_div_2exp(uflowmin.get_mpf_t(), one.get_mpf_t(), exp2);
            mtmp1 = uflowmin/2.0;
            gmp_printf("p %d, uflowmin: %16.20Fe uflowmin/2 %16.20Fe\n", p, uflowmin.get_mpf_t(), mtmp1.get_mpf_t());
            mtmp2 = mtmp1 * 2.0;
            gmp_printf("mtmp2: %16.20Fe %lu\n", mtmp2.get_mpf_t(),exp2);
            if (uflowmin != mtmp2 ) { printf("underflow\n"); exit(1); }
          }
        }
    */
}

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchU_gmp(void) {
    REAL underflowmin;
    REAL one = 1.0;
    unsigned long exp2;
    exp2 = (1UL << (mp_bits_per_limb - 8)) - 1; // 6 seems to be the smallest on amd64 but for safty
    mpf_div_2exp(underflowmin.get_mpf_t(), one.get_mpf_t(), exp2);
    return underflowmin;
}

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchL_gmp(void) {
    REAL maxexp;
    unsigned long exp2;
    exp2 = (1UL << (mp_bits_per_limb - 8)) - 1; // 6 seems to be the smallest on amd64 but for safty
    maxexp = exp2;
    return maxexp;
}

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchO_gmp(void) {
    REAL overflowmax;
    REAL one = 1.0;
    unsigned long exp2;

    exp2 = (1UL << (mp_bits_per_limb - 8)) - 1; // 6 seems to be the smallest on amd64 but for safty
    mpf_mul_2exp(overflowmax.get_mpf_t(), one.get_mpf_t(), exp2);

    return overflowmax;
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
REAL RlamchZ_gmp(void) {
    REAL mtemp = 0.0;
    return mtemp;
}

REAL Rlamch_gmp(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE_gmp();
    if (Mlsame(cmach, "S"))
        return RlamchS_gmp();
    if (Mlsame(cmach, "B"))
        return RlamchB_gmp();
    if (Mlsame(cmach, "P"))
        return RlamchP_gmp();
    if (Mlsame(cmach, "N"))
        return RlamchN_gmp();
    if (Mlsame(cmach, "R"))
        return RlamchR_gmp();
    if (Mlsame(cmach, "M"))
        return RlamchM_gmp();
    if (Mlsame(cmach, "U"))
        return RlamchU_gmp();
    if (Mlsame(cmach, "L"))
        return RlamchL_gmp();
    if (Mlsame(cmach, "O"))
        return RlamchO_gmp();

    Mxerbla("Rlamch", 1);
    return RlamchZ_gmp();
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
qd_real RlamchE_qd(void) {
    // 2^-209 = 1.21e-63
    return qd_real::_eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchS_qd(void) {
    // 2^(-1022+3*53) = 1.626e-260
    return qd_real::_min_normalized;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchB_qd(void) {
    qd_real two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchP_qd(void) {
    qd_real base, eps, prec;

    base = RlamchB_qd();
    eps = RlamchE_qd();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchN_qd(void) {
    return (qd_real)208.0; // 52*4
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchR_qd(void) {
    qd_real mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchM_qd(void) { return qd_real(-1022.0 + 3.0 * 53.0); }

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchU_qd(void) { return qd_real::_min_normalized; }

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchL_qd(void) { return (qd_real)1024.0; }

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchO_qd(void) {
    return qd_real::_max; // approx 1.7976931348623157E+308 in float.h
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchZ_qd(void) {
    qd_real mtemp = 0.0;
    return mtemp;
}

qd_real Rlamch_qd(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE_qd();
    if (Mlsame(cmach, "S"))
        return RlamchS_qd();
    if (Mlsame(cmach, "B"))
        return RlamchB_qd();
    if (Mlsame(cmach, "P"))
        return RlamchP_qd();
    if (Mlsame(cmach, "N"))
        return RlamchN_qd();
    if (Mlsame(cmach, "R"))
        return RlamchR_qd();
    if (Mlsame(cmach, "M"))
        return RlamchM_qd();
    if (Mlsame(cmach, "U"))
        return RlamchU_qd();
    if (Mlsame(cmach, "L"))
        return RlamchL_qd();
    if (Mlsame(cmach, "O"))
        return RlamchO_qd();

    Mxerbla("Rlamch", 1);
    return RlamchZ_qd();
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
dd_real RlamchE_dd(void) {
    // 2^(-52-52) = 2^-104 = 4.93e-32
    return dd_real::_eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchS_dd(void) {
    // 2^(-1022+53) = 2.0042e-292
    return dd_real::_min_normalized;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchB_dd(void) {
    dd_real two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchP_dd(void) {
    dd_real base, eps, prec;

    base = RlamchB_dd();
    eps = RlamchE_dd();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchN_dd(void) {
    return (dd_real)104.0; // 52*2
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchR_dd(void) {
    dd_real mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchM_dd(void) { return dd_real(-1022.0 + 53.0); }

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchU_dd(void) { return dd_real::_min_normalized; }

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchL_dd(void) { return (dd_real)1024.0; }

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchO_dd(void) {
    return dd_real::_max; // approx 1.7976931348623157E+308 in float.h
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
dd_real RlamchZ_dd(void) {
    dd_real mtemp = 0.0;
    return mtemp;
}

dd_real Rlamch_dd(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE_dd();
    if (Mlsame(cmach, "S"))
        return RlamchS_dd();
    if (Mlsame(cmach, "B"))
        return RlamchB_dd();
    if (Mlsame(cmach, "P"))
        return RlamchP_dd();
    if (Mlsame(cmach, "N"))
        return RlamchN_dd();
    if (Mlsame(cmach, "R"))
        return RlamchR_dd();
    if (Mlsame(cmach, "M"))
        return RlamchM_dd();
    if (Mlsame(cmach, "U"))
        return RlamchU_dd();
    if (Mlsame(cmach, "L"))
        return RlamchL_dd();
    if (Mlsame(cmach, "O"))
        return RlamchO_dd();

    Mxerbla("Rlamch", 1);
    return RlamchZ_dd();
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
double RlamchE_double(void) {
    static double eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    // We all know double is the IEEE 754 2008 binary64 format has 53bit significant digits
    for (int i = 0; i < DBL_MANT_DIG; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchS_double(void) {
    return DBL_MIN;

    // IEEE 754 2008 binary64: emin = -1022
    // 2^{-1022} = 2.225073858507201e-308 = DBL_MIN in float.h
    static double eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    // We all know double is the IEEE 754 2008 binary64 format has 53bit significant digits
    for (int i = 0; i < 1022; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchB_double(void) {
    double two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchP_double(void) {
    double base, eps, prec;

    base = RlamchB_double();
    eps = RlamchE_double();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchN_double(void) {
    // IEEE 754 2008 binary64 has 53 (52+1) bit significant digits
    return (double)DBL_MANT_DIG;
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchR_double(void) {
    double mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchM_double(void) {
    // the exponent of IEEE 754 2008 binary64 is -1022.
    // then -1022 + 1 = -1021.
    return (double)DBL_MIN_EXP;
}

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchU_double(void) {
    return DBL_MIN;

    // 2^{-1021-1} minimum exponent
    static double eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    // We all know double is the IEEE 754 2008 binary64 format has 53bit significant digits
    for (int i = 0; i < 1022; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchL_double(void) {
    //+1023 in IEEE 754 2008 binary64
    // then 1023 + 1 = 1024.
    return DBL_MAX_EXP;
}

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchO_double(void) {
    // 1.7976931348623157E+308 in IEEE 754 2008 binary64.
    return DBL_MAX;
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
double RlamchZ_double(void) {
    double mtemp = 0.0;
    return mtemp;
}

double Rlamch_double(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE_double();
    if (Mlsame(cmach, "S"))
        return RlamchS_double();
    if (Mlsame(cmach, "B"))
        return RlamchB_double();
    if (Mlsame(cmach, "P"))
        return RlamchP_double();
    if (Mlsame(cmach, "N"))
        return RlamchN_double();
    if (Mlsame(cmach, "R"))
        return RlamchR_double();
    if (Mlsame(cmach, "M"))
        return RlamchM_double();
    if (Mlsame(cmach, "U"))
        return RlamchU_double();
    if (Mlsame(cmach, "L"))
        return RlamchL_double();
    if (Mlsame(cmach, "O"))
        return RlamchO_double();
    Mxerbla("Rlamch", 1);
    return RlamchZ_double();
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
_Float64x RlamchE__Float64x(void) {
    static _Float64x eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    //_Float64x is the 80-bit extended precision format with 64bit significant digits
    for (int i = 0; i < LDBL_MANT_DIG; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchS__Float64x(void) {
    return LDBL_MIN;

    static _Float64x eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    /* We all know double is the _Float64x 80 bit format has 64bit significant digits and 15 for the exponent (2^15)/2 - 1 = 16383 */
    for (int i = 0; i < 16383; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchB__Float64x(void) {
    _Float64x two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchP__Float64x(void) {
    _Float64x base, eps, prec;

    base = RlamchB__Float64x();
    eps = RlamchE__Float64x();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchN__Float64x(void) {
    // _Float64x with 80 bits has 64 bit significant digits
    return (_Float64x)LDBL_MANT_DIG;
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchR__Float64x(void) {
    _Float64x mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchM__Float64x(void) { return (_Float64x)LDBL_MIN_EXP; }

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchU__Float64x(void) {
    return LDBL_MIN;

    static double eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    /* We all know double is the _Float64x 80 bit format has 64bit significant digits and 15 for the exponent (2^15)/2 - 1 = 16383 */
    for (int i = 0; i < 16383; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchL__Float64x(void) { return LDBL_MAX_EXP; }

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchO__Float64x(void) { return LDBL_MAX; }

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
_Float64x RlamchZ__Float64x(void) {
    _Float64x mtemp = 0.0;
    return mtemp;
}

_Float64x Rlamch__Float64x(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE__Float64x();
    if (Mlsame(cmach, "S"))
        return RlamchS__Float64x();
    if (Mlsame(cmach, "B"))
        return RlamchB__Float64x();
    if (Mlsame(cmach, "P"))
        return RlamchP__Float64x();
    if (Mlsame(cmach, "N"))
        return RlamchN__Float64x();
    if (Mlsame(cmach, "R"))
        return RlamchR__Float64x();
    if (Mlsame(cmach, "M"))
        return RlamchM__Float64x();
    if (Mlsame(cmach, "U"))
        return RlamchU__Float64x();
    if (Mlsame(cmach, "L"))
        return RlamchL__Float64x();
    if (Mlsame(cmach, "O"))
        return RlamchO__Float64x();
    Mxerbla("Rlamch", 1);
    return RlamchZ__Float64x();
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
//"E" denots we always calculate relative machine precision (e).
// where 1+e = 1, minimum of e.
_Float128 RlamchE__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return FLT128_EPSILON;
#else // XXX
    return 1.92592994438723585305597794258492732e-34Q;
#endif
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchS__Float128(void) {
// IEEE 754 2008 binary128: emin = -16382
// 2^{-16382} = 3.36210314311209350626267781732175260e-4932Q
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return FLT128_MIN;
#else // XXX
    return 3.36210314311209350626267781732175260e-4932Q;
#endif

    static _Float128 eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    // We all know double is the IEEE 754 2008 binary128 format has 113bit significant digits
    for (int i = 0; i < 16383; i++) {
        eps = eps / 2.0Q;
    }
    called = 1;
}

//"B" base  = base of the machine
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchB__Float128(void) {
    _Float128 two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchP__Float128(void) {
    _Float128 base, eps, prec;

    base = RlamchB__Float128();
    eps = RlamchE__Float128();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchN__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_MANT_DIG; // 113
#else
    return (_Float128)113;
#endif
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchR__Float128(void) {
    _Float128 mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchM__Float128(void) {
    // the exponent of IEEE 754 2008 binary64 is -16382.
    // then -16382 + 1 = -16381
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return FLT128_MIN_EXP;
#else
    return (-16381);
#endif
}

//"U"
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchU__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return FLT128_MIN;
#else
    return 3.36210314311209350626267781732175260e-4932Q;
#endif

    // 2^{-16382+1} minimum exponent
    static double eps;
    static int called = 0;
    if (called)
        return eps;
    eps = 1.0;
    // We all know double is the IEEE 754 2008 binary128 format has 113bit significant digits
    for (int i = 0; i < 16382; i++) {
        eps = eps / 2.0;
    }
    called = 1;
    return eps;
}

//"L"
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchL__Float128(void) {
    //+16383 in IEEE 754 2008 binary128
    // then 16383 + 1 = 16384
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return FLT128_MAX_EXP;
#else
    return 16384;
#endif
}

//"O"
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchO__Float128(void) {
    // 1.18973149535723176508575932662800702e4932Q in IEEE 754 2008 binary128.
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return FLT128_MAX;
#else
    return 1.18973149535723176508575932662800702e4932Q;
#endif
}

//"Z" :dummy
// cf.http://www.netlib.org/blas/dlamch.f
_Float128 RlamchZ__Float128(void) {
    _Float128 mtemp = 0.0;
    return mtemp;
}

_Float128 Rlamch__Float128(const char *cmach) {
    if (Mlsame(cmach, "E"))
        return RlamchE__Float128();
    if (Mlsame(cmach, "S"))
        return RlamchS__Float128();
    if (Mlsame(cmach, "B"))
        return RlamchB__Float128();
    if (Mlsame(cmach, "P"))
        return RlamchP__Float128();
    if (Mlsame(cmach, "N"))
        return RlamchN__Float128();
    if (Mlsame(cmach, "R"))
        return RlamchR__Float128();
    if (Mlsame(cmach, "M"))
        return RlamchM__Float128();
    if (Mlsame(cmach, "U"))
        return RlamchU__Float128();
    if (Mlsame(cmach, "L"))
        return RlamchL__Float128();
    if (Mlsame(cmach, "O"))
        return RlamchO__Float128();

    Mxerbla("Rlamch", 1);
    return RlamchZ__Float128();
}
#endif

REAL Rlamc3(REAL a, REAL b) { return a + b; }
