/*
 * Copyright (c) 2008-2026
 *      Nakata, Maho
 *      All rights reserved.
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

REAL RlamchE_mpfr(void);
REAL RlamchS_mpfr(void);
REAL RlamchB_mpfr(void);
REAL RlamchP_mpfr(void);
REAL RlamchN_mpfr(void);
REAL RlamchR_mpfr(void);
REAL RlamchM_mpfr(void);
REAL RlamchU_mpfr(void);
REAL RlamchL_mpfr(void);
REAL RlamchO_mpfr(void);
REAL RlamchZ_mpfr(void);

// "E" : relative machine precision (unit roundoff), assume rounding
REAL RlamchE_mpfr(void) {
    REAL one = 1.0;
    const mp_prec_t prec = one.get_prec();

    REAL eps;
    eps.set_prec(prec);

    // eps = 2^(-prec)  (DLAMCH('E') semantics: EPSILON * 0.5)
    eps = mul_2si(one, -static_cast<mp_exp_t>(prec));
    return eps;
}

// "S" : safe minimum such that 1/sfmin does not overflow (netlib logic)
// sfmin = rmin; small = 1/rmax; if (small >= sfmin) sfmin = small*(1+eps)
REAL RlamchS_mpfr(void) {
    REAL one = 1.0;
    const mp_prec_t prec = one.get_prec();

    const REAL eps = RlamchE_mpfr();
    const REAL rmin = RlamchU_mpfr();
    const REAL rmax = RlamchO_mpfr();

    REAL sfmin;
    sfmin.set_prec(prec);

    sfmin = rmin;
    const REAL small = one / rmax;
    if (small >= sfmin) {
        sfmin = small * (one + eps);
    }
    return sfmin;
}

// "B" : base of the machine
REAL RlamchB_mpfr(void) { return REAL(2.0); }

// "P" : precision = eps*base
REAL RlamchP_mpfr(void) { return RlamchE_mpfr() * RlamchB_mpfr(); }

// "N" : number of (base) digits in the mantissa
REAL RlamchN_mpfr(void) {
    REAL one = 1.0;
    return REAL(one.get_prec());
}

// "R" : rounding mode (1.0 when rounding occurs in addition, 0.0 otherwise)
REAL RlamchR_mpfr(void) {
    // DLAMCH('R') reports whether the arithmetic rounds (1) or chops/truncates (0).
    // For MPFR, this is controlled by the current default rounding mode.
    // We interpret MPFR_RNDZ (round toward zero) as chopping; other modes imply rounding.
    const mpfr_rnd_t rnd = mpfr_get_default_rounding_mode();
    return (rnd == MPFR_RNDZ) ? REAL(0.0) : REAL(1.0);
}

// "M" : minimum exponent (emin) used by DLAMCH
REAL RlamchM_mpfr(void) {
    REAL one = 1.0;
    return REAL(one.get_emin());
}

// "U" : underflow threshold rmin = base^(emin-1)
// MPFR: smallest positive number is 0.5 * 2^emin => rmin = 2^(emin-1)
REAL RlamchU_mpfr(void) {
    REAL one = 1.0;
    const mp_prec_t prec = one.get_prec();
    const mp_exp_t emin = one.get_emin();

    REAL rmin;
    rmin.set_prec(prec);
    rmin = mul_2si(one, emin - 1);
    return rmin;
}

// "L" : largest exponent before overflow (emax) used by DLAMCH
REAL RlamchL_mpfr(void) {
    REAL one = 1.0;
    return REAL(one.get_emax());
}

// "O" : overflow threshold rmax = (base^emax) * (1 - eps)
// MPFR: largest value is (1 - eps) * 2^emax
REAL RlamchO_mpfr(void) {
    REAL one = 1.0;
    const mp_prec_t prec = one.get_prec();
    const mp_exp_t emax = one.get_emax();

    const REAL eps = RlamchE_mpfr();

    REAL rmax;
    rmax.set_prec(prec);

    // rmax = (1 - eps) * 2^emax
    rmax = mul_2si(one - eps, emax);
    return rmax;
}

// "Z" : dummy
REAL RlamchZ_mpfr(void) { return REAL(0.0); }

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
REAL RlamchE_gmp(void);
REAL RlamchS_gmp(void);
REAL RlamchB_gmp(void);
REAL RlamchP_gmp(void);
REAL RlamchN_gmp(void);
REAL RlamchR_gmp(void);
REAL RlamchM_gmp(void);
REAL RlamchU_gmp(void);
REAL RlamchL_gmp(void);
REAL RlamchO_gmp(void);
REAL RlamchZ_gmp(void);

// Helper: choose a conservative exponent limit for GMP mpf.
//
// mp_exp_t is a signed type (typically long). We use half of its maximum
// to keep headroom for internal GMP exponent arithmetic and match the
// historical constants (0x3FFFFFFF for 32-bit, 0x3FFFFFFFFFFFFFFF for 64-bit).
static inline mp_bitcnt_t get_max_safe_exponent() {
    const mp_exp_t max_e = std::numeric_limits<mp_exp_t>::max();
    return static_cast<mp_bitcnt_t>(max_e / 2);
}

// "E": Relative machine precision for chopping/truncation.
// Option B for GMP mpf: chopping.
// Define eps as the smallest E such that 1 + E > 1 (gap near 1).
// For this mpf setting, eps = 2^(-prec).
REAL RlamchE_gmp(void) {
    const REAL one = 1.0;
    const mp_bitcnt_t prec = mpf_get_prec(one.get_mpf_t());

    REAL eps;
    eps.set_prec(prec);

    // eps = 2^(-prec)
    mpf_div_2exp(eps.get_mpf_t(), one.get_mpf_t(), prec);

    return eps;
}

// "S": Safe minimum such that 1/sfmin does not overflow.
//
// Netlib DLAMCH logic:
//   sfmin = rmin
//   small = 1/rmax
//   if (small >= sfmin) sfmin = small * (1 + eps)
// Dependency: RlamchE_gmp, RlamchU_gmp, RlamchO_gmp
REAL RlamchS_gmp(void) {
    const REAL one = 1.0;
    const mp_bitcnt_t prec = mpf_get_prec(one.get_mpf_t());

    REAL sfmin;
    sfmin.set_prec(prec);

    const REAL eps = RlamchE_gmp();
    const REAL rmin = RlamchU_gmp();
    const REAL rmax = RlamchO_gmp();

    sfmin = rmin;
    const REAL small = one / rmax;
    if (small >= sfmin) {
        sfmin = small * (one + eps);
    }
    return sfmin;
}

// "B": Base of the machine (radix).
REAL RlamchB_gmp(void) { return REAL(2.0); }

// "P": Precision = eps * base.
// Dependency: RlamchE_gmp, RlamchB_gmp
REAL RlamchP_gmp(void) { return RlamchE_gmp() * RlamchB_gmp(); }

// "N": Number of base-B digits in the mantissa.
REAL RlamchN_gmp(void) {
    const REAL one = 1.0;
    const mp_bitcnt_t prec = mpf_get_prec(one.get_mpf_t());
    return REAL(prec);
}

// "R": Rounding mode indicator.
//      Returns 1.0 when rounding occurs in addition, 0.0 otherwise.
//      GMP mpf arithmetic uses truncation (chopping), not round-to-nearest.
REAL RlamchR_gmp(void) { return REAL(0.0); }

// "M": Minimum exponent (emin).
//      DLAMCH defines rmin = base^(emin-1).
//      With base=2 and rmin = 2^(-E), we have emin = 1 - E.
//
// Dependency: get_max_safe_exponent (not other Rlamch functions)
REAL RlamchM_gmp(void) {
    const mp_bitcnt_t E = get_max_safe_exponent();
    return REAL(1.0) - REAL(E);
}

// "U": Underflow threshold (rmin), minimum positive normalized number.
//      rmin = 2^(-E) where E = get_max_safe_exponent()
// Dependency: None (base quantity)
REAL RlamchU_gmp(void) {
    const REAL one = 1.0;
    const mp_bitcnt_t prec = mpf_get_prec(one.get_mpf_t());
    const mp_bitcnt_t E = get_max_safe_exponent();

    REAL rmin;
    rmin.set_prec(prec);
    // rmin = 2^(-E)
    mpf_div_2exp(rmin.get_mpf_t(), one.get_mpf_t(), E);

    return rmin;
}

// "L": Maximum exponent (emax).
//      DLAMCH defines rmax ~ base^emax.
//      With rmax ~ 2^E, we have emax = E.
REAL RlamchL_gmp(void) {
    const mp_bitcnt_t E = get_max_safe_exponent();
    return REAL(E);
}

// "O": Overflow threshold (rmax), maximum finite number.
//      rmax = (1 - eps) * 2^emax  (DLAMCH model)
//
// Dependency: RlamchE_gmp
REAL RlamchO_gmp(void) {
    const REAL one = 1.0;
    const mp_bitcnt_t prec = mpf_get_prec(one.get_mpf_t());
    const mp_bitcnt_t E = get_max_safe_exponent();

    REAL rmax;
    rmax.set_prec(prec);

    const REAL eps = RlamchE_gmp();
    const REAL mant = one - eps;

    // rmax = (1 - eps) * 2^E
    mpf_mul_2exp(rmax.get_mpf_t(), mant.get_mpf_t(), E);

    return rmax;
}

// "Z": Dummy/placeholder (not used in standard LAPACK).
REAL RlamchZ_gmp(void) { return REAL(0.0); }

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
qd_real RlamchE_qd(void);
qd_real RlamchS_qd(void);
qd_real RlamchB_qd(void);
qd_real RlamchP_qd(void);
qd_real RlamchN_qd(void);
qd_real RlamchR_qd(void);
qd_real RlamchM_qd(void);
qd_real RlamchU_qd(void);
qd_real RlamchL_qd(void);
qd_real RlamchO_qd(void);
qd_real RlamchZ_qd(void);

// "E" : relative machine precision (unit roundoff), assume rounding
qd_real RlamchE_qd(void) {
    // QD library provides a constant for machine precision.
    // We use it as DLAMCH('E') (unit roundoff) for qd_real.
    return qd_real::_eps;
}

// "S" : safe minimum such that 1/sfmin does not overflow (netlib logic)
qd_real RlamchS_qd(void) {
    const qd_real one(1.0);
    const qd_real eps = RlamchE_qd();

    // netlib: sfmin = rmin; small = 1/rmax; if (small >= sfmin) sfmin = small*(1+eps)
    qd_real sfmin = RlamchU_qd();
    const qd_real small = one / RlamchO_qd();

    if (small >= sfmin) {
        sfmin = small * (one + eps);
    }
    return sfmin;
}

// "B" : base of the machine
qd_real RlamchB_qd(void) { return qd_real(2.0); }

// "P" : precision = eps*base
qd_real RlamchP_qd(void) { return RlamchE_qd() * RlamchB_qd(); }

// "N" : number of (base) digits in the mantissa
qd_real RlamchN_qd(void) {
    // Keep N consistent with E/P used above.
    // If E  2^(-N) (unit roundoff) and base=2, then P = 2*E.
    return qd_real(209.0);
}

// "R" : rounding mode (1.0 when rounding occurs in addition, 0.0 otherwise)
qd_real RlamchR_qd(void) {
    // QD operations are performed via IEEE double arithmetic; rounding occurs.
    return qd_real(1.0);
}

// "M" : minimum exponent (emin) used by DLAMCH
qd_real RlamchM_qd(void) {
    // With rmin = _min_normalized  2^(-1022 + 3*53) = 2^-863,
    // DLAMCH defines rmin = base^(emin-1) with base=2, hence emin = -862.
    return qd_real(-862.0);
}

// "U" : underflow threshold (rmin; minimum positive "full-precision" normal)
qd_real RlamchU_qd(void) { return qd_real::_min_normalized; }

// "L" : largest exponent (emax) used by DLAMCH
qd_real RlamchL_qd(void) {
    // Exponent range is effectively that of IEEE double.
    return qd_real(1024.0);
}

// "O" : overflow threshold (rmax; maximum finite)
qd_real RlamchO_qd(void) {
    // Use _max as the overflow threshold.
    // If you prefer extra headroom for intermediate operations, consider _safe_max.
    return qd_real::_max;
}

// "Z" : dummy
qd_real RlamchZ_qd(void) { return qd_real(0.0); }

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
dd_real RlamchE_dd(void);
dd_real RlamchS_dd(void);
dd_real RlamchB_dd(void);
dd_real RlamchP_dd(void);
dd_real RlamchN_dd(void);
dd_real RlamchR_dd(void);
dd_real RlamchM_dd(void);
dd_real RlamchU_dd(void);
dd_real RlamchL_dd(void);
dd_real RlamchO_dd(void);
dd_real RlamchZ_dd(void);

// "E" : relative machine precision (unit roundoff), assume rounding
dd_real RlamchE_dd(void) {
    // dd_real::_eps is already unit roundoff (2^-104)
    return dd_real::_eps;
}

// "S" : safe minimum such that 1/sfmin does not overflow (netlib logic)
dd_real RlamchS_dd(void) {
    const dd_real one(1.0);
    const dd_real eps = RlamchE_dd();
    dd_real sfmin = RlamchU_dd();
    const dd_real small = one / RlamchO_dd();
    if (small >= sfmin) {
        sfmin = small * (one + eps);
    }
    return sfmin;
}

// "B" : base of the machine
dd_real RlamchB_dd(void) { return dd_real(2.0); }

// "P" : precision = eps*base
dd_real RlamchP_dd(void) { return RlamchE_dd() * RlamchB_dd(); }

// "N" : number of (base) digits in the mantissa
dd_real RlamchN_dd(void) {
    return dd_real((double)std::numeric_limits<dd_real>::digits); // 104
}

// "R" : rounding mode (1.0 when rounding occurs in addition)
dd_real RlamchR_dd(void) { return dd_real(1.0); }

// "M" : minimum exponent (emin) where rmin = base^(emin-1)
dd_real RlamchM_dd(void) {
    int emin = 0;
    (void)std::frexp(dd_real::_min_normalized, &emin);
    return dd_real((double)emin);
}

// "U" : underflow threshold (rmin; minimum positive normalized)
dd_real RlamchU_dd(void) { return dd_real::_min_normalized; }

// "L" : largest exponent (emax) where rmax  base^emax
dd_real RlamchL_dd(void) { return dd_real(1024.0); }

// "O" : overflow threshold (rmax; maximum finite)
dd_real RlamchO_dd(void) { return dd_real::_max; }

// "Z" : dummy
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
double RlamchE_double(void);
double RlamchS_double(void);
double RlamchB_double(void);
double RlamchP_double(void);
double RlamchN_double(void);
double RlamchR_double(void);
double RlamchM_double(void);
double RlamchU_double(void);
double RlamchL_double(void);
double RlamchO_double(void);
double RlamchZ_double(void);

// "E" : relative machine precision (unit roundoff), assume rounding
double RlamchE_double(void) {
    static double eps;
    static int called = 0;

    if (called)
        return eps;

    eps = 1.0;
    // IEEE 754 binary64 has DBL_MANT_DIG (=53) significant bits.
    // This computes eps = 2^(-DBL_MANT_DIG) = 2^-53 (unit roundoff).
    for (int i = 0; i < DBL_MANT_DIG; i++) {
        eps = eps / 2.0;
    }

    called = 1;
    return eps;
}

// "S" : safe minimum such that 1/sfmin does not overflow
double RlamchS_double(void) {
    // For IEEE binary64, DBL_MIN is safe since 1/DBL_MIN does not overflow.
    return DBL_MIN;
}

// "B" : base of the machine
double RlamchB_double(void) { return 2.0; }

// "P" : precision = eps*base
double RlamchP_double(void) { return RlamchE_double() * RlamchB_double(); }

// "N" : number of digits in mantissa
double RlamchN_double(void) {
    return (double)DBL_MANT_DIG; // 53
}

// "R" : 1.0 when rounding occurs in addition, 0.0 otherwise
double RlamchR_double(void) { return 1.0; }

// "M" : minimum exponent (emin) used by DLAMCH
double RlamchM_double(void) {
    // For IEEE binary64, DBL_MIN_EXP is -1021 (so emin-1 = -1022 gives DBL_MIN).
    return (double)DBL_MIN_EXP;
}

// "U" : underflow threshold (rmin; minimum positive normal)
double RlamchU_double(void) { return DBL_MIN; }

// "L" : largest exponent (emax) used by DLAMCH
double RlamchL_double(void) {
    return (double)DBL_MAX_EXP; // 1024
}

// "O" : overflow threshold (rmax; maximum finite)
double RlamchO_double(void) { return DBL_MAX; }

// "Z" : dummy
double RlamchZ_double(void) { return 0.0; }
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
_Float64x RlamchE__Float64x(void);
_Float64x RlamchS__Float64x(void);
_Float64x RlamchB__Float64x(void);
_Float64x RlamchP__Float64x(void);
_Float64x RlamchN__Float64x(void);
_Float64x RlamchR__Float64x(void);
_Float64x RlamchM__Float64x(void);
_Float64x RlamchU__Float64x(void);
_Float64x RlamchL__Float64x(void);
_Float64x RlamchO__Float64x(void);
_Float64x RlamchZ__Float64x(void);

// "E" : relative machine precision (unit roundoff), assume rounding
_Float64x RlamchE__Float64x(void) {
    static _Float64x eps;
    static int called = 0;

    if (called)
        return eps;

    eps = (_Float64x)1.0;

#if defined(FLT64X_MANT_DIG)
    // Compute eps = 2^(-FLT64X_MANT_DIG)
    for (int i = 0; i < FLT64X_MANT_DIG; i++)
        eps /= (_Float64x)2.0;
#else
    // Fallback: assume _Float64x matches long double characteristics.
    for (int i = 0; i < LDBL_MANT_DIG; i++)
        eps /= (_Float64x)2.0;
#endif

    called = 1;
    return eps;
}

// "S" : safe minimum such that 1/sfmin does not overflow
_Float64x RlamchS__Float64x(void) {
#if defined(FLT64X_MIN)
    return (_Float64x)FLT64X_MIN;
#else
    // Fallback: assume _Float64x matches long double characteristics.
    return (_Float64x)LDBL_MIN;
#endif
}

// "B" : base of the machine
_Float64x RlamchB__Float64x(void) { return (_Float64x)2.0; }

// "P" : precision = eps*base
_Float64x RlamchP__Float64x(void) { return RlamchE__Float64x() * RlamchB__Float64x(); }

// "N" : number of digits in mantissa
_Float64x RlamchN__Float64x(void) {
#if defined(FLT64X_MANT_DIG)
    return (_Float64x)FLT64X_MANT_DIG;
#else
    return (_Float64x)LDBL_MANT_DIG;
#endif
}

// "R" : 1.0 when rounding occurs in addition, 0.0 otherwise
_Float64x RlamchR__Float64x(void) { return (_Float64x)1.0; }

// "M" : minimum exponent (emin) used by DLAMCH
_Float64x RlamchM__Float64x(void) {
#if defined(FLT64X_MIN_EXP)
    return (_Float64x)FLT64X_MIN_EXP;
#else
    return (_Float64x)LDBL_MIN_EXP;
#endif
}

// "U" : underflow threshold (rmin; minimum positive normal)
_Float64x RlamchU__Float64x(void) {
#if defined(FLT64X_MIN)
    return (_Float64x)FLT64X_MIN;
#else
    return (_Float64x)LDBL_MIN;
#endif
}

// "L" : largest exponent (emax) used by DLAMCH
_Float64x RlamchL__Float64x(void) {
#if defined(FLT64X_MAX_EXP)
    return (_Float64x)FLT64X_MAX_EXP;
#else
    return (_Float64x)LDBL_MAX_EXP;
#endif
}

// "O" : overflow threshold (rmax; maximum finite)
_Float64x RlamchO__Float64x(void) {
#if defined(FLT64X_MAX)
    return (_Float64x)FLT64X_MAX;
#else
    return (_Float64x)LDBL_MAX;
#endif
}

// "Z" : dummy
_Float64x RlamchZ__Float64x(void) { return (_Float64x)0.0; }

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
_Float128 RlamchE__Float128(void);
_Float128 RlamchS__Float128(void);
_Float128 RlamchB__Float128(void);
_Float128 RlamchP__Float128(void);
_Float128 RlamchN__Float128(void);
_Float128 RlamchR__Float128(void);
_Float128 RlamchM__Float128(void);
_Float128 RlamchU__Float128(void);
_Float128 RlamchL__Float128(void);
_Float128 RlamchO__Float128(void);
_Float128 RlamchZ__Float128(void);

// "E" : relative machine precision (unit roundoff), assume rounding
// DLAMCH('E') corresponds to 0.5 * (gap near 1) when rounding occurs.
_Float128 RlamchE__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_EPSILON * (_Float128)0.5;
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    // 2^-113
    return 9.6296497219361792652798897129246365926905e-35L;
#else
    // 2^-113
    return 9.6296497219361792652798897129246365926905e-35Q;
#endif
}

// "S" : safe minimum such that 1/sfmin does not overflow (netlib logic)
// sfmin = rmin; small = 1/rmax; if (small >= sfmin) sfmin = small*(1+eps)
_Float128 RlamchS__Float128(void) {
    const _Float128 one = (_Float128)1.0;
    const _Float128 eps = RlamchE__Float128();

    _Float128 sfmin = RlamchU__Float128();
    const _Float128 small = one / RlamchO__Float128();

    if (small >= sfmin) {
        sfmin = small * (one + eps);
    }
    return sfmin;
}

// "B" : base of the machine
_Float128 RlamchB__Float128(void) { return (_Float128)2.0; }

// "P" : precision = eps*base
_Float128 RlamchP__Float128(void) { return RlamchE__Float128() * RlamchB__Float128(); }

// "N" : number of digits in mantissa
_Float128 RlamchN__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_MANT_DIG; // 113
#else
    return (_Float128)113;
#endif
}

// "R" : 1.0 when rounding occurs in addition, 0.0 otherwise
_Float128 RlamchR__Float128(void) { return (_Float128)1.0; }

// "M" : minimum exponent (emin) used by DLAMCH
// For IEEE binary128, FLT128_MIN is 2^(-16382) and DLAMCH uses rmin = base^(emin-1).
// Thus emin = -16381 so that emin-1 = -16382.
_Float128 RlamchM__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_MIN_EXP; // -16381
#else
    return (_Float128)(-16381);
#endif
}

// "U" : underflow threshold (rmin; minimum positive normal)
_Float128 RlamchU__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_MIN;
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    return 3.3621031431120935062626778173217526025981e-4932L;
#else
    return 3.3621031431120935062626778173217526025981e-4932Q;
#endif
}

// "L" : largest exponent (emax) used by DLAMCH
_Float128 RlamchL__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_MAX_EXP; // 16384
#else
    return (_Float128)16384;
#endif
}

// "O" : overflow threshold (rmax; maximum finite)
_Float128 RlamchO__Float128(void) {
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    return (_Float128)FLT128_MAX;
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    return 1.1897314953572317650857593266280070161965e+4932L;
#else
    return 1.1897314953572317650857593266280070161965e+4932Q;
#endif
}

// "Z" : dummy
_Float128 RlamchZ__Float128(void) { return (_Float128)0.0; }

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
