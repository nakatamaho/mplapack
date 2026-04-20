/*
 * Copyright (c) 2012-2021
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

#define ___MPLAPACK_MPLAPACK_INIT___

#include <mpblas.h>
#include <stdio.h>

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void __attribute__((constructor)) mplapack_initialize_gmp(void);
void __attribute__((destructor)) mplapack_finalize_gmp(void);
void mplapack_initialize_gmp(void) {
    mpf_set_default_prec(___MPLAPACK_GMP_DEFAULT_PRECISION___);
    char *p = getenv("MPLAPACK_GMP_PRECISION");
    if (p) {
        mpf_set_default_prec(atoi(p));
    }
}
void mplapack_finalize_gmp(void) {
    // no finalization needed
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
// Parse an integer exponent from an environment variable.
// Abort on invalid syntax / overflow.
static mpfr_exp_t parse_env_exp_or_abort(const char *name, bool *present) {
    const char *p = getenv(name);
    if (!p || !*p) {
        *present = false;
        return 0;
    }
    *present = true;

    errno = 0;
    char *end = nullptr;
    long v = strtol(p, &end, 10);

    if (errno != 0 || end == p || *end != '\0') {
        std::fprintf(stderr, "mplapack: invalid %s='%s'\n", name, p);
        std::abort();
    }
    return static_cast<mpfr_exp_t>(v);
}

// Parse an MPFR precision from an environment variable.
// Abort on invalid syntax / out-of-range values.
static mpfr_prec_t parse_env_prec_or_abort(const char *name, bool *present) {
    const char *p = getenv(name);
    if (!p || !*p) {
        *present = false;
        return 0;
    }
    *present = true;

    errno = 0;
    char *end = nullptr;
    long v = strtol(p, &end, 10);

    if (errno != 0 || end == p || *end != '\0') {
        std::fprintf(stderr, "mplapack: invalid %s='%s'\n", name, p);
        std::abort();
    }

    // MPFR requires precision within [MPFR_PREC_MIN, MPFR_PREC_MAX].
    if (v < (long)MPFR_PREC_MIN || v > (long)MPFR_PREC_MAX) {
        std::fprintf(stderr, "mplapack: %s out of range (%ld not in [%ld, %ld])\n", name, v, (long)MPFR_PREC_MIN, (long)MPFR_PREC_MAX);
        std::abort();
    }

    return static_cast<mpfr_prec_t>(v);
}

// Optionally set MPFR exponent range (emin/emax) from environment variables.
// - MPLAPACK_MPFR_EMIN
// - MPLAPACK_MPFR_EMAX
// Abort on invalid settings to avoid silent misconfiguration.
static bool mpfr_set_exp_range_from_env_or_abort(void) {
    bool have_emin = false, have_emax = false;
    const mpfr_exp_t env_emin = parse_env_exp_or_abort("MPLAPACK_MPFR_EMIN", &have_emin);
    const mpfr_exp_t env_emax = parse_env_exp_or_abort("MPLAPACK_MPFR_EMAX", &have_emax);
    if (!have_emin && !have_emax)
        return false;

    const mpfr_exp_t cur_emin = mpfr_get_emin();
    const mpfr_exp_t cur_emax = mpfr_get_emax();

    const mpfr_exp_t target_emin = have_emin ? env_emin : cur_emin;
    const mpfr_exp_t target_emax = have_emax ? env_emax : cur_emax;

    const mpfr_exp_t emin_min = mpfr_get_emin_min();
    const mpfr_exp_t emin_max = mpfr_get_emin_max();
    const mpfr_exp_t emax_min = mpfr_get_emax_min();
    const mpfr_exp_t emax_max = mpfr_get_emax_max();

    if (target_emin < emin_min || target_emin > emin_max) {
        std::fprintf(stderr, "mplapack: %s out of range (%ld not in [%ld, %ld])\n", "MPLAPACK_MPFR_EMIN", (long)target_emin, (long)emin_min, (long)emin_max);
        std::abort();
    }
    if (target_emax < emax_min || target_emax > emax_max) {
        std::fprintf(stderr, "mplapack: %s out of range (%ld not in [%ld, %ld])\n", "MPLAPACK_MPFR_EMAX", (long)target_emax, (long)emax_min, (long)emax_max);
        std::abort();
    }
    if (target_emin >= target_emax) {
        std::fprintf(stderr, "mplapack: invalid exp range (emin=%ld, emax=%ld): emin must be < emax\n", (long)target_emin, (long)target_emax);
        std::abort();
    }

    // Avoid transient (emin >= emax) during updates.
    if (target_emax < cur_emin) {
        if (mpfr_set_emin(target_emin) != 0) {
            std::fprintf(stderr, "mplapack: mpfr_set_emin(%ld) failed\n", (long)target_emin);
            std::abort();
        }
        if (mpfr_set_emax(target_emax) != 0) {
            std::fprintf(stderr, "mplapack: mpfr_set_emax(%ld) failed\n", (long)target_emax);
            std::abort();
        }
    } else if (target_emin > cur_emax) {
        if (mpfr_set_emax(target_emax) != 0) {
            std::fprintf(stderr, "mplapack: mpfr_set_emax(%ld) failed\n", (long)target_emax);
            std::abort();
        }
        if (mpfr_set_emin(target_emin) != 0) {
            std::fprintf(stderr, "mplapack: mpfr_set_emin(%ld) failed\n", (long)target_emin);
            std::abort();
        }
    } else {
        if (mpfr_set_emin(target_emin) != 0) {
            std::fprintf(stderr, "mplapack: mpfr_set_emin(%ld) failed\n", (long)target_emin);
            std::abort();
        }
        if (mpfr_set_emax(target_emax) != 0) {
            std::fprintf(stderr, "mplapack: mpfr_set_emax(%ld) failed\n", (long)target_emax);
            std::abort();
        }
    }
    return true;
}

void __attribute__((constructor)) mplapack_initialize_mpfr(void);
void mplapack_initialize_mpfr(void) {
    // Apply MPFR exponent range from environment variables (if set).
    mpfr_set_exp_range_from_env_or_abort();

    // Default settings.
    mpreal::default_rnd = mpfr_get_default_rounding_mode();
    mpreal::default_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpreal::default_base = 10;
    mpreal::double_bits = -1;
    mpcomplex::default_rnd = MPC_RND(mpfr_get_default_rounding_mode(), mpfr_get_default_rounding_mode());
    mpcomplex::default_real_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpcomplex::default_imag_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpcomplex::default_base = 10;
    mpcomplex::double_bits = -1;

    // Optional: set global precision for both mpreal and mpcomplex.
    {
        bool have = false;
        const mpfr_prec_t prec = parse_env_prec_or_abort("MPLAPACK_MPFR_PRECISION", &have);
        if (have) {
            mpreal::default_prec = prec;
            mpcomplex::default_real_prec = prec;
            mpcomplex::default_imag_prec = prec;
        }
    }

    // Optional: override mpcomplex real/imag precision independently.
    // These take precedence over MPLAPACK_MPFR_PRECISION if both are set.
    {
        bool have_r = false;
        const mpfr_prec_t rprec = parse_env_prec_or_abort("MPLAPACK_MPC_REAL_PRECISION", &have_r);
        if (have_r) {
            mpcomplex::default_real_prec = rprec;
        }

        bool have_i = false;
        const mpfr_prec_t iprec = parse_env_prec_or_abort("MPLAPACK_MPC_IMAG_PRECISION", &have_i);
        if (have_i) {
            mpcomplex::default_imag_prec = iprec;
        }
    }
}
void __attribute__((destructor)) mplapack_finalize_mpfr(void);
void mplapack_finalize_mpfr(void) {}
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
void __attribute__((constructor)) mplapack_initialize_qd(void);
void __attribute__((destructor)) mplapack_finalize_qd(void);
static unsigned int oldcw_qd;
void mplapack_initialize_qd(void) { fpu_fix_start(&oldcw_qd); }
void mplapack_finalize_qd(void) { fpu_fix_end(&oldcw_qd); }
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
void __attribute__((constructor)) mplapack_initialize_dd(void);
void __attribute__((destructor)) mplapack_finalize_dd(void);
static unsigned int oldcw_dd;
void mplapack_initialize_dd(void) { fpu_fix_start(&oldcw_dd); }
void mplapack_finalize_dd(void) { fpu_fix_end(&oldcw_dd); }
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
void __attribute__((constructor)) mplapack_initialize_double(void);
void __attribute__((destructor)) mplapack_finalize_double(void);
void mplapack_initialize_double(void) {
    // no initializization needed
}

void mplapack_finalize_double(void) {
    // no finalization needed
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
void __attribute__((constructor)) mplapack_initialize_binary80(void);
void __attribute__((destructor)) mplapack_finalize_binary80(void);
void mplapack_initialize_binary80(void) {
    // no initializization needed
}

void mplapack_finalize_binary80(void) {
    // no finalization needed
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
void __attribute__((constructor)) mplapack_initialize_binary128(void);
void __attribute__((destructor)) mplapack_finalize_binary128(void);
void mplapack_initialize_binary128(void) {
    // no initializization needed
}

void mplapack_finalize_binary128(void) {
    // no finalization needed
}
#endif
