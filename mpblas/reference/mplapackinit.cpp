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
#include <gmpfrxx_mkII/detail/gmp_default_context.hpp>

namespace {

gmpxx_mkII_default_context_v1 make_initial_gmp_context() noexcept {
    gmpxx_mkII_default_context_v1 context{};
    context.abi_version = gmpfrxx_mkII::detail::gmp_default_context_abi_version;
    context.struct_size = sizeof(context);
    context.mpf_precision_bits =
        static_cast<std::uint64_t>(gmpfrxx_mkII::detail::read_default_mpf_precision_from_environment());
    return context;
}

const gmpxx_mkII_default_context_v1 &initial_gmp_context() noexcept {
    static const gmpxx_mkII_default_context_v1 context = make_initial_gmp_context();
    return context;
}

gmpxx_mkII_default_context_v1 &thread_gmp_context() noexcept {
    static thread_local gmpxx_mkII_default_context_v1 context = initial_gmp_context();
    return context;
}

int gmp_context_provider_token = 0;

} // namespace

extern "C" GMPXX_MKII_API void gmpxx_mkII_get_current_default_context_v1(gmpxx_mkII_default_context_v1 *out) noexcept {
    if (out != nullptr) {
        *out = thread_gmp_context();
    }
}

extern "C" GMPXX_MKII_API void gmpxx_mkII_set_thread_default_context_v1(const gmpxx_mkII_default_context_v1 *context) noexcept {
    if (context != nullptr) {
        gmpfrxx_mkII::detail::validate_default_context_or_abort(*context);
        thread_gmp_context() = *context;
    }
}

extern "C" GMPXX_MKII_API void gmpxx_mkII_reset_thread_default_context_v1() noexcept {
    thread_gmp_context() = initial_gmp_context();
}

extern "C" GMPXX_MKII_API const void *gmpxx_mkII_default_context_provider_token_v1() noexcept {
    return &gmp_context_provider_token;
}

extern "C" GMPXX_MKII_API int gmpxx_mkII_default_context_mode_v1() noexcept {
    return GMPXX_MKII_DEFAULT_CONTEXT_EXTERNAL_PROVIDER;
}

void __attribute__((constructor)) mplapack_initialize_gmp(void);
void __attribute__((destructor)) mplapack_finalize_gmp(void);
void mplapack_initialize_gmp(void) {
    // The GMP backend DSO owns the gmpfrxx_mkII default-precision context.
}
void mplapack_finalize_gmp(void) {}
#endif

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
void __attribute__((constructor)) mplapack_initialize_mpfr(void);
void __attribute__((destructor)) mplapack_finalize_mpfr(void);
void mplapack_initialize_mpfr(void) {
    // MPFR/MPC default state is initialized and owned by gmpfrxx_mkII.
}
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
