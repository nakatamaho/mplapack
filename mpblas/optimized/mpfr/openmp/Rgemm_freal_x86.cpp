/*
 * Copyright (c) 2026
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

// Rgemm_freal_engine.h with the x86-64 mulx/adcx/adox multiply of
// ../lgpl/Rfreal_host.h: the functions of this file are compiled for
// BMI2 + ADX (and -O3) whatever the flags of the library, and
// Rgemm_freal_omp.cpp calls them only on a CPU that has both.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <mpblas_mpfr.h>
#include "Rgemm_freal.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(MPLAPACK_NO_FREAL_X86)
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3")
#pragma GCC target("bmi2,adx")
#endif
#define MPFREAL_FORCE_X86_MULX 1
#define MPFREAL_NS mplapack_freal_x86
#define MPFREAL_ENTRY Rgemm_freal_x86
#include "Rgemm_freal_engine.h"
bool Rgemm_freal_x86_available() { return true; }
#else
bool Rgemm_freal_x86(int, bool, bool, long, long, long, const mpfr_class &, const mpfr_class *, long, const mpfr_class *, long, const mpfr_class &, int, mpfr_class *, long, bool) { return false; }
bool Rgemm_freal_x86_available() { return false; }
#endif
