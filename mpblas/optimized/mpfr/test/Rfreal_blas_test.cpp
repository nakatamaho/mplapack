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

// The fixed-precision MPFR BLAS (openmp/Rfreal_blas.h) against the reference
// code (<name>_reference_cpu), bit for bit -- signed zeros and precisions
// included -- for every routine, several precisions, increments, options and
// scalars, with operands that contain signed zeros, small integers (exact
// cancellations, ties) and spread exponents, with OpenMP on and off; that
// the kernels ran (mplapack_freal_done) where they should, and that the
// calls they decline (NaN, Inf, another precision, a complex division,
// overflow, rounding mode) give the MPFR result.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <mpblas_mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../openmp/Rfreal_state.h"

typedef mpfr_class R;
typedef mpc_class C;
typedef mplapackint I;

#define MPLAPACK_FREAL_TEST_REF(name, ret, params) ret name##_reference_cpu params;
// the reference routines of the library (../../common/blas/<name>_omp.cpp)
MPLAPACK_FREAL_TEST_REF(Rsyrk, void, (const char *, const char *, I, I, R, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rsyr2k, void, (const char *, const char *, I, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rsymm, void, (const char *, const char *, I, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtrmm, void, (const char *, const char *, const char *, const char *, I, I, R, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtrsm, void, (const char *, const char *, const char *, const char *, I, I, R, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rgemmtr, void, (const char *, const char *, const char *, I, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rgemv, void, (const char *, I, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rsymv, void, (const char *, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtrmv, void, (const char *, const char *, const char *, I, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtrsv, void, (const char *, const char *, const char *, I, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rger, void, (I, I, R, R *, I, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rsyr, void, (const char *, I, R, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rsyr2, void, (const char *, I, R, R *, I, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rgbmv, void, (const char *, I, I, I, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rsbmv, void, (const char *, I, I, R, R *, I, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rspmv, void, (const char *, I, R, R *, R *, I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rspr, void, (const char *, I, R, R *, I, R *))
MPLAPACK_FREAL_TEST_REF(Rspr2, void, (const char *, I, R, R *, I, R *, I, R *))
MPLAPACK_FREAL_TEST_REF(Rtbmv, void, (const char *, const char *, const char *, I, I, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtbsv, void, (const char *, const char *, const char *, I, I, R *, I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtpmv, void, (const char *, const char *, const char *, I, R *, R *, I))
MPLAPACK_FREAL_TEST_REF(Rtpsv, void, (const char *, const char *, const char *, I, R *, R *, I))
MPLAPACK_FREAL_TEST_REF(Rscal, void, (I, R, R *, I))
MPLAPACK_FREAL_TEST_REF(Rrot, void, (I, R *, I, R *, I, R, R))
MPLAPACK_FREAL_TEST_REF(Rasum, R, (I, R *, I))
MPLAPACK_FREAL_TEST_REF(Rdot, R, (I, R *, I, R *, I))
void Raxpy_ref(mplapackint n, mpfr_class da, mpfr_class *dx, mplapackint incx, mpfr_class *dy, mplapackint incy);
MPLAPACK_FREAL_TEST_REF(Cgemm, void, (const char *, const char *, I, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Csymm, void, (const char *, const char *, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Chemm, void, (const char *, const char *, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctrmm, void, (const char *, const char *, const char *, const char *, I, I, C, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctrsm, void, (const char *, const char *, const char *, const char *, I, I, C, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Csyrk, void, (const char *, const char *, I, I, C, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Cherk, void, (const char *, const char *, I, I, R, C *, I, R, C *, I))
MPLAPACK_FREAL_TEST_REF(Csyr2k, void, (const char *, const char *, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Cher2k, void, (const char *, const char *, I, I, C, C *, I, C *, I, R, C *, I))
MPLAPACK_FREAL_TEST_REF(Cgemmtr, void, (const char *, const char *, const char *, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Cgemv, void, (const char *, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Cgbmv, void, (const char *, I, I, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Chemv, void, (const char *, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Chbmv, void, (const char *, I, I, C, C *, I, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Chpmv, void, (const char *, I, C, C *, C *, I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(Cgeru, void, (I, I, C, C *, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Cgerc, void, (I, I, C, C *, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Cher, void, (const char *, I, R, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Cher2, void, (const char *, I, C, C *, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Chpr, void, (const char *, I, R, C *, I, C *))
MPLAPACK_FREAL_TEST_REF(Chpr2, void, (const char *, I, C, C *, I, C *, I, C *))
MPLAPACK_FREAL_TEST_REF(Ctrmv, void, (const char *, const char *, const char *, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctbmv, void, (const char *, const char *, const char *, I, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctpmv, void, (const char *, const char *, const char *, I, C *, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctrsv, void, (const char *, const char *, const char *, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctbsv, void, (const char *, const char *, const char *, I, I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Ctpsv, void, (const char *, const char *, const char *, I, C *, C *, I))
MPLAPACK_FREAL_TEST_REF(Caxpy, void, (I, C, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Cscal, void, (I, C, C *, I))
MPLAPACK_FREAL_TEST_REF(CRscal, void, (I, R, C *, I))
MPLAPACK_FREAL_TEST_REF(CRrot, void, (I, C *, I, C *, I, R, R))
MPLAPACK_FREAL_TEST_REF(Cdotc, C, (I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(Cdotu, C, (I, C *, I, C *, I))
MPLAPACK_FREAL_TEST_REF(RCasum, R, (I, C *, I))

namespace {

std::mt19937_64 rng(20261123);
gmp_randstate_t gstate;
long failures = 0, checks = 0;
long prec = 512;
std::map<std::string, long> used, calls;

// zeros of both signs, small integers, random mantissas with exponents near
// each other or spread; "tame" keeps the magnitudes in [1/2, 2] (solves)
R rv(bool tame = false)
{
    const int r = (int)(rng() % 16);
    if (!tame) {
        if (r == 0)
            return R(0.0);
        if (r == 1)
            return -R(0.0);
        if (r < 5)
            return R((double)((long)(rng() % 9) - 4));
    }
    R x;
    mpfr_urandomb(x.get_mpfr_t(), gstate);
    if (rng() % 2)
        x = -x;
    long e = tame ? 0 : (long)(rng() % 7) - 3;
    if (!tame && rng() % 10 == 0)
        e = (long)(rng() % 200) - 100;
    mpfr_mul_2si(x.get_mpfr_t(), x.get_mpfr_t(), e, MPFR_RNDN);
    if (rng() % 6 == 0) { // short mantissa: exact cancellations, ties
        mpfr_prec_round(x.get_mpfr_t(), 1 + (long)(rng() % 12), MPFR_RNDN);
        mpfr_prec_round(x.get_mpfr_t(), prec, MPFR_RNDN);
    }
    return x;
}
C cv()
{
    const int r = (int)(rng() % 8);
    if (r == 0)
        return C(rv(), R(0.0)); // real (mul_real of MPC)
    if (r == 1)
        return C(-R(0.0), rv()); // imaginary (mul_imag)
    return C(rv(), rv());
}
R rs(int kind) { return kind == 0 ? R(0.0) : kind == 1 ? R(1.0) : kind == 2 ? R(-0.0) : rv(); }
C cs(int kind) { return kind == 0 ? C(0.0, 0.0) : kind == 1 ? C(1.0, 0.0) : kind == 2 ? C(R(0.0), -R(0.0)) : cv(); }

bool same(const R &a, const R &b)
{
    mpfr_srcptr x = a.get_mpfr_t(), y = b.get_mpfr_t();
    if (mpfr_get_prec(x) != mpfr_get_prec(y))
        return false;
    if (mpfr_nan_p(x) || mpfr_nan_p(y))
        return mpfr_nan_p(x) && mpfr_nan_p(y);
    return mpfr_equal_p(x, y) && mpfr_signbit(x) == mpfr_signbit(y);
}
bool same(const C &a, const C &b) { return same(R(a.real()), R(b.real())) && same(R(a.imag()), R(b.imag())); }

struct bufs {
    std::vector<R> r[4];
    std::vector<C> c[4];
    R rres;
    C cres;
};

// runs f(b, true) (the library: the kernels where they apply) and f(b, false)
// (the reference) on copies of the same operands; compares every operand
void run_case(const std::string &name, const std::string &what, std::initializer_list<long> rsz, std::initializer_list<long> csz, bool tame,
              const std::function<void(bufs &, bool)> &f, bool expect_used = true)
{
    bufs a;
    int t = 0;
    for (long s : rsz) {
        a.r[t].resize((size_t)s);
        for (auto &x : a.r[t])
            x = rv(tame);
        t++;
    }
    t = 0;
    for (long s : csz) {
        a.c[t].resize((size_t)s);
        for (auto &x : a.c[t])
            x = tame ? C(rv(true), rv(true)) : cv();
        t++;
    }
    bufs g = a, r = a;
    const unsigned long d0 = mplapack_freal_done.load();
    f(g, true);
    const bool u = mplapack_freal_done.load() != d0;
    f(r, false);
    checks++;
    calls[name]++;
    if (u)
        used[name]++;
    bool ok = same(g.rres, r.rres) && same(g.cres, r.cres);
    for (int q = 0; q < 4 && ok; q++) {
        for (size_t i = 0; i < g.r[q].size() && ok; i++)
            ok = same(g.r[q][i], r.r[q][i]);
        for (size_t i = 0; i < g.c[q].size() && ok; i++)
            ok = same(g.c[q][i], r.c[q][i]);
    }
    if (!ok || (expect_used && !u && !tame)) {
        if (!ok || std::getenv("FREAL_TEST_VERBOSE"))
            std::printf("%s %ld %s: %s\n", ok ? "NOT USED" : "FAIL", prec, name.c_str(), what.c_str());
        if (!ok)
            failures++;
    }
}

std::string tag(std::initializer_list<const char *> o, std::initializer_list<long> v)
{
    std::string s;
    for (const char *x : o)
        s += std::string(x) + " ";
    for (long x : v)
        s += std::to_string(x) + " ";
    return s;
}

#define CALL(name, args)                                                                                                                                       \
    do {                                                                                                                                                       \
        if (lib)                                                                                                                                               \
            name args;                                                                                                                                         \
        else                                                                                                                                                   \
            name##_reference_cpu args;                                                                                                                         \
    } while (0)

const char *UL[] = {"U", "L"}, *TR[] = {"N", "T", "C"}, *SD[] = {"L", "R"}, *DG[] = {"N", "U"};

void level3(long m, long n, long k)
{
    for (int ak = 0; ak < 4; ak++)
        for (int bk = 0; bk < 4; bk++) {
            if ((ak == 1 || ak == 2) && bk == 2)
                continue;
            const R ra = rs(ak), rb = rs(bk);
            const C ca = cs(ak), cb = cs(bk);
            for (const char *u : UL)
                for (const char *t : TR) {
                    const bool nt = t[0] == 'N';
                    const long ld = (nt ? n : k) + 1, ca_ = nt ? k : n, ldc = n + 2;
                    const std::string w = tag({u, t}, {m, n, k, ak, bk});
                    if (t[0] != 'C') {
                        run_case("Rsyrk", w, {ld * ca_, ldc * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rsyrk, (u, t, n, k, ra, b.r[0].data(), ld, rb, b.r[1].data(), ldc)); });
                        run_case("Rsyr2k", w, {ld * ca_, ld * ca_, ldc * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rsyr2k, (u, t, n, k, ra, b.r[0].data(), ld, b.r[1].data(), ld, rb, b.r[2].data(), ldc)); });
                        run_case("Csyrk", w, {}, {ld * ca_, ldc * n}, false, [&](bufs &b, bool lib) { CALL(Csyrk, (u, t, n, k, ca, b.c[0].data(), ld, cb, b.c[1].data(), ldc)); });
                        run_case("Csyr2k", w, {}, {ld * ca_, ld * ca_, ldc * n}, false, [&](bufs &b, bool lib) { CALL(Csyr2k, (u, t, n, k, ca, b.c[0].data(), ld, b.c[1].data(), ld, cb, b.c[2].data(), ldc)); });
                    }
                    if (t[0] != 'T') {
                        run_case("Cherk", w, {}, {ld * ca_, ldc * n}, false, [&](bufs &b, bool lib) { CALL(Cherk, (u, t, n, k, ra, b.c[0].data(), ld, rb, b.c[1].data(), ldc)); });
                        run_case("Cher2k", w, {}, {ld * ca_, ld * ca_, ldc * n}, false, [&](bufs &b, bool lib) { CALL(Cher2k, (u, t, n, k, ca, b.c[0].data(), ld, b.c[1].data(), ld, rb, b.c[2].data(), ldc)); });
                    }
                    for (const char *t2 : TR) {
                        const bool nt2 = t2[0] == 'N';
                        const long lda = (nt ? n : k) + 1, ldb = (nt2 ? k : n) + 1;
                        const std::string w2 = tag({u, t, t2}, {n, k, ak, bk});
                        if (t[0] != 'C' && t2[0] != 'C')
                            run_case("Rgemmtr", w2, {lda * (nt ? k : n), ldb * (nt2 ? n : k), ldc * n}, {}, false,
                                     [&](bufs &b, bool lib) { CALL(Rgemmtr, (u, t, t2, n, k, ra, b.r[0].data(), lda, b.r[1].data(), ldb, rb, b.r[2].data(), ldc)); });
                        run_case("Cgemmtr", w2, {}, {lda * (nt ? k : n), ldb * (nt2 ? n : k), ldc * n}, false,
                                 [&](bufs &b, bool lib) { CALL(Cgemmtr, (u, t, t2, n, k, ca, b.c[0].data(), lda, b.c[1].data(), ldb, cb, b.c[2].data(), ldc)); });
                        if (u[0] == 'U') {
                            const long lda3 = (nt ? m : k) + 1, ldc3 = m + 1;
                            run_case("Cgemm", tag({t, t2}, {m, n, k, ak, bk}), {}, {lda3 * (nt ? k : m), ldb * (nt2 ? n : k), ldc3 * n}, false,
                                     [&](bufs &b, bool lib) { CALL(Cgemm, (t, t2, m, n, k, ca, b.c[0].data(), lda3, b.c[1].data(), ldb, cb, b.c[2].data(), ldc3)); });
                        }
                    }
                }
            for (const char *s : SD)
                for (const char *u : UL) {
                    const long na = s[0] == 'L' ? m : n, lda = na + 1, ldb = m + 2;
                    const std::string w = tag({s, u}, {m, n, ak, bk});
                    run_case("Rsymm", w, {lda * na, ldb * n, ldb * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rsymm, (s, u, m, n, ra, b.r[0].data(), lda, b.r[1].data(), ldb, rb, b.r[2].data(), ldb)); });
                    run_case("Csymm", w, {}, {lda * na, ldb * n, ldb * n}, false, [&](bufs &b, bool lib) { CALL(Csymm, (s, u, m, n, ca, b.c[0].data(), lda, b.c[1].data(), ldb, cb, b.c[2].data(), ldb)); });
                    run_case("Chemm", w, {}, {lda * na, ldb * n, ldb * n}, false, [&](bufs &b, bool lib) { CALL(Chemm, (s, u, m, n, ca, b.c[0].data(), lda, b.c[1].data(), ldb, cb, b.c[2].data(), ldb)); });
                    if (bk != 0)
                        continue;
                    for (const char *t : TR)
                        for (const char *d : DG) {
                            const std::string w3 = tag({s, u, t, d}, {m, n, ak});
                            run_case("Rtrmm", w3, {lda * na, ldb * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rtrmm, (s, u, t, d, m, n, ra, b.r[0].data(), lda, b.r[1].data(), ldb)); });
                            run_case("Ctrmm", w3, {}, {lda * na, ldb * n}, false, [&](bufs &b, bool lib) { CALL(Ctrmm, (s, u, t, d, m, n, ca, b.c[0].data(), lda, b.c[1].data(), ldb)); });
                            // solves: a triangle of magnitudes in [1/2, 2] and small multipliers
                            auto tri = [&](bufs &b, bool lib, bool cplx) {
                                (void)lib;
                                for (long j = 0; j < na; j++)
                                    for (long i = 0; i < na; i++)
                                        if (i != j) {
                                            if (cplx)
                                                b.c[0][i + j * lda] = b.c[0][i + j * lda] * R(1.0 / (4 * na));
                                            else
                                                b.r[0][i + j * lda] = b.r[0][i + j * lda] / R(4.0 * na);
                                        }
                            };
                            run_case("Rtrsm", w3, {lda * na, ldb * n}, {}, false, [&](bufs &b, bool lib) {
                                tri(b, lib, false);
                                CALL(Rtrsm, (s, u, t, d, m, n, ra, b.r[0].data(), lda, b.r[1].data(), ldb));
                            });
                            run_case("Ctrsm", w3, {}, {lda * na, ldb * n}, false, [&](bufs &b, bool lib) {
                                tri(b, lib, true);
                                CALL(Ctrsm, (s, u, t, d, m, n, ca, b.c[0].data(), lda, b.c[1].data(), ldb));
                            }, d[0] == 'U');
                        }
                }
        }
}

void level2(long m, long n, long kl, long ku)
{
    const long incs[] = {1, 2, -3};
    for (int ak = 0; ak < 4; ak++)
        for (int bk = 0; bk < 4; bk++) {
            if ((ak == 1 || ak == 2) && bk == 2)
                continue;
            const R ra = rs(ak), rb = rs(bk);
            const C ca = cs(ak), cb = cs(bk);
            for (long ix : incs)
                for (long iy : incs) {
                    const long lx = std::labs(ix) * std::max(m, n) + 1, ly = std::labs(iy) * std::max(m, n) + 1, lda = m + 1;
                    for (const char *t : TR) {
                        const std::string w = tag({t}, {m, n, ix, iy, ak, bk});
                        if (t[0] != 'C') {
                            run_case("Rgemv", w, {lda * n, lx, ly}, {}, false, [&](bufs &b, bool lib) { CALL(Rgemv, (t, m, n, ra, b.r[0].data(), lda, b.r[1].data(), ix, rb, b.r[2].data(), iy)); });
                            run_case("Rgbmv", w, {(kl + ku + 2) * n, lx, ly}, {}, false,
                                     [&](bufs &b, bool lib) { CALL(Rgbmv, (t, m, n, kl, ku, ra, b.r[0].data(), kl + ku + 2, b.r[1].data(), ix, rb, b.r[2].data(), iy)); });
                        }
                        run_case("Cgemv", w, {}, {lda * n, lx, ly}, false, [&](bufs &b, bool lib) { CALL(Cgemv, (t, m, n, ca, b.c[0].data(), lda, b.c[1].data(), ix, cb, b.c[2].data(), iy)); });
                        run_case("Cgbmv", w, {}, {(kl + ku + 2) * n, lx, ly}, false,
                                 [&](bufs &b, bool lib) { CALL(Cgbmv, (t, m, n, kl, ku, ca, b.c[0].data(), kl + ku + 2, b.c[1].data(), ix, cb, b.c[2].data(), iy)); });
                    }
                    const std::string w = tag({}, {m, n, ix, iy, ak});
                    if (bk == 0) {
                        run_case("Rger", w, {lx, ly, lda * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rger, (m, n, ra, b.r[0].data(), ix, b.r[1].data(), iy, b.r[2].data(), lda)); });
                        run_case("Cgeru", w, {}, {lx, ly, lda * n}, false, [&](bufs &b, bool lib) { CALL(Cgeru, (m, n, ca, b.c[0].data(), ix, b.c[1].data(), iy, b.c[2].data(), lda)); });
                        run_case("Cgerc", w, {}, {lx, ly, lda * n}, false, [&](bufs &b, bool lib) { CALL(Cgerc, (m, n, ca, b.c[0].data(), ix, b.c[1].data(), iy, b.c[2].data(), lda)); });
                    }
                    for (const char *u : UL) {
                        const long ldn = n + 1, np = n * (n + 1) / 2, k = kl;
                        const std::string wu = tag({u}, {n, k, ix, iy, ak, bk});
                        run_case("Rsymv", wu, {ldn * n, lx, ly}, {}, false, [&](bufs &b, bool lib) { CALL(Rsymv, (u, n, ra, b.r[0].data(), ldn, b.r[1].data(), ix, rb, b.r[2].data(), iy)); });
                        run_case("Rsbmv", wu, {(k + 2) * n, lx, ly}, {}, false, [&](bufs &b, bool lib) { CALL(Rsbmv, (u, n, k, ra, b.r[0].data(), k + 2, b.r[1].data(), ix, rb, b.r[2].data(), iy)); });
                        run_case("Rspmv", wu, {np, lx, ly}, {}, false, [&](bufs &b, bool lib) { CALL(Rspmv, (u, n, ra, b.r[0].data(), b.r[1].data(), ix, rb, b.r[2].data(), iy)); });
                        run_case("Chemv", wu, {}, {ldn * n, lx, ly}, false, [&](bufs &b, bool lib) { CALL(Chemv, (u, n, ca, b.c[0].data(), ldn, b.c[1].data(), ix, cb, b.c[2].data(), iy)); });
                        run_case("Chbmv", wu, {}, {(k + 2) * n, lx, ly}, false, [&](bufs &b, bool lib) { CALL(Chbmv, (u, n, k, ca, b.c[0].data(), k + 2, b.c[1].data(), ix, cb, b.c[2].data(), iy)); });
                        run_case("Chpmv", wu, {}, {np, lx, ly}, false, [&](bufs &b, bool lib) { CALL(Chpmv, (u, n, ca, b.c[0].data(), b.c[1].data(), ix, cb, b.c[2].data(), iy)); });
                        if (bk != 0)
                            continue;
                        run_case("Rsyr", wu, {lx, ldn * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rsyr, (u, n, ra, b.r[0].data(), ix, b.r[1].data(), ldn)); });
                        run_case("Rsyr2", wu, {lx, ly, ldn * n}, {}, false, [&](bufs &b, bool lib) { CALL(Rsyr2, (u, n, ra, b.r[0].data(), ix, b.r[1].data(), iy, b.r[2].data(), ldn)); });
                        run_case("Rspr", wu, {lx, np}, {}, false, [&](bufs &b, bool lib) { CALL(Rspr, (u, n, ra, b.r[0].data(), ix, b.r[1].data())); });
                        run_case("Rspr2", wu, {lx, ly, np}, {}, false, [&](bufs &b, bool lib) { CALL(Rspr2, (u, n, ra, b.r[0].data(), ix, b.r[1].data(), iy, b.r[2].data())); });
                        run_case("Cher", wu, {}, {lx, ldn * n}, false, [&](bufs &b, bool lib) { CALL(Cher, (u, n, ra, b.c[0].data(), ix, b.c[1].data(), ldn)); });
                        run_case("Cher2", wu, {}, {lx, ly, ldn * n}, false, [&](bufs &b, bool lib) { CALL(Cher2, (u, n, ca, b.c[0].data(), ix, b.c[1].data(), iy, b.c[2].data(), ldn)); });
                        run_case("Chpr", wu, {}, {lx, np}, false, [&](bufs &b, bool lib) { CALL(Chpr, (u, n, ra, b.c[0].data(), ix, b.c[1].data())); });
                        run_case("Chpr2", wu, {}, {lx, ly, np}, false, [&](bufs &b, bool lib) { CALL(Chpr2, (u, n, ca, b.c[0].data(), ix, b.c[1].data(), iy, b.c[2].data())); });
                        if (ak != 3 || iy != 1)
                            continue;
                        for (const char *t : TR)
                            for (const char *d : DG) {
                                const std::string wt = tag({u, t, d}, {n, k, ix});
                                // triangles of magnitudes in [1/2, 2] with small off-diagonal elements
                                auto shrink = [&](bufs &b, int kind) {
                                    const long rows = kind == 0 ? ldn : kind == 1 ? k + 2 : 0;
                                    if (kind == 2) {
                                        for (long q = 0; q < np; q++) {
                                            b.r[0].size() ? (void)(b.r[0][q] = b.r[0][q] / R(4.0 * n)) : (void)(b.c[0][q] = b.c[0][q] * R(1.0 / (4 * n)));
                                        }
                                        return;
                                    }
                                    for (long q = 0; q < rows * n; q++)
                                        if (!(kind == 0 && q % rows == q / rows) && !(kind == 1 && q % rows == (u[0] == 'U' ? k : 0))) {
                                            if (b.r[0].size())
                                                b.r[0][q] = b.r[0][q] / R(4.0 * n);
                                            else
                                                b.c[0][q] = b.c[0][q] * R(1.0 / (4 * n));
                                        }
                                };
                                run_case("Rtrmv", wt, {ldn * n, lx}, {}, false, [&](bufs &b, bool lib) { CALL(Rtrmv, (u, t, d, n, b.r[0].data(), ldn, b.r[1].data(), ix)); });
                                run_case("Rtbmv", wt, {(k + 2) * n, lx}, {}, false, [&](bufs &b, bool lib) { CALL(Rtbmv, (u, t, d, n, k, b.r[0].data(), k + 2, b.r[1].data(), ix)); });
                                run_case("Rtpmv", wt, {np, lx}, {}, false, [&](bufs &b, bool lib) { CALL(Rtpmv, (u, t, d, n, b.r[0].data(), b.r[1].data(), ix)); });
                                run_case("Ctrmv", wt, {}, {ldn * n, lx}, false, [&](bufs &b, bool lib) { CALL(Ctrmv, (u, t, d, n, b.c[0].data(), ldn, b.c[1].data(), ix)); });
                                run_case("Ctbmv", wt, {}, {(k + 2) * n, lx}, false, [&](bufs &b, bool lib) { CALL(Ctbmv, (u, t, d, n, k, b.c[0].data(), k + 2, b.c[1].data(), ix)); });
                                run_case("Ctpmv", wt, {}, {np, lx}, false, [&](bufs &b, bool lib) { CALL(Ctpmv, (u, t, d, n, b.c[0].data(), b.c[1].data(), ix)); });
                                if (t[0] != 'C') {
                                    run_case("Rtrsv", wt, {ldn * n, lx}, {}, false, [&](bufs &b, bool lib) { shrink(b, 0); CALL(Rtrsv, (u, t, d, n, b.r[0].data(), ldn, b.r[1].data(), ix)); });
                                    run_case("Rtbsv", wt, {(k + 2) * n, lx}, {}, false, [&](bufs &b, bool lib) { shrink(b, 1); CALL(Rtbsv, (u, t, d, n, k, b.r[0].data(), k + 2, b.r[1].data(), ix)); });
                                    run_case("Rtpsv", wt, {np, lx}, {}, false, [&](bufs &b, bool lib) { shrink(b, 2); CALL(Rtpsv, (u, t, d, n, b.r[0].data(), b.r[1].data(), ix)); });
                                }
                                const bool unit = d[0] == 'U';
                                run_case("Ctrsv", wt, {}, {ldn * n, lx}, false, [&](bufs &b, bool lib) { shrink(b, 0); CALL(Ctrsv, (u, t, d, n, b.c[0].data(), ldn, b.c[1].data(), ix)); }, unit);
                                run_case("Ctbsv", wt, {}, {(k + 2) * n, lx}, false, [&](bufs &b, bool lib) { shrink(b, 1); CALL(Ctbsv, (u, t, d, n, k, b.c[0].data(), k + 2, b.c[1].data(), ix)); }, unit);
                                run_case("Ctpsv", wt, {}, {np, lx}, false, [&](bufs &b, bool lib) { shrink(b, 2); CALL(Ctpsv, (u, t, d, n, b.c[0].data(), b.c[1].data(), ix)); }, unit);
                            }
                    }
                }
        }
}

void level1(long n)
{
    const long incs[] = {1, 2, -3};
    for (int ak = 0; ak < 4; ak++)
        for (long ix : incs)
            for (long iy : incs) {
                const R ra = rs(ak), rb = rs((ak + 1) % 4);
                const C ca = cs(ak);
                const long lx = std::labs(ix) * n + 1, ly = std::labs(iy) * n + 1;
                const std::string w = tag({}, {n, ix, iy, ak});
                run_case("Raxpy", w, {lx, ly}, {}, false, [&](bufs &b, bool lib) {
                    if (lib)
                        Raxpy(n, ra, b.r[0].data(), ix, b.r[1].data(), iy);
                    else
                        Raxpy_ref(n, ra, b.r[0].data(), ix, b.r[1].data(), iy);
                }, ak != 0 && ak != 2);
                run_case("Rrot", w, {lx, ly}, {}, false, [&](bufs &b, bool lib) { CALL(Rrot, (n, b.r[0].data(), ix, b.r[1].data(), iy, ra, rb)); });
                run_case("Rdot", w, {lx, ly}, {}, false, [&](bufs &b, bool lib) {
                    if (lib)
                        b.rres = Rdot(n, b.r[0].data(), ix, b.r[1].data(), iy);
                    else
                        b.rres = Rdot_reference_cpu(n, b.r[0].data(), ix, b.r[1].data(), iy);
                });
                run_case("Caxpy", w, {}, {lx, ly}, false, [&](bufs &b, bool lib) { CALL(Caxpy, (n, ca, b.c[0].data(), ix, b.c[1].data(), iy)); }, ak != 0 && ak != 2);
                run_case("CRrot", w, {}, {lx, ly}, false, [&](bufs &b, bool lib) { CALL(CRrot, (n, b.c[0].data(), ix, b.c[1].data(), iy, ra, rb)); });
                run_case("Cdotc", w, {}, {lx, ly}, false, [&](bufs &b, bool lib) {
                    if (lib)
                        b.cres = Cdotc(n, b.c[0].data(), ix, b.c[1].data(), iy);
                    else
                        b.cres = Cdotc_reference_cpu(n, b.c[0].data(), ix, b.c[1].data(), iy);
                });
                run_case("Cdotu", w, {}, {lx, ly}, false, [&](bufs &b, bool lib) {
                    if (lib)
                        b.cres = Cdotu(n, b.c[0].data(), ix, b.c[1].data(), iy);
                    else
                        b.cres = Cdotu_reference_cpu(n, b.c[0].data(), ix, b.c[1].data(), iy);
                });
                if (iy != 1)
                    continue;
                const bool pos = ix > 0;
                run_case("Rscal", w, {lx}, {}, false, [&](bufs &b, bool lib) { CALL(Rscal, (n, ra, b.r[0].data(), ix)); }, pos && ak != 1);
                run_case("Cscal", w, {}, {lx}, false, [&](bufs &b, bool lib) { CALL(Cscal, (n, ca, b.c[0].data(), ix)); }, pos);
                run_case("CRscal", w, {}, {lx}, false, [&](bufs &b, bool lib) { CALL(CRscal, (n, ra, b.c[0].data(), ix)); }, pos);
                run_case("Rasum", w, {lx}, {}, false, [&](bufs &b, bool lib) {
                    if (lib)
                        b.rres = Rasum(n, b.r[0].data(), ix);
                    else
                        b.rres = Rasum_reference_cpu(n, b.r[0].data(), ix);
                }, pos);
                run_case("RCasum", w, {}, {lx}, false, [&](bufs &b, bool lib) {
                    if (lib)
                        b.rres = RCasum(n, b.c[0].data(), ix);
                    else
                        b.rres = RCasum_reference_cpu(n, b.c[0].data(), ix);
                }, pos);
            }
}

// calls the kernels must decline; the result is still the MPFR one
void declined()
{
    const long n = 24, ld = n + 1;
    auto expect_declined = [&](const char *what, const std::function<void(std::vector<R> &, std::vector<R> &, bool)> &f, bool compare = true) {
        std::vector<R> A((size_t)(ld * n)), B((size_t)(ld * n));
        for (auto &x : A)
            x = rv();
        for (auto &x : B)
            x = rv();
        std::vector<R> A2 = A, B2 = B;
        const unsigned long d0 = mplapack_freal_done.load();
        f(A, B, true);
        const bool u = mplapack_freal_done.load() != d0;
        f(A2, B2, false);
        checks++;
        bool ok = !u;
        for (size_t i = 0; i < B.size() && ok && compare; i++)
            ok = same(B[i], B2[i]);
        if (!ok) {
            failures++;
            std::printf("FAIL declined %s (%s)\n", what, u ? "kernels ran" : "result differs");
        }
    };
    auto syrk = [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        const R al = R(0.75), be = R(-1.25);
        if (lib)
            Rsyrk("U", "N", n, n, al, A.data(), ld, be, B.data(), ld);
        else
            Rsyrk_reference_cpu("U", "N", n, n, al, A.data(), ld, be, B.data(), ld);
    };
    expect_declined("NaN", [&](std::vector<R> &A, std::vector<R> &B, bool lib) { mpfr_set_nan(A[3].get_mpfr_t()); syrk(A, B, lib); });
    expect_declined("Inf", [&](std::vector<R> &A, std::vector<R> &B, bool lib) { mpfr_set_inf(A[5].get_mpfr_t(), -1); syrk(A, B, lib); });
    expect_declined("precision", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        mpfr_prec_round(A[7].get_mpfr_t(), prec + 64, MPFR_RNDN);
        syrk(A, B, lib);
    });
    expect_declined("overflow", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        mpfr_mul_2si(A[0].get_mpfr_t(), A[0].get_mpfr_t(), mpfr_get_emax() - 2, MPFR_RNDN);
        mpfr_set_d(A[0].get_mpfr_t(), 1.0, MPFR_RNDN);
        mpfr_mul_2si(A[0].get_mpfr_t(), A[0].get_mpfr_t(), mpfr_get_emax() - 1, MPFR_RNDN);
        syrk(A, B, lib);
    });
    expect_declined("underflow", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        mpfr_set_d(A[0].get_mpfr_t(), 1.0, MPFR_RNDN);
        mpfr_mul_2si(A[0].get_mpfr_t(), A[0].get_mpfr_t(), mpfr_get_emin() + 1, MPFR_RNDN);
        syrk(A, B, lib);
    });
    expect_declined("solve overflow", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        // B(:,j) grows by a factor 2^(emax/8) per step of the solve
        R big = R(1.0);
        mpfr_mul_2si(big.get_mpfr_t(), big.get_mpfr_t(), mpfr_get_emax() / 8, MPFR_RNDN);
        for (long j = 0; j < n; j++)
            for (long i = 0; i < n; i++)
                A[i + j * ld] = i == j ? R(1.0) : i > j ? -big : R(0.0);
        if (lib)
            Rtrsm("L", "L", "N", "U", n, n, R(1.0), A.data(), ld, B.data(), ld);
        else
            Rtrsm_reference_cpu("L", "L", "N", "U", n, n, R(1.0), A.data(), ld, B.data(), ld);
    });
    expect_declined("division by zero", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        A[4 + 4 * ld] = R(0.0);
        if (lib)
            Rtrsm("L", "U", "N", "N", n, n, R(1.0), A.data(), ld, B.data(), ld);
        else
            Rtrsm_reference_cpu("L", "U", "N", "N", n, n, R(1.0), A.data(), ld, B.data(), ld);
    });
    expect_declined("RNDZ", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        mpfr_set_default_rounding_mode(MPFR_RNDZ);
        syrk(A, B, lib);
        mpfr_set_default_rounding_mode(MPFR_RNDN);
    });
    expect_declined("overlap", [&](std::vector<R> &A, std::vector<R> &B, bool lib) {
        (void)B;
        // C overlaps A (invalid for the parallel MPFR code too: no comparison)
        if (lib)
            Rsyrk("U", "N", n - 2, n - 2, R(0.5), A.data(), ld, R(2.0), A.data() + 2 * ld, ld);
        else
            Rsyrk_reference_cpu("U", "N", n - 2, n - 2, R(0.5), A.data(), ld, R(2.0), A.data() + 2 * ld, ld);
    }, false);
}

// operands that are blocks of one matrix, as in LAPACK: the kernels run
void blocks()
{
    const long n = 40, ld = n + 3, nb = 12;
    auto expect_used = [&](const char *what, const std::function<void(std::vector<R> &, bool)> &f) {
        std::vector<R> A((size_t)(ld * n));
        for (auto &x : A)
            x = rv();
        for (long i = 0; i < n; i++)
            A[i + i * ld] = R(2.0) + abs(A[i + i * ld]);
        std::vector<R> A2 = A;
        const unsigned long d0 = mplapack_freal_done.load();
        f(A, true);
        const bool u = mplapack_freal_done.load() != d0;
        f(A2, false);
        checks++;
        bool ok = u;
        for (size_t i = 0; i < A.size() && ok; i++)
            ok = same(A[i], A2[i]);
        if (!ok) {
            failures++;
            std::printf("FAIL blocks %s (%s)\n", what, u ? "result differs" : "kernels did not run");
        }
    };
    // Rgetrf / Rpotrf steps: A11 \ A12, A22 - A21 A21^T, y = A22 x with x a column
    expect_used("Rtrsm", [&](std::vector<R> &A, bool lib) {
        if (lib)
            Rtrsm("L", "L", "N", "U", nb, n - nb, R(1.0), A.data(), ld, A.data() + nb * ld, ld);
        else
            Rtrsm_reference_cpu("L", "L", "N", "U", nb, n - nb, R(1.0), A.data(), ld, A.data() + nb * ld, ld);
    });
    expect_used("Rsyrk", [&](std::vector<R> &A, bool lib) {
        if (lib)
            Rsyrk("L", "N", n - nb, nb, R(-1.0), A.data() + nb, ld, R(1.0), A.data() + nb + nb * ld, ld);
        else
            Rsyrk_reference_cpu("L", "N", n - nb, nb, R(-1.0), A.data() + nb, ld, R(1.0), A.data() + nb + nb * ld, ld);
    });
    if (prec > 512) // Level 2: up to MPLAPACK_MPFR_CPU_FREAL_L12_PREC
        return;
    expect_used("Rgemv", [&](std::vector<R> &A, bool lib) {
        if (lib)
            Rgemv("T", n - nb, n - nb - 1, R(1.0), A.data() + nb + (nb + 1) * ld, ld, A.data() + nb + nb * ld, 1, R(0.0), A.data() + nb - 1 + (nb + 1) * ld, ld);
        else
            Rgemv_reference_cpu("T", n - nb, n - nb - 1, R(1.0), A.data() + nb + (nb + 1) * ld, ld, A.data() + nb + nb * ld, 1, R(0.0), A.data() + nb - 1 + (nb + 1) * ld, ld);
    });
}

} // namespace

int main()
{
    gmp_randinit_default(gstate);
    gmp_randseed_ui(gstate, 20261123);
    setenv("MPLAPACK_MPFR_CPU_FREAL_MIN", "1", 1);
    setenv("MPLAPACK_BLAS_OMP_MIN", "1", 1);
    const long precs[] = {512, 128, 192, 448, 1024, 2048};
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#endif
    for (long p : precs) {
        prec = p;
        mpfrxx::set_default_precision_bits(p);
        mpfrxx::set_default_mpc_precision_bits(p);
        for (int omp = 0; omp < 2; omp++) {
#ifdef _OPENMP
            omp_set_num_threads(omp ? threads : 1);
#endif
            const bool big = p == 512 && omp;
            level3(7, 5, 3);
            level3(9, 12, 6);
            if (big)
                level3(33, 31, 17);
            level2(9, 7, 2, 1);
            if (big)
                level2(40, 33, 3, 5);
            level1(25);
            if (big)
                level1(301);
            declined();
            blocks();
        }
    }
    long notused = 0;
    for (auto &c : calls)
        if (used[c.first] == 0) {
            std::printf("FAIL %s: the kernels never ran\n", c.first.c_str());
            notused++;
        }
    std::printf("Rfreal_blas: %ld checks, %ld failures, %ld routines not run\n", checks, failures, notused);
    return failures || notused ? 1 : 0;
}
