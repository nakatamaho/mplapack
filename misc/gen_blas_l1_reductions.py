#!/usr/bin/env python3
"""Generate the parallel versions of the reference Level 1 reductions.

The reference reductions are sequential sums (or comparison chains), and a
reordered sum is not the reference result.  What is independent is the
work on each element: the products of Cdotu/Cdotc, RCabs1 of RCasum, the
six-element partial sums of the unrolled Rasum loop, the class and term of
every element in Rnrm2/RCnrm2, |x(i)| (RCabs1) in iRamax/iCamax.
mpblas/optimized/common/blas/<R>_omp.cpp computes those in parallel
(Rblas_split_common.h) and then performs the reference accumulation in the
reference order on them, so the result is the reference one bit for bit.
Unit increments only; <R>.cpp is the dispatcher.  The accumulation code of
Rnrm2/RCnrm2 is the reference text with its loop rewritten (rerun when the
reference changes).

usage: misc/gen_blas_l1_reductions.py   (from the top source directory)
"""
import os
import re

OUT = "mpblas/optimized/common/blas"
LIC = open(os.path.join(OUT, "Cher.cpp")).read().split("*/", 1)[0] + "*/\n"


def head(r, desc):
    return LIC + """
// The reference %(r)s (mpblas/reference/%(r)s.cpp, unchanged) under the name
// %(r)s_reference_cpu, and %(r)s_omp: %(desc)s
// (misc/gen_blas_l1_reductions.py).  Unit increment only.

#define %(r)s %(r)s_reference_cpu
#include "../../../reference/%(r)s.cpp"
#undef %(r)s

#include <vector>
#include "../Rblas_split_common.h"
""" % dict(r=r, desc=desc)


def dispatcher(r, rtype, sig, args, cond):
    return LIC + """
// %(r)s of the optimized libraries: large calls with unit increment compute
// the independent per-element work in parallel (%(r)s_omp.cpp), with the
// reference result bit for bit; the others run the reference code.

#include <mpblas.h>

bool %(r)s_omp(%(sig)s, %(rtype)s &result);
%(rtype)s %(r)s_reference_cpu(%(sig)s);

%(rtype)s %(r)s(%(sig)s)
{
    %(rtype)s result%(init)s;
    if (%(cond)s%(r)s_omp(%(args)s, result))
        return result;
    return %(r)s_reference_cpu(%(args)s);
}
""" % dict(r=r, rtype=rtype, sig=sig, args=args, cond=cond, init=" = 0" if rtype == "INTEGER" else "")


files = {}

# ---- Cdotu, Cdotc: ztemp += (product), the products in parallel ----
for r, prod in [("Cdotu", "zx[i] * zy[i]"), ("Cdotc", "conj(zx[i]) * zy[i]")]:
    sig = "INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy"
    files[r + "_omp.cpp"] = head(r, "the products in parallel, then the sum in the\n// reference order") + """
bool %(r)s_omp(%(sig)s, COMPLEX &result)
{
    using namespace mplapack_blas;
    if (n <= 0 || incx != 1 || incy != 1)
        return false;
    std::vector<COMPLEX> p((size_t)n);
    if (!blas_split(8.0 * n, n, [&](long i0, long i1) {
            for (long i = i0; i < i1; i++)
                p[(size_t)i] = %(prod)s;
        }))
        return false;
    COMPLEX ztemp = COMPLEX(0.0, 0.0);
    for (INTEGER i = 1; i <= n; i = i + 1) {
        ztemp += p[(size_t)(i - 1)];
    }
    result = ztemp;
    return true;
}
""" % dict(r=r, sig=sig, prod=prod)
    files[r + ".cpp"] = dispatcher(r, "COMPLEX", sig, "n, zx, incx, zy, incy", "")

# ---- RCasum: stemp += RCabs1(zx(i)) ----
sig = "INTEGER const n, COMPLEX *zx, INTEGER const incx"
files["RCasum_omp.cpp"] = head("RCasum", "RCabs1 of the elements in parallel, then the sum in\n// the reference order") + """
bool RCasum_omp(%s, REAL &result)
{
    using namespace mplapack_blas;
    if (n <= 0 || incx != 1)
        return false;
    std::vector<REAL> p((size_t)n);
    if (!blas_split(2.0 * n, n, [&](long i0, long i1) {
            for (long i = i0; i < i1; i++)
                p[(size_t)i] = RCabs1(zx[i]);
        }))
        return false;
    REAL stemp = 0.0;
    for (INTEGER i = 1; i <= n; i = i + 1) {
        stemp += p[(size_t)(i - 1)];
    }
    result = stemp;
    return true;
}
""" % sig
files["RCasum.cpp"] = dispatcher("RCasum", "REAL", sig, "n, zx, incx", "")

# ---- Rasum: the six-element sums of the unrolled loop in parallel ----
sig = "INTEGER const n, REAL *dx, INTEGER const incx"
files["Rasum_omp.cpp"] = head("Rasum", "the sums of six of the unrolled reference loop in\n// parallel, then their sum in the reference order") + """
bool Rasum_omp(%s, REAL &result)
{
    using namespace mplapack_blas;
    if (n < 12 || incx != 1)
        return false;
    const INTEGER m = n %% 6, groups = (n - m) / 6;
    std::vector<REAL> g((size_t)groups);
    if (!blas_split(2.0 * n, groups, [&](long q0, long q1) {
            for (long q = q0; q < q1; q++) {
                const INTEGER i = m + 1 + 6 * q;
                g[(size_t)q] = abs(dx[i - 1]) + abs(dx[(i + 1) - 1]) + abs(dx[(i + 2) - 1]) + abs(dx[(i + 3) - 1]) + abs(dx[(i + 4) - 1]) + abs(dx[(i + 5) - 1]);
            }
        }))
        return false;
    REAL dtemp = 0.0;
    for (INTEGER i = 1; i <= m; i = i + 1) {
        dtemp += abs(dx[i - 1]);
    }
    for (INTEGER q = 0; q < groups; q++) {
        dtemp += g[(size_t)q];
    }
    result = dtemp;
    return true;
}
""" % sig
files["Rasum.cpp"] = dispatcher("Rasum", "REAL", sig, "n, dx, incx", "")

# ---- iRamax, iCamax: the first index of the largest value ----
for r, T, val in [("iRamax", "REAL", "abs(dx[i])"), ("iCamax", "COMPLEX", "RCabs1(zx[i])")]:
    xv = "dx" if r == "iRamax" else "zx"
    sig = "INTEGER const n, %s *%s, INTEGER const incx" % (T, xv)
    files[r + "_omp.cpp"] = head(r, "the values in parallel; the first index of their\n// maximum, as the reference finds it (comparisons are exact; values that\n// compare false with everything, NaN, go to the reference)") + """
bool %(r)s_omp(%(sig)s, INTEGER &result)
{
    using namespace mplapack_blas;
    if (n < 2 || incx != 1)
        return false;
    std::vector<REAL> v((size_t)n);
    std::vector<char> ok((size_t)n);
    if (!blas_split(2.0 * n, n, [&](long i0, long i1) {
            for (long i = i0; i < i1; i++) {
                v[(size_t)i] = %(val)s;
                ok[(size_t)i] = v[(size_t)i] == v[(size_t)i];
            }
        }))
        return false;
    for (INTEGER i = 0; i < n; i++)
        if (!ok[(size_t)i])
            return false;
    // the reference: dmax = v(1); for i = 2..n: if (v(i) > dmax) take i
    INTEGER best = 1;
    REAL dmax = v[0];
    for (INTEGER i = 2; i <= n; i = i + 1) {
        if (v[(size_t)(i - 1)] > dmax) {
            best = i;
            dmax = v[(size_t)(i - 1)];
        }
    }
    result = best;
    return true;
}
""" % dict(r=r, sig=sig, val=val)
    files[r + ".cpp"] = dispatcher(r, "INTEGER", sig, "n, %s, incx" % xv, "")

# ---- Rnrm2, RCnrm2: class and term of each value in parallel ----
for r, T, xs in [("Rnrm2", "REAL", ["x[k]"]), ("RCnrm2", "COMPLEX", ["x[k].real()", "x[k].imag()"])]:
    src = open("mpblas/reference/%s.cpp" % r).read()
    m = re.search(r"^REAL %s\((.*?)\) \{\n" % r, src, re.M)
    sig = m.group(1)
    end = src.index("\n}\n", m.end())
    body = src[m.end():end + 3]
    loop = re.search(r"    for \(i = 1; i <= n; i = i \+ 1\) \{\n.*?\n    \}\n", body, re.S)
    w = len(xs)
    newloop = """    for (i = 1; i <= n; i = i + 1) {
        for (int c = 0; c < %(w)d; c++) {
            const size_t t = (size_t)((ix - 1) * %(w)d + c);
            if (cls[t] == 0) {
                abig += term[t];
                notbig = false;
            } else if (cls[t] == 1) {
                if (notbig) {
                    asml += term[t];
                }
            } else {
                amed += term[t];
            }
        }
        ix += incx;
    }
""" % dict(w=w)
    # the values in the order of the reference loop: x(k).real(), x(k).imag()
    body = body[:loop.start()] + newloop + body[loop.end():]
    pre = "\n".join("""                ax = abs(%s);
                if (ax > tbig) {
                    cls[t] = 0, term[t] = pow2((ax * sbig));
                } else if (ax < tsml) {
                    cls[t] = 1, term[t] = pow2((ax * ssml));
                } else {
                    cls[t] = 2, term[t] = pow2(ax);
                }
                t++;""" % e for e in xs)
    files[r + "_omp.cpp"] = head(r, "the class (big, small, medium) and the term of every\n// value in parallel, then the reference accumulation in the reference order") + """
namespace {

// The reference %(r)s with its loop over x replaced by the accumulation of
// the precomputed classes and terms (the reference computes the term of a
// small value only while no big value has occurred; the term itself does
// not depend on that).
REAL %(r)s_terms(%(sig)s, const signed char *cls, const REAL *term) {
%(body)s
} // namespace

bool %(r)s_omp(%(sig)s, REAL &result)
{
    using namespace mplapack_blas;
    if (n <= 0 || incx != 1)
        return false;
    std::vector<signed char> cls((size_t)n * %(w)d);
    std::vector<REAL> term((size_t)n * %(w)d);
    const auto &ap = mplapack::get_arithmetic_params<REAL>();
    const auto bp = mplapack::make_blue_scaling_params(ap);
    const REAL tbig = bp.tbig;
    const REAL sbig = bp.sbig;
    const REAL tsml = bp.tsml;
    const REAL ssml = bp.ssml;
    if (!blas_split(4.0 * n * %(w)d, n, [&](long k0, long k1) {
            REAL ax = 0.0;
            for (long k = k0; k < k1; k++) {
                size_t t = (size_t)k * %(w)d;
%(pre)s
            }
        }))
        return false;
    result = %(r)s_terms(n, x, incx, cls.data(), term.data());
    return true;
}
""" % dict(r=r, sig=sig, body=body.rstrip("\n"), w=w, pre=pre)
    files[r + ".cpp"] = dispatcher(r, "REAL", sig, "n, x, incx", "")

for f, s in files.items():
    open(os.path.join(OUT, f), "w").write(s)
    print("wrote", os.path.join(OUT, f))
