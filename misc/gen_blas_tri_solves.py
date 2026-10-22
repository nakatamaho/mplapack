#!/usr/bin/env python3
"""Generate the dd/td/qd versions of the triangular solves Rtrsv, Ctrsv,
Rtbsv, Ctbsv, Rtpsv and Ctpsv (mpblas/optimized/common/blas/<R>_omp.cpp).

The columns are solved one after the other, as in the reference.  trans = N:
the updates x(i) = x(i) - temp*a(i,j) of a column are independent and use
the AVX2 loops of Rblas_simd_common.h / Cblas_simd_common.h.  trans = T, C:
the products a(i,j)*x(i) (conj(a(i,j))*x(i)) of a column are computed with
the AVX2 loops first, then subtracted from temp in the reference order.
The operations are those of the reference, so the result is the reference
one bit for bit.  Unit increment only.

usage: misc/gen_blas_tri_solves.py   (from the top source directory)
"""
import os

OUT = "mpblas/optimized/common/blas"
LIC = open(os.path.join(OUT, "Cher.cpp")).read().split("*/", 1)[0] + "*/\n"

for n, T, kind in [("Rtrsv", "R", "full"), ("Ctrsv", "C", "full"), ("Rtbsv", "R", "band"), ("Ctbsv", "C", "band"), ("Rtpsv", "R", "packed"), ("Ctpsv", "C", "packed")]:
    TY = "REAL" if T == "R" else "COMPLEX"
    Z = "0.0" if T == "R" else "COMPLEX(0.0, 0.0)"
    simd = "MPLAPACK_BLAS_SIMD_AVAILABLE" if T == "R" else "MPLAPACK_CBLAS_SIMD_AVAILABLE"
    hdr = "Rblas_simd_common.h" if T == "R" else "Cblas_simd_common.h"
    ax = "bk_axpy" if T == "R" else "cbk_axpy"
    vm = "bk_vmul(cnt, p.data(), col, x + r0)" if T == "R" else "cbk_vmul(cnt, p.data(), col, x + r0, !noconj)"
    if kind == "full":
        sig = "const char *uplo, const char *trans, const char *diag, INTEGER const n, %s *a, INTEGER const lda, %s *x, INTEGER const incx" % (TY, TY)
        chk = "n <= 0 || lda < std::max((INTEGER)1, n) || incx != 1"
        work = "(double)n * n"
        geo = """        const %(TY)s *d = a + j + j * lda; // a(j, j)
        const long r0 = upper ? 0 : j + 1;
        const long cnt = upper ? j : (long)n - j - 1;
        const %(TY)s *col = a + r0 + j * lda; // a(r0, j)""" % dict(TY=TY)
        un = "(void)uplo, (void)trans, (void)diag, (void)n, (void)a, (void)lda, (void)x, (void)incx;"
        store = "column"
    elif kind == "band":
        sig = "const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, %s *a, INTEGER const lda, %s *x, INTEGER const incx" % (TY, TY)
        chk = "n <= 0 || k < 0 || lda < k + 1 || incx != 1"
        work = "(double)n * (k + 1)"
        geo = """        const %(TY)s *d = upper ? a + k + j * lda : a + j * lda; // a(j, j)
        const long r0 = upper ? std::max(0L, j - (long)k) : j + 1;
        const long cnt = upper ? j - r0 : std::min((long)n - 1, j + (long)k) - j;
        const %(TY)s *col = upper ? a + (k + r0 - j) + j * lda : a + 1 + j * lda; // a(r0, j)""" % dict(TY=TY)
        un = "(void)uplo, (void)trans, (void)diag, (void)n, (void)k, (void)a, (void)lda, (void)x, (void)incx;"
        store = "band column"
    else:
        sig = "const char *uplo, const char *trans, const char *diag, INTEGER const n, %s *ap, %s *x, INTEGER const incx" % (TY, TY)
        chk = "n <= 0 || incx != 1"
        work = "(double)n * n / 2"
        geo = """        // column j of the packed triangle: upper from a(0, j), lower from a(j, j)
        const long s = upper ? j * (j + 1) / 2 : j * (long)n - j * (j - 1) / 2;
        const %(TY)s *d = upper ? ap + s + j : ap + s;
        const long r0 = upper ? 0 : j + 1;
        const long cnt = upper ? j : (long)n - j - 1;
        const %(TY)s *col = upper ? ap + s : ap + s + 1;""" % dict(TY=TY)
        un = "(void)uplo, (void)trans, (void)diag, (void)n, (void)ap, (void)x, (void)incx;"
        store = "packed column"
    if T == "R":
        div = "temp = temp / *d;"
        trck = 'notr = Mlsame(trans, "N");\n    if (!notr && !Mlsame(trans, "T") && !Mlsame(trans, "C"))\n        return false;'
        nc = ""
    else:
        div = "temp = noconj ? temp / *d : temp / conj(*d);"
        trck = 'notr = Mlsame(trans, "N");\n    if (!notr && !Mlsame(trans, "T") && !Mlsame(trans, "C"))\n        return false;'
        nc = '    const bool noconj = Mlsame(trans, "T");\n'
    s = LIC + """
// The reference %(n)s (mpblas/reference/%(n)s.cpp, unchanged) under the name
// %(n)s_reference_cpu, and %(n)s_omp: for dd/td/qd and incx = 1 the loops of the
// reference with AVX2 (misc/gen_blas_tri_solves.py): trans = N, the updates
// of x in a %(store)s; trans = T, C, the products of a %(store)s, then
// subtracted in the reference order.  The columns are solved one after the
// other, as in the reference.

#define %(n)s %(n)s_reference_cpu
#include "../../../reference/%(n)s.cpp"
#undef %(n)s

#include <algorithm>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../%(hdr)s"
#endif

bool %(n)s_omp(%(sig)s)
{
#ifdef %(simd)s
    using namespace mplapack_blas;
    const bool upper = Mlsame(uplo, "U"), nounit = Mlsame(diag, "N"), %(trck)s
    if ((!upper && !Mlsame(uplo, "L")) || (!nounit && !Mlsame(diag, "U")) || %(chk)s)
        return false;
    if (blas_env_long("MPLAPACK_BLAS_SIMD", 1) == 0 || %(work)s < (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768))
        return false;
    const %(TY)s zero = %(Z)s;
%(nc)s    if (notr) {
        for (long jj = 0; jj < n; jj++) {
            const long j = upper ? n - 1 - jj : jj;
%(geo)s
            if (x[j] != zero) {
                if (nounit)
                    x[j] = x[j] / *d;
                const %(TY)s temp = x[j];
                %(ax)s(cnt, x + r0, 1, temp, col, 1, true); // x(i) = x(i) - temp*a(i,j)
            }
        }
        return true;
    }
    // x(j) = (x(j) - a(r0, j)*x(r0) - ...) / a(j, j): upper i ascending, lower
    // i descending, as the reference
    std::vector<%(TY)s> p((size_t)n);
    for (long jj = 0; jj < n; jj++) {
        const long j = upper ? jj : n - 1 - jj;
%(geo)s
        %(vm)s;
        %(TY)s temp = x[j];
        if (upper)
            for (long t = 0; t < cnt; t++)
                temp = temp - p[(size_t)t];
        else
            for (long t = cnt - 1; t >= 0; t--)
                temp = temp - p[(size_t)t];
        if (nounit)
            %(div)s
        x[j] = temp;
    }
    return true;
#else
    %(un)s
    return false;
#endif
}
""" % dict(n=n, store=store, hdr=hdr, sig=sig, simd=simd, trck=trck, chk=chk, work=work, TY=TY, Z=Z, nc=nc, geo=geo.replace("\n        ", "\n            ").replace("        const", "            const", 1), ax=ax, vm=vm, div=div, un=un)
    # the second use of geo keeps the 8-space indent
    s = s.replace("        const long j = upper ? jj : n - 1 - jj;\n" + geo.replace("\n        ", "\n            ").replace("        const", "            const", 1),
                  "        const long j = upper ? jj : n - 1 - jj;\n" + geo)
    open(os.path.join(OUT, n + "_omp.cpp"), "w").write(s)
    print("wrote", n)
