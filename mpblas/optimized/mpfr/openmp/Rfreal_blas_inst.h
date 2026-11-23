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

// One precision of the fixed-precision MPFR BLAS (Rfreal_blas.h), included
// by Rfreal_blas_<PB>.cpp with FRB_PB defined: in the namespace frb<PB>,
// REAL and COMPLEX are fr<PB> and fc<PB> (Rfreal_types.h) and the files of
// the library's OpenMP BLAS (../../common/blas/<name>_omp.cpp with the
// reference code) are compiled again for them, each in its namespace
// <name>_ns; <name>_freal() copies the
// operands in, runs <name>_omp of these types (else the reference code),
// checks that nothing was reported and copies the results back.
//
// On x86-64 the file is compiled for BMI2 and ADX (the mulx/adcx/adox
// kernels); Rfreal_blas.cpp calls it only on a CPU that has them.  Every
// header that the included files use is included before that.

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <deque>
#include <vector>
#include <mpblas.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../../common/Rblas_split_common.h"
#include "../../common/Rblas_cpu_common.h"
#include "../../common/Rblas1_fast_common.h"
#define MPLAPACK_GEMM_NO_OZAKI
#include "../../common/Rgemm_cpu_common.h"
#include "Rgemm_freal.h"
#include "Rfreal_blas.h"

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#pragma GCC optimize("O3")
#pragma GCC target("bmi2,adx")
#define MPFREAL_FORCE_X86_MULX 1
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC optimize("O3")
#endif

#define FRB_CAT2(a, b) a##b
#define FRB_CAT(a, b) FRB_CAT2(a, b)
#define FRB_NS FRB_CAT(frb, FRB_PB)
#define MPFREAL_NS FRB_NS
#include "../lgpl/Rfreal_host.h"
#include "../lgpl/Rfround_host.h"
#include "../lgpl/Rfcomplex_host.h"
#include "Rfreal_types.h"

#define MPLAPACK_FREAL_INSTANCE 1

namespace FRB_NS {
typedef fr<FRB_PB> REAL;
typedef fc<FRB_PB> COMPLEX;
using std::abs;
using std::sqrt;
// invalid arguments: the MPFR code runs and reports them
inline void freal_xerbla(const char *, int) { mplapack_freal_fail(); }
} // namespace FRB_NS

#undef Mxerbla
#define Mxerbla FRB_NS::freal_xerbla

namespace FRB_NS {
#include "../../../reference/RCabs1.cpp"
namespace Rsyrk_ns {
#include "../../common/blas/Rsyrk_omp.cpp"
}
namespace Rsyr2k_ns {
#include "../../common/blas/Rsyr2k_omp.cpp"
}
namespace Rsymm_ns {
#include "../../common/blas/Rsymm_omp.cpp"
}
namespace Rtrmm_ns {
#include "../../common/blas/Rtrmm_omp.cpp"
}
namespace Rtrsm_ns {
#include "../../common/blas/Rtrsm_omp.cpp"
}
namespace Rgemmtr_ns {
#include "../../common/blas/Rgemmtr_omp.cpp"
}
namespace Rgemv_ns {
#include "../../common/blas/Rgemv_omp.cpp"
}
namespace Rsymv_ns {
#include "../../common/blas/Rsymv_omp.cpp"
}
namespace Rtrmv_ns {
#include "../../common/blas/Rtrmv_omp.cpp"
}
namespace Rger_ns {
#include "../../common/blas/Rger_omp.cpp"
}
namespace Rsyr_ns {
#include "../../common/blas/Rsyr_omp.cpp"
}
namespace Rsyr2_ns {
#include "../../common/blas/Rsyr2_omp.cpp"
}
namespace Rgbmv_ns {
#include "../../common/blas/Rgbmv_omp.cpp"
}
namespace Rsbmv_ns {
#include "../../common/blas/Rsbmv_omp.cpp"
}
namespace Rspmv_ns {
#include "../../common/blas/Rspmv_omp.cpp"
}
namespace Rspr_ns {
#include "../../common/blas/Rspr_omp.cpp"
}
namespace Rspr2_ns {
#include "../../common/blas/Rspr2_omp.cpp"
}
namespace Rtbmv_ns {
#include "../../common/blas/Rtbmv_omp.cpp"
}
namespace Rtpmv_ns {
#include "../../common/blas/Rtpmv_omp.cpp"
}
namespace Rtrsv_ns {
#include "../../common/blas/Rtrsv_omp.cpp"
}
namespace Rtbsv_ns {
#include "../../common/blas/Rtbsv_omp.cpp"
}
namespace Rtpsv_ns {
#include "../../common/blas/Rtpsv_omp.cpp"
}
namespace Rscal_ns {
#include "../../common/blas/Rscal_omp.cpp"
}
namespace Rrot_ns {
#include "../../common/blas/Rrot_omp.cpp"
}
namespace Rasum_ns {
#include "../../common/blas/Rasum_omp.cpp"
}
namespace Rdot_ns {
#include "../../common/blas/Rdot_par.cpp"
}
namespace Raxpy_ns {
#define Raxpy Raxpy_reference_cpu
#include "../../../reference/Raxpy.cpp"
#undef Raxpy
}
namespace Cgemm_ns {
#include "../../common/blas/Cgemm_omp.cpp"
}
namespace Csymm_ns {
#include "../../common/blas/Csymm_omp.cpp"
}
namespace Chemm_ns {
#include "../../common/blas/Chemm_omp.cpp"
}
namespace Ctrmm_ns {
#include "../../common/blas/Ctrmm_omp.cpp"
}
namespace Ctrsm_ns {
#include "../../common/blas/Ctrsm_omp.cpp"
}
namespace Csyrk_ns {
#include "../../common/blas/Csyrk_omp.cpp"
}
namespace Cherk_ns {
#include "../../common/blas/Cherk_omp.cpp"
}
namespace Csyr2k_ns {
#include "../../common/blas/Csyr2k_omp.cpp"
}
namespace Cher2k_ns {
#include "../../common/blas/Cher2k_omp.cpp"
}
namespace Cgemmtr_ns {
#include "../../common/blas/Cgemmtr_omp.cpp"
}
namespace Cgemv_ns {
#include "../../common/blas/Cgemv_omp.cpp"
}
namespace Cgbmv_ns {
#include "../../common/blas/Cgbmv_omp.cpp"
}
namespace Chemv_ns {
#include "../../common/blas/Chemv_omp.cpp"
}
namespace Chbmv_ns {
#include "../../common/blas/Chbmv_omp.cpp"
}
namespace Chpmv_ns {
#include "../../common/blas/Chpmv_omp.cpp"
}
namespace Cgeru_ns {
#include "../../common/blas/Cgeru_omp.cpp"
}
namespace Cgerc_ns {
#include "../../common/blas/Cgerc_omp.cpp"
}
namespace Cher_ns {
#include "../../common/blas/Cher_omp.cpp"
}
namespace Cher2_ns {
#include "../../common/blas/Cher2_omp.cpp"
}
namespace Chpr_ns {
#include "../../common/blas/Chpr_omp.cpp"
}
namespace Chpr2_ns {
#include "../../common/blas/Chpr2_omp.cpp"
}
namespace Ctrmv_ns {
#include "../../common/blas/Ctrmv_omp.cpp"
}
namespace Ctbmv_ns {
#include "../../common/blas/Ctbmv_omp.cpp"
}
namespace Ctpmv_ns {
#include "../../common/blas/Ctpmv_omp.cpp"
}
namespace Ctrsv_ns {
#include "../../common/blas/Ctrsv_omp.cpp"
}
namespace Ctbsv_ns {
#include "../../common/blas/Ctbsv_omp.cpp"
}
namespace Ctpsv_ns {
#include "../../common/blas/Ctpsv_omp.cpp"
}
namespace Caxpy_ns {
#include "../../common/blas/Caxpy_omp.cpp"
}
namespace Cscal_ns {
#include "../../common/blas/Cscal_omp.cpp"
}
namespace CRscal_ns {
#include "../../common/blas/CRscal_omp.cpp"
}
namespace CRrot_ns {
#include "../../common/blas/CRrot_omp.cpp"
}
namespace Cdotc_ns {
#include "../../common/blas/Cdotc_omp.cpp"
}
namespace Cdotu_ns {
#include "../../common/blas/Cdotu_omp.cpp"
}
namespace RCasum_ns {
#include "../../common/blas/RCasum_omp.cpp"
}
} // namespace FRB_NS

#undef Mxerbla
#define Mxerbla Mxerbla_mpfr

namespace FRB_NS {

// ---- copies of the operands ----

// the referenced part of a matrix (rows x cols): all, a triangle, a band
// (row r of column j is element (r - ku + j, j) of an m-row matrix), without
// the diagonal for diag = U
struct fregion {
    int kind; // 0 all, 1 upper, 2 lower, 3 band, 4 row 0
    long kl, ku, m;
    bool nodiag;
    bool ref(long i, long j) const
    {
        switch (kind) {
        case 1:
            return i < j || (i == j && !nodiag);
        case 2:
            return i > j || (i == j && !nodiag);
        case 3: {
            const long r = i - ku + j;
            return r >= 0 && r < m && !(nodiag && i == ku);
        }
        case 4: // a vector: the first row
            return i == 0;
        default:
            return true;
        }
    }
};
inline fregion rg_all() { return fregion{0, 0, 0, 0, false}; }
inline fregion rg_tri(bool upper, bool nodiag = false) { return fregion{upper ? 1 : 2, 0, 0, 0, nodiag}; }
inline fregion rg_band(long kl, long ku, long m, bool nodiag = false) { return fregion{3, kl, ku, m, nodiag}; }

enum { F_IN = 1, F_OUT = 2, F_INOUT = 3 };

inline bool load(REAL &f, mpfr_srcptr x, long emin, long emax)
{
    if (mpfr_get_prec(x) != FRB_PB)
        return false;
    if (mpfr_zero_p(x)) {
        f.v.set_zero();
        f.v.sign = mpfr_signbit(x) ? -1 : 1;
        return true;
    }
    if (!mpfr_regular_p(x))
        return false;
    const long e = mpfr_get_exp(x);
    if (e > emax || e < emin)
        return false;
    f.v.sign = mpfr_sgn(x) < 0 ? -1 : 1;
    f.v.exp = e;
    std::memcpy(f.v.m, x->_mpfr_d, sizeof f.v.m);
    return true;
}
inline void store(mpfr_ptr x, const REAL &f)
{
    if (f.is_zero()) {
        mpfr_set_zero(x, f.v.sign);
        return;
    }
    std::memcpy(x->_mpfr_d, f.v.m, sizeof f.v.m);
    x->_mpfr_exp = f.v.exp;
    x->_mpfr_sign = f.v.sign;
}
inline bool load(REAL &f, const mpfr_class &x, long emin, long emax) { return load(f, x.get_mpfr_t(), emin, emax); }
inline bool load(COMPLEX &f, const mpc_class &x, long emin, long emax) { return load(f.re, mpc_realref(x.mpc_data()), emin, emax) && load(f.im, mpc_imagref(x.mpc_data()), emin, emax); }
inline void store(mpfr_class &x, const REAL &f) { store(x.get_mpfr_t(), f); }
inline void store(mpc_class &x, const COMPLEX &f)
{
    store(mpc_realref(x.mpc_data()), f.re);
    store(mpc_imagref(x.mpc_data()), f.im);
}
inline bool prec_ok(const mpfr_class &x) { return mpfr_get_prec(x.get_mpfr_t()) == FRB_PB; }
inline bool prec_ok(const mpc_class &x) { return mpfr_get_prec(mpc_realref(x.mpc_data())) == FRB_PB && mpfr_get_prec(mpc_imagref(x.mpc_data())) == FRB_PB; }

template <class S> struct value_of;
template <> struct value_of<mpfr_class> {
    typedef REAL type;
};
template <> struct value_of<mpc_class> {
    typedef COMPLEX type;
};

// the referenced rows of column j: [lo, hi] (empty if lo > hi)
inline void ref_rows(const fregion &rg, long rows, long j, long &lo, long &hi)
{
    lo = 0, hi = rows - 1;
    switch (rg.kind) {
    case 1:
        hi = std::min(rg.nodiag ? j - 1 : j, rows - 1);
        break;
    case 2:
        lo = rg.nodiag ? j + 1 : j;
        break;
    case 3:
        lo = std::max(0L, rg.ku - j), hi = std::min(rows - 1, rg.m - 1 + rg.ku - j);
        break;
    case 4:
        hi = 0;
        break;
    }
}

typedef std::vector<std::pair<const char *, const char *> > byte_ranges;

struct operand_base {
    byte_ranges bytes; // the referenced objects, ascending, merged
    int mode;
    virtual void write_back(bool par) = 0;
    virtual ~operand_base() {}
};

// rows x cols of X (leading dimension ld) in a buffer with leading
// dimension rows; a vector is row 0 of |inc| x n.
template <class S> struct operand : operand_base {
    typedef typename value_of<S>::type V;
    S *X;
    long rows, cols, xld;
    fregion rg;
    std::vector<V> buf;
    V *p;
    long ld;
    void write_back(bool par)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (par)
#endif
        for (long j = 0; j < cols; j++) {
            long lo, hi;
            ref_rows(rg, rows, j, lo, hi);
            for (long i = lo; i <= hi; i++)
                store(X[i + j * xld], buf[i + j * rows]);
        }
    }
};

// the copies in parallel for large operands
inline bool copy_parallel(double elements)
{
#ifdef _OPENMP
    return elements >= 4096 && omp_get_max_threads() >= 2 && !omp_in_parallel();
#else
    (void)elements;
    return false;
#endif
}

inline bool ranges_meet(const byte_ranges &a, const byte_ranges &b)
{
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i].second <= b[j].first)
            i++;
        else if (b[j].second <= a[i].first)
            j++;
        else
            return true;
    }
    return false;
}

struct fcall {
    double work;
    bool ok;
    std::atomic<bool> bad;
    mplapack_freal_state saved;
    long emin, emax;
    std::deque<operand_base *> ops;
    // level 3: always; levels 1 and 2 (one or two operations per element
    // copied): up to MPLAPACK_MPFR_CPU_FREAL_L12_PREC bits (512), and level 1
    // only when it runs on one thread (the OpenMP MPFR code is as fast)
    fcall(double w, int level) : work(w), ok(w >= (double)mplapack_freal_min_work()), bad(false), saved(mplapack_freal_tls)
    {
        if (level < 3 && FRB_PB > mplapack_freal_l12_prec())
            ok = false;
#ifdef _OPENMP
        if (level == 1 && omp_get_max_threads() >= 2 && !omp_in_parallel())
            ok = false;
#endif
        // results must stay in MPFR's exponent range (and far from the
        // limits of long)
        const long lim = 1L << 40;
        emin = std::max((long)mpfr_get_emin(), -lim);
        emax = std::min((long)mpfr_get_emax(), lim);
        mplapack_freal_tls.emin = emin;
        mplapack_freal_tls.emax = emax;
        mplapack_freal_tls.bad = &bad;
    }
    ~fcall()
    {
        mplapack_freal_tls = saved;
        for (size_t i = 0; i < ops.size(); i++)
            delete ops[i];
    }
    REAL re(const mpfr_class &x)
    {
        REAL f;
        ok = ok && load(f, x, emin, emax);
        return f;
    }
    COMPLEX cx(const mpc_class &x)
    {
        COMPLEX f;
        ok = ok && load(f, x, emin, emax);
        return f;
    }
    template <class S> operand<S> &mat(S *X, long rows, long cols, long ld, fregion rg, int mode)
    {
        operand<S> *o = new operand<S>;
        ops.push_back(o);
        o->X = X, o->rows = rows, o->cols = cols, o->xld = ld, o->rg = rg, o->mode = mode;
        o->ld = std::max(rows, 1L);
        o->p = 0;
        if (!ok)
            return *o;
        if (rows < 0 || cols < 0 || ld < std::max(rows, 1L) || X == 0) {
            ok = false;
            return *o;
        }
        o->buf.resize((size_t)std::max(rows * cols, 1L));
        o->p = o->buf.data();
        for (long j = 0; j < cols; j++) {
            long lo, hi;
            ref_rows(rg, rows, j, lo, hi);
            if (lo > hi)
                continue;
            const char *b0 = (const char *)&X[lo + j * ld], *b1 = (const char *)(&X[hi + j * ld] + 1);
            if (!o->bytes.empty() && o->bytes.back().second == b0)
                o->bytes.back().second = b1;
            else
                o->bytes.push_back(std::make_pair(b0, b1));
        }
        long fails = 0;
        const long e0 = emin, e1 = emax;
        typename operand<S>::V *buf = o->p;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : fails) if (copy_parallel((double)rows * cols))
#endif
        for (long j = 0; j < cols; j++) {
            long lo, hi;
            ref_rows(rg, rows, j, lo, hi);
            for (long i = lo; i <= hi; i++) {
                const S &x = X[i + j * ld];
                if (!((mode & F_IN) ? load(buf[i + j * rows], x, e0, e1) : prec_ok(x)))
                    fails++;
            }
        }
        ok = fails == 0;
        return *o;
    }
    // n elements, increment inc, in a buffer with the same increment (the
    // reference code takes other loops, other orders of the sums, for inc = 1)
    template <class S> operand<S> &vec(S *x, long n, long inc, int mode)
    {
        if (inc == 0)
            ok = false;
        const long a = inc == 0 ? 1 : std::labs(inc);
        operand<S> &o = mat(x, a, n, a, fregion{4, 0, 0, 0, false}, mode);
        o.ld = inc;
        return o;
    }
    template <class S> operand<S> &packed(S *ap, long n, int mode) { return vec(ap, n * (n + 1) / 2, 1, mode); }

    // f() runs the routine on the copies; true: the results are copied back
    template <class F> bool run(F f)
    {
        if (!ok)
            return false;
        // an operand written must not share an object with another operand
        // (blocks of one matrix are fine)
        for (size_t i = 0; i < ops.size(); i++)
            for (size_t j = 0; j < ops.size(); j++)
                if (i != j && (ops[i]->mode & F_OUT) && ranges_meet(ops[i]->bytes, ops[j]->bytes))
                    return false;
        const unsigned long orphan = mplapack_freal_orphan.load();
        f();
        if (bad.load() || mplapack_freal_orphan.load() != orphan)
            return false;
        for (size_t i = 0; i < ops.size(); i++)
            if (ops[i]->mode & F_OUT)
                ops[i]->write_back(copy_parallel(work));
        mplapack_freal_done.fetch_add(1, std::memory_order_relaxed);
        return true;
    }
};

inline bool lsm(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline int rw(const REAL &beta) { return beta.is_zero() ? F_OUT : F_INOUT; }
inline int rw(const COMPLEX &beta) { return (beta.re.is_zero() && beta.im.is_zero()) ? F_OUT : F_INOUT; }

} // namespace FRB_NS

#include "Rfreal_blas_drivers.h"
