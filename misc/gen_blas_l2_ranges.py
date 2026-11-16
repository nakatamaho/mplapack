#!/usr/bin/env python3
"""Generate the range variants of reference Level 2 BLAS routines.

mpblas/optimized/common/blas/<R>_omp.cpp holds <R>_rows, the unit-increment
branch of mpblas/reference/<R>.cpp with its loops restricted to a range of
rows of y (Chemv, Chbmv, Chpmv, Rsbmv, Rspmv, Rgbmv, Cgbmv), of elements of
x (Ctrmv, Rtbmv, Ctbmv, Rtpmv, Ctpmv: read from a copy of the original x)
or of columns of AP (Rspr, Rspr2, Chpr, Chpr2), and <R>_omp, which runs
disjoint ranges in parallel; <R>.cpp is the dispatcher.  The statements are
those of the reference, in the same order for every element, so the result
is the reference one bit for bit.  The text below follows the reference
sources; check it against them when they change.

usage: misc/gen_blas_l2_ranges.py   (from the top source directory)
"""
import os
import shutil
import subprocess
OUT = 'mpblas/optimized/common/blas'
# the license header of the generated files
lic = open(os.path.join(OUT, 'Cher.cpp')).read().split('*/', 1)[0] + '*/\n'

def emit(name, sig, args, desc, check, work, parts, rows_fn, extra_includes=''):
    # sig: parameter list string of the routine; args: argument names list
    omp = lic + '''
// The reference %(n)s (mpblas/reference/%(n)s.cpp, unchanged) under the name
// %(n)s_reference_cpu, and %(n)s_omp: %(desc)s
// Every element receives the operations of the reference code in the same
// order, so the result is the reference one bit for bit.  Unit increments
// only; the other calls run the reference code.

#define %(n)s %(n)s_reference_cpu
#include "../../../reference/%(n)s.cpp"
#undef %(n)s

#include <algorithm>
#include <vector>
#include "../Rblas_split_common.h"

namespace {

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline INTEGER imax(INTEGER a, INTEGER b) { return a > b ? a : b; }
inline INTEGER imin(INTEGER a, INTEGER b) { return a < b ? a : b; }

%(rows)s
} // namespace

bool %(n)s_omp(%(sig)s)
{
    using namespace mplapack_blas;
%(check)s
    return blas_split(%(work)s, %(parts)s, [&](long p0, long p1) { %(n)s_rows(%(args)s, (INTEGER)p0, (INTEGER)p1); });
}
''' % dict(n=name, desc=desc, sig=sig, rows=rows_fn, check=check, work=work, parts=parts, args=', '.join(args))
    open(os.path.join(OUT, name + '_omp.cpp'), 'w').write(omp)
    disp = lic + '''
// %(n)s of the optimized libraries: large calls with unit increments run
// in parallel (%(n)s_omp.cpp), with the reference result bit for bit; the
// others run the reference code.

#include <mpblas.h>

bool %(n)s_omp(%(sig)s);
void %(n)s_reference_cpu(%(sig)s);

void %(n)s(%(sig)s)
{
    if (%(n)s_omp(%(a)s))
        return;
    %(n)s_reference_cpu(%(a)s);
}
''' % dict(n=name, sig=sig, a=', '.join(args))
    open(os.path.join(OUT, name + '.cpp'), 'w').write(disp)

ROWDOC = '// %s restricted to the elements r0 < i <= r1 of %s (1-based).\n'

def beta_block(T, ylen='n'):
    z = 'COMPLEX(0.0, 0.0)' if T == 'COMPLEX' else '0.0'
    o = 'COMPLEX(1.0, 0.0)' if T == 'COMPLEX' else '1.0'
    return '''    const %(T)s zero = %(z)s;
    const %(T)s one = %(o)s;
    INTEGER i = 0;
    if (beta != one) {
        if (beta == zero) {
            for (i = r0 + 1; i <= r1; i = i + 1) {
                y[i - 1] = zero;
            }
        } else {
            for (i = r0 + 1; i <= r1; i = i + 1) {
                y[i - 1] = beta * y[i - 1];
            }
        }
    }
''' % dict(T=T, z=z, o=o)

def mvcheck(T, extra):
    z = 'COMPLEX(0.0, 0.0)' if T == 'COMPLEX' else '0.0'
    return '''    const %s zero = %s;
    if (%s || incx != 1 || incy != 1 || alpha == zero) // alpha == 0: only y = beta y
        return false;''' % (T, z, extra)

# ---------------- Hermitian / symmetric matrix-vector ----------------
def symv(name, T, kind):
    cj = (lambda e: 'conj(%s)' % e) if T == 'COMPLEX' else (lambda e: e)
    rl = (lambda e: '%s.real()' % e) if T == 'COMPLEX' else (lambda e: e)
    if kind == 'full':
        sig = 'const char *uplo, INTEGER const n, %(T)s const alpha, %(T)s *a, INTEGER const lda, %(T)s *x, INTEGER const incx, %(T)s const beta, %(T)s *y, INTEGER const incy' % dict(T=T)
        args = ['uplo', 'n', 'alpha', 'a', 'lda', 'x', 'incx', 'beta', 'y', 'incy']
        A = lambda i, j: 'a[(%s - 1) + (%s - 1) * lda]' % (i, j)
        lo_u = lambda j: '1'; hi_l = lambda j: 'n'
        jlo_u = 'r0 + 1'; jhi_u = 'n'; jlo_l = '1'; jhi_l = 'r1'
        pre = ''
        extra = '(!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || lda < imax(1, n)'
        work = '4.0 * n * n'
    elif kind == 'band':
        sig = 'const char *uplo, INTEGER const n, INTEGER const k, %(T)s const alpha, %(T)s *a, INTEGER const lda, %(T)s *x, INTEGER const incx, %(T)s const beta, %(T)s *y, INTEGER const incy' % dict(T=T)
        args = ['uplo', 'n', 'k', 'alpha', 'a', 'lda', 'x', 'incx', 'beta', 'y', 'incy']
        # upper: element (i,j) = a[((kplus1 - j + i) - 1) + (j-1)*lda]; lower: a[((1 - j + i) - 1) + (j-1)*lda]
        A = None
        extra = '(!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || k < 0 || lda < k + 1'
        work = '4.0 * n * (k + 1)'
    else:  # packed
        sig = 'const char *uplo, INTEGER const n, %(T)s const alpha, %(T)s *ap, %(T)s *x, INTEGER const incx, %(T)s const beta, %(T)s *y, INTEGER const incy' % dict(T=T)
        args = ['uplo', 'n', 'alpha', 'ap', 'x', 'incx', 'beta', 'y', 'incy']
        extra = '(!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0'
        work = '4.0 * n * n'
    if kind == 'full':
        Au = Al = A
        du = A('j', 'j'); dl = A('j', 'j')
        ulo = '1'; lhi = 'n'
        uj0, uj1 = 'r0 + 1', 'n'; lj0, lj1 = '1', 'r1'
        decl = ''
    elif kind == 'band':
        Au = lambda i, j: 'a[((l + %s) - 1) + (%s - 1) * lda]' % (i, j)
        Al = Au
        du = 'a[(kplus1 - 1) + (j - 1) * lda]'; dl = 'a[(j - 1) * lda]'
        ulo = 'imax((INTEGER)1, j - k)'; lhi = 'imin(n, j + k)'
        uj0, uj1 = 'r0 + 1', 'imin(n, r1 + k)'; lj0, lj1 = 'imax((INTEGER)1, r0 + 1 - k)', 'r1'
        decl = '    const INTEGER kplus1 = k + 1;\n    INTEGER l = 0;\n'
    else:
        Au = lambda i, j: 'ap[(%s * (%s - 1) / 2 + %s) - 1]' % (j, j, i)
        Al = lambda i, j: 'ap[((%s - 1) * n - (%s - 1) * (%s - 2) / 2 + %s - %s + 1) - 1]' % (j, j, j, i, j)
        du = Au('j', 'j'); dl = Al('j', 'j')
        ulo = '1'; lhi = 'n'
        uj0, uj1 = 'r0 + 1', 'n'; lj0, lj1 = '1', 'r1'
        decl = ''
    lset_u = '            l = kplus1 - j;\n' if kind == 'band' else ''
    lset_l = '            l = 1 - j;\n' if kind == 'band' else ''
    rows = ROWDOC % (name, 'y') + '''void %(n)s_rows(%(sig)s, INTEGER const r0, INTEGER const r1)
{
%(beta)s    INTEGER j = 0;
    %(T)s temp1 = 0.0;
    %(T)s temp2 = 0.0;
%(decl)s    if (Mlsame(uplo, "U")) {
        // row i: column i (with the whole column above i), then columns j > i
        for (j = %(uj0)s; j <= %(uj1)s; j = j + 1) {
            temp1 = alpha * x[j - 1];
%(lset_u)s            for (i = imax(%(ulo)s, r0 + 1); i <= imin(j - 1, r1); i = i + 1) {
                y[i - 1] += temp1 * %(Aij_u)s;
            }
            if (j <= r1) {
                temp2 = zero;
                for (i = %(ulo)s; i <= j - 1; i = i + 1) {
                    temp2 += %(cAij_u)s * x[i - 1];
                }
                y[j - 1] += temp1 * %(du)s + alpha * temp2;
            }
        }
    } else {
        // row i: columns j < i, then column i (with the whole column below i)
        for (j = %(lj0)s; j <= %(lj1)s; j = j + 1) {
            temp1 = alpha * x[j - 1];
%(lset_l)s            if (j > r0) {
                y[j - 1] += temp1 * %(dl)s;
            }
            for (i = imax(j + 1, r0 + 1); i <= imin(%(lhi)s, r1); i = i + 1) {
                y[i - 1] += temp1 * %(Aij_l)s;
            }
            if (j > r0) {
                temp2 = zero;
                for (i = j + 1; i <= %(lhi)s; i = i + 1) {
                    temp2 += %(cAij_l)s * x[i - 1];
                }
                y[j - 1] += alpha * temp2;
            }
        }
    }
}
''' % dict(n=name, sig=sig, beta=beta_block(T), T=T, decl=decl, uj0=uj0, uj1=uj1, lj0=lj0, lj1=lj1,
           lset_u=lset_u, lset_l=lset_l, ulo=ulo, lhi=lhi, Aij_u=Au('i', 'j'), cAij_u=cj(Au('i', 'j')),
           du=rl(du), dl=rl(dl), Aij_l=Al('i', 'j'), cAij_l=cj(Al('i', 'j')))
    emit(name, sig, args, 'ranges of rows of y in parallel.', mvcheck(T, extra), work, 'n', rows)

symv('Chemv', 'COMPLEX', 'full')
symv('Chbmv', 'COMPLEX', 'band')
symv('Rsbmv', 'REAL', 'band')
symv('Chpmv', 'COMPLEX', 'packed')
symv('Rspmv', 'REAL', 'packed')

# ---------------- general band matrix-vector ----------------
def gbmv(name, T):
    sig = 'const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, %(T)s const alpha, %(T)s *a, INTEGER const lda, %(T)s *x, INTEGER const incx, %(T)s const beta, %(T)s *y, INTEGER const incy' % dict(T=T)
    args = ['trans', 'm', 'n', 'kl', 'ku', 'alpha', 'a', 'lda', 'x', 'incx', 'beta', 'y', 'incy']
    A = 'a[((k + i) - 1) + (j - 1) * lda]'
    if T == 'COMPLEX':
        tbody = '''            if (noconj) {
                for (i = imax((INTEGER)1, j - ku); i <= imin(m, j + kl); i = i + 1) {
                    temp += %s * x[i - 1];
                }
            } else {
                for (i = imax((INTEGER)1, j - ku); i <= imin(m, j + kl); i = i + 1) {
                    temp += conj(%s) * x[i - 1];
                }
            }
''' % (A, A)
        nc = '    bool noconj = Mlsame(trans, "T");\n'
        tr = '(!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C"))'
    else:
        tbody = '''            for (i = imax((INTEGER)1, j - ku); i <= imin(m, j + kl); i = i + 1) {
                temp += %s * x[i - 1];
            }
''' % A
        nc = ''
        tr = '(!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C"))'
    rows = ROWDOC % (name, 'y') + '''void %(n)s_rows(%(sig)s, INTEGER const r0, INTEGER const r1)
{
%(beta)s%(nc)s    const INTEGER kup1 = ku + 1;
    INTEGER j = 0;
    %(T)s temp = 0.0;
    INTEGER k = 0;
    if (Mlsame(trans, "N")) {
        // row i: the columns j of its band, in order
        for (j = imax((INTEGER)1, r0 + 1 - kl); j <= imin(n, r1 + ku); j = j + 1) {
            temp = alpha * x[j - 1];
            k = kup1 - j;
            for (i = imax(imax((INTEGER)1, j - ku), r0 + 1); i <= imin(imin(m, j + kl), r1); i = i + 1) {
                y[i - 1] += temp * %(A)s;
            }
        }
    } else {
        // y(j): column j
        for (j = r0 + 1; j <= r1; j = j + 1) {
            temp = zero;
            k = kup1 - j;
%(tbody)s            y[j - 1] += alpha * temp;
        }
    }
}
''' % dict(n=name, sig=sig, beta=beta_block(T), nc=nc, T=T, A=A, tbody=tbody)
    check = mvcheck(T, '%s || m <= 0 || n <= 0 || kl < 0 || ku < 0 || lda < kl + ku + 1' % tr)
    emit(name, sig, args, 'ranges of rows of y in parallel.', check, '4.0 * (ls(trans, "N") ? m : n) * (kl + ku + 1)', 'ls(trans, "N") ? m : n', rows)

gbmv('Rgbmv', 'REAL')
gbmv('Cgbmv', 'COMPLEX')

# ---------------- triangular matrix-vector (x := op(A) x) ----------------
def trmv(name, T, kind):
    cplx = T == 'COMPLEX'
    z = 'COMPLEX(0.0, 0.0)' if cplx else '0.0'
    if kind == 'full':
        sig = 'const char *uplo, const char *trans, const char *diag, INTEGER const n, %(T)s *a, INTEGER const lda, %(T)s *x, INTEGER const incx' % dict(T=T)
        args = ['uplo', 'trans', 'diag', 'n', 'a', 'lda', 'x', 'incx']
        Au = Al = lambda i, j: 'a[(%s - 1) + (%s - 1) * lda]' % (i, j)
        du = dl = Au('j', 'j')
        ulo = '1'; lhi = 'n'
        nu_j0, nu_j1 = 'r0 + 1', 'n'          # upper N columns
        nl_j1, nl_j0 = 'r1', '1'              # lower N columns (descending)
        decl = ''; lu = ''; ll = ''
        extra = 'n <= 0 || lda < imax(1, n)'
        work = '4.0 * n * n'
    elif kind == 'band':
        sig = 'const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, %(T)s *a, INTEGER const lda, %(T)s *x, INTEGER const incx' % dict(T=T)
        args = ['uplo', 'trans', 'diag', 'n', 'k', 'a', 'lda', 'x', 'incx']
        Au = Al = lambda i, j: 'a[((l + %s) - 1) + (%s - 1) * lda]' % (i, j)
        du = 'a[(kplus1 - 1) + (j - 1) * lda]'; dl = 'a[(j - 1) * lda]'
        ulo = 'imax((INTEGER)1, j - k)'; lhi = 'imin(n, j + k)'
        nu_j0, nu_j1 = 'r0 + 1', 'imin(n, r1 + k)'
        nl_j1, nl_j0 = 'r1', 'imax((INTEGER)1, r0 + 1 - k)'
        decl = '    const INTEGER kplus1 = k + 1;\n    INTEGER l = 0;\n'
        lu = '            l = kplus1 - j;\n'; ll = '            l = 1 - j;\n'
        extra = 'n <= 0 || k < 0 || lda < k + 1'
        work = '4.0 * n * (k + 1)'
    else:
        sig = 'const char *uplo, const char *trans, const char *diag, INTEGER const n, %(T)s *ap, %(T)s *x, INTEGER const incx' % dict(T=T)
        args = ['uplo', 'trans', 'diag', 'n', 'ap', 'x', 'incx']
        Au = lambda i, j: 'ap[(%s * (%s - 1) / 2 + %s) - 1]' % (j, j, i)
        Al = lambda i, j: 'ap[((%s - 1) * n - (%s - 1) * (%s - 2) / 2 + %s - %s + 1) - 1]' % (j, j, j, i, j)
        du = Au('j', 'j'); dl = Al('j', 'j')
        ulo = '1'; lhi = 'n'
        nu_j0, nu_j1 = 'r0 + 1', 'n'
        nl_j1, nl_j0 = 'r1', '1'
        decl = ''; lu = ''; ll = ''
        extra = 'n <= 0'
        work = '4.0 * n * n'
    lu2 = lu.replace('            ', '            ')
    def tdot(A, d, lo_loop):
        # lo_loop: the i loop header
        if cplx:
            return '''            if (noconj) {
                if (nounit) {
                    temp = temp * %(d)s;
                }
                %(loop)s {
                    temp += %(A)s * x0[i - 1];
                }
            } else {
                if (nounit) {
                    temp = temp * conj(%(d)s);
                }
                %(loop)s {
                    temp += conj(%(A)s) * x0[i - 1];
                }
            }
''' % dict(d=d, A=A, loop=lo_loop)
        return '''            if (nounit) {
                temp = temp * %(d)s;
            }
            %(loop)s {
                temp += %(A)s * x0[i - 1];
            }
''' % dict(d=d, A=A, loop=lo_loop)
    nc = '    bool noconj = Mlsame(trans, "T");\n' if cplx else ''
    rows = '''// %(n)s restricted to the elements r0 < i <= r1 of x (1-based); x0 is the
// original x (the reference reads x(j) before any update of it, so it is x0(j)).
void %(n)s_rows(%(sig)s, const %(T)s *x0, INTEGER const r0, INTEGER const r1)
{
%(nc)s    bool nounit = Mlsame(diag, "N");
    INTEGER j = 0;
    const %(T)s zero = %(z)s;
    %(T)s temp = 0.0;
    INTEGER i = 0;
%(decl)s    if (Mlsame(trans, "N")) {
        if (Mlsame(uplo, "U")) {
            for (j = %(nu_j0)s; j <= %(nu_j1)s; j = j + 1) {
%(lu)s                if (x0[j - 1] != zero) {
                    temp = x0[j - 1];
                    for (i = imax(%(ulo)s, r0 + 1); i <= imin(j - 1, r1); i = i + 1) {
                        x[i - 1] += temp * %(Aiju)s;
                    }
                    if (nounit && j <= r1) {
                        x[j - 1] = x[j - 1] * %(du)s;
                    }
                }
            }
        } else {
            for (j = %(nl_j1)s; j >= %(nl_j0)s; j = j - 1) {
%(ll)s                if (x0[j - 1] != zero) {
                    temp = x0[j - 1];
                    for (i = imin(%(lhi)s, r1); i >= imax(j + 1, r0 + 1); i = i - 1) {
                        x[i - 1] += temp * %(Aijl)s;
                    }
                    if (nounit && j > r0) {
                        x[j - 1] = x[j - 1] * %(dl)s;
                    }
                }
            }
        }
    } else {
        // x(j) from column j and the original x
        if (Mlsame(uplo, "U")) {
            for (j = r1; j >= r0 + 1; j = j - 1) {
                temp = x0[j - 1];
%(lu)s%(tu)s                x[j - 1] = temp;
            }
        } else {
            for (j = r0 + 1; j <= r1; j = j + 1) {
                temp = x0[j - 1];
%(ll)s%(tl)s                x[j - 1] = temp;
            }
        }
    }
}
''' % dict(n=name, sig=sig, T=T, z=z, nc=nc, decl=decl, nu_j0=nu_j0, nu_j1=nu_j1, nl_j0=nl_j0, nl_j1=nl_j1,
           lu=lu.replace('            ', '                ') if False else lu, ll=ll, ulo=ulo, lhi=lhi,
           Aiju=Au('i', 'j'), Aijl=Al('i', 'j'), du=du, dl=dl,
           tu=tdot(Au('i', 'j'), du, 'for (i = j - 1; i >= %s; i = i - 1)' % ulo),
           tl=tdot(Al('i', 'j'), dl, 'for (i = j + 1; i <= %s; i = i + 1)' % lhi))
    # indentation fix for the l = lines inside the N loops (16 spaces)
    rows = rows.replace('            for (j = %s; j <= %s; j = j + 1) {\n%s' % (nu_j0, nu_j1, lu),
                        '            for (j = %s; j <= %s; j = j + 1) {\n%s' % (nu_j0, nu_j1, lu.replace('            ', '                ')))
    rows = rows.replace('            for (j = %s; j >= %s; j = j - 1) {\n%s' % (nl_j1, nl_j0, ll),
                        '            for (j = %s; j >= %s; j = j - 1) {\n%s' % (nl_j1, nl_j0, ll.replace('            ', '                ')))
    check = '''    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || %s || incx != 1)
        return false;
    const std::vector<%s> x0(x, x + n);''' % (extra, T)
    sig_rows_args = args + ['x0.data()']
    emit(name, sig, sig_rows_args, 'ranges of x in parallel, from a copy of the original x.', check, work, 'n', rows)
    # the dispatcher passes the routine's own arguments
    disp = os.path.join(OUT, name + '.cpp')
    s = open(disp).read().replace(', x0.data()', '')
    open(disp, 'w').write(s)

trmv('Ctrmv', 'COMPLEX', 'full')
trmv('Rtbmv', 'REAL', 'band')
trmv('Ctbmv', 'COMPLEX', 'band')
trmv('Rtpmv', 'REAL', 'packed')
trmv('Ctpmv', 'COMPLEX', 'packed')

# ---------------- packed rank-1/2 updates: ranges of columns ----------------
KKU = '        kk = j * (j - 1) / 2 + 1;\n'
KKL = '        kk = (j - 1) * n - (j - 1) * (j - 2) / 2 + 1;\n'
def spr(name, T, two):
    cplx = T == 'COMPLEX'
    at = 'REAL' if (cplx and not two) else T
    if two:
        sig = 'const char *uplo, INTEGER const n, %(A)s const alpha, %(T)s *x, INTEGER const incx, %(T)s *y, INTEGER const incy, %(T)s *ap' % dict(A=at, T=T)
        args = ['uplo', 'n', 'alpha', 'x', 'incx', 'y', 'incy', 'ap']
    else:
        sig = 'const char *uplo, INTEGER const n, %(A)s const alpha, %(T)s *x, INTEGER const incx, %(T)s *ap' % dict(A=at, T=T)
        args = ['uplo', 'n', 'alpha', 'x', 'incx', 'ap']
    z = 'COMPLEX(0.0, 0.0)' if cplx else '0.0'
    if not cplx and not two:   # Rspr
        tdecl = '    REAL temp = 0.0;\n'
        up = '''            if (x[j - 1] != zero) {
                temp = alpha * x[j - 1];
                k = kk;
                for (i = 1; i <= j; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp;
                    k++;
                }
            }
'''
        lo = '''            if (x[j - 1] != zero) {
                temp = alpha * x[j - 1];
                k = kk;
                for (i = j; i <= n; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp;
                    k++;
                }
            }
'''
    elif not cplx:             # Rspr2
        tdecl = '    REAL temp1 = 0.0;\n    REAL temp2 = 0.0;\n'
        up = '''            if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                temp1 = alpha * y[j - 1];
                temp2 = alpha * x[j - 1];
                k = kk;
                for (i = 1; i <= j; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                    k++;
                }
            }
'''
        lo = up.replace('for (i = 1; i <= j; i = i + 1)', 'for (i = j; i <= n; i = i + 1)')
    elif not two:              # Chpr
        tdecl = '    COMPLEX temp = 0.0;\n'
        up = '''            if (x[j - 1] != zero) {
                temp = alpha * conj(x[j - 1]);
                k = kk;
                for (i = 1; i <= j - 1; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp;
                    k++;
                }
                ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1].real() + (x[j - 1] * temp).real();
            } else {
                ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1].real();
            }
'''
        lo = '''            if (x[j - 1] != zero) {
                temp = alpha * conj(x[j - 1]);
                ap[kk - 1] = ap[kk - 1].real() + (temp * x[j - 1]).real();
                k = kk + 1;
                for (i = j + 1; i <= n; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp;
                    k++;
                }
            } else {
                ap[kk - 1] = ap[kk - 1].real();
            }
'''
    else:                      # Chpr2
        tdecl = '    COMPLEX temp1 = 0.0;\n    COMPLEX temp2 = 0.0;\n'
        up = '''            if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                temp1 = alpha * conj(y[j - 1]);
                temp2 = conj(alpha * x[j - 1]);
                k = kk;
                for (i = 1; i <= j - 1; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                    k++;
                }
                ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1].real() + (x[j - 1] * temp1 + y[j - 1] * temp2).real();
            } else {
                ap[(kk + j - 1) - 1] = ap[(kk + j - 1) - 1].real();
            }
'''
        lo = '''            if ((x[j - 1] != zero) || (y[j - 1] != zero)) {
                temp1 = alpha * conj(y[j - 1]);
                temp2 = conj(alpha * x[j - 1]);
                ap[kk - 1] = ap[kk - 1].real() + (x[j - 1] * temp1 + y[j - 1] * temp2).real();
                k = kk + 1;
                for (i = j + 1; i <= n; i = i + 1) {
                    ap[k - 1] += x[i - 1] * temp1 + y[i - 1] * temp2;
                    k++;
                }
            } else {
                ap[kk - 1] = ap[kk - 1].real();
            }
'''
    rows = '''// %(n)s restricted to the columns r0 < j <= r1 of AP (1-based); kk is the
// position of the first element of column j, as the reference counts it.
void %(n)s_rows(%(sig)s, INTEGER const r0, INTEGER const r1)
{
    const %(T)s zero = %(z)s;
    INTEGER kk = 1;
    INTEGER j = 0;
%(tdecl)s    INTEGER k = 0;
    INTEGER i = 0;
    if (Mlsame(uplo, "U")) {
        for (j = r0 + 1; j <= r1; j = j + 1) {
    %(kku)s%(up)s        }
    } else {
        for (j = r0 + 1; j <= r1; j = j + 1) {
    %(kkl)s%(lo)s        }
    }
}
''' % dict(n=name, sig=sig, T=T, z=z, tdecl=tdecl, kku=KKU, kkl=KKL,
           up=up.replace('\n            ', '\n            ').replace('            ', '            '), lo=lo)
    azero = 'zero.real()' if (cplx and not two) else 'zero'
    check = '''    const %s zero = %s;
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || alpha == %s || incx != 1%s)
        return false;''' % (T, z, azero, ' || incy != 1' if two else '')
    emit(name, sig, args, 'ranges of columns of AP in parallel.', check, ('4.0' if two else '2.0') + ' * n * n', 'n', rows)

spr('Rspr', 'REAL', False)
spr('Rspr2', 'REAL', True)
spr('Chpr', 'COMPLEX', False)
spr('Chpr2', 'COMPLEX', True)

# the layout of the other sources
STYLE = "{BasedOnStyle: LLVM, IndentWidth: 4, ColumnLimit: 0, BreakBeforeBraces: Custom, BraceWrapping: {AfterFunction: true}, AllowShortFunctionsOnASingleLine: All, SortIncludes: false}"
if shutil.which("clang-format"):
    for r in "Chemv Chbmv Rsbmv Chpmv Rspmv Rgbmv Cgbmv Ctrmv Rtbmv Ctbmv Rtpmv Ctpmv Rspr Rspr2 Chpr Chpr2".split():
        subprocess.run(["clang-format", "-style=" + STYLE, "-i", os.path.join(OUT, r + "_omp.cpp")], check=True)

# The same range variants as host/device templates for the GPU
# (blas_dev/l2_rows_dev.h, instantiated with the device types of
# ../Rblas_dev_types.h or mpfr/cuda/Cblas_mpfr_dev_types.h): every GPU thread
# computes one element (one column for the packed updates).
import re
dev = ["// Generated by misc/gen_blas_l2_ranges.py; do not edit.",
       "// <R>_rows_dev: <R>_rows of mpblas/optimized/common/blas/<R>_omp.cpp as a",
       "// host/device template (Mlsame -> lsame, INTEGER -> mplapackint).", "",
       "#ifndef _MPLAPACK_L2_ROWS_DEV_H_", "#define _MPLAPACK_L2_ROWS_DEV_H_", "",
       "namespace mplapack_blas_dev {", ""]
for r in "Chemv Chbmv Rsbmv Chpmv Rspmv Rgbmv Cgbmv Ctrmv Rtbmv Ctbmv Rtpmv Ctpmv Rspr Rspr2 Chpr Chpr2".split():
    s = open(os.path.join(OUT, r + "_omp.cpp")).read()
    m = re.search(r"^void %s_rows\(.*?^}\n" % r, s, re.M | re.S)
    f = m.group(0)
    f = f.replace("void %s_rows(" % r, "template <class REAL, class COMPLEX> MPLAPACK_DEV_HD void %s_rows_dev(" % r, 1)
    f = re.sub(r"\bINTEGER\b", "mplapackint", f)
    f = re.sub(r"\bMlsame\(", "lsame(", f)
    f = re.sub(r"\bimax\(", "max<mplapackint>(", f)
    f = re.sub(r"\bimin\(", "min<mplapackint>(", f)
    dev.append(f)
dev += ["} // namespace mplapack_blas_dev", "", "#endif", ""]
open("mpblas/optimized/common/blas_dev/l2_rows_dev.h", "w").write("\n".join(dev))
