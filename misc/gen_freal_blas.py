#!/usr/bin/env python3
# Generates the drivers of the fixed-precision MPFR BLAS
# (mpblas/optimized/mpfr/openmp/Rfreal_blas_drivers.h) and one source file
# per precision (Rfreal_blas_<PB>.cpp) from the routine table
# MPLAPACK_FREAL_ROUTINES of Rfreal_blas.h and the operand descriptions
# below.  A driver copies the operands (the referenced region, with its access:
# F_IN read, F_OUT written only, F_INOUT both) and runs the routine of the
# instance (Rfreal_blas_inst.h) on the copies.
#   usage: misc/gen_freal_blas.py   (from the top directory)
import re

D = 'mpblas/optimized/mpfr/openmp/'
PRECISIONS = list(range(128, 1025, 64)) + [1536, 2048]

# body of <name>_freal; CALL(args) runs <name>_omp, else the reference code
B = {}
# ---- real Level 3
B['Rsyrk'] = '''const bool notr = lsm(trans, "N");
fcall cl((double)n * n * k / 2);
const REAL al = cl.re(alpha), be = cl.re(beta);
auto &A = cl.mat(a, notr ? n : k, notr ? k : n, lda, rg_all(), F_IN);
auto &C = cl.mat(cc, n, n, ldc, rg_tri(lsm(uplo, "U")), rw(be));
return cl.run([&] { CALL((uplo, trans, n, k, al, A.p, A.ld, be, C.p, C.ld)); });'''
B['Rsyr2k'] = '''const bool notr = lsm(trans, "N");
fcall cl((double)n * n * k);
const REAL al = cl.re(alpha), be = cl.re(beta);
auto &A = cl.mat(a, notr ? n : k, notr ? k : n, lda, rg_all(), F_IN);
auto &Bm = cl.mat(b, notr ? n : k, notr ? k : n, ldb, rg_all(), F_IN);
auto &C = cl.mat(cc, n, n, ldc, rg_tri(lsm(uplo, "U")), rw(be));
return cl.run([&] { CALL((uplo, trans, n, k, al, A.p, A.ld, Bm.p, Bm.ld, be, C.p, C.ld)); });'''
SYMM = '''const long na = lsm(side, "L") ? m : n;
fcall cl((double)m * n * na);
const T al = cl.S(alpha), be = cl.S(beta);
auto &A = cl.mat(a, na, na, lda, rg_tri(lsm(uplo, "U")), F_IN);
auto &Bm = cl.mat(b, m, n, ldb, rg_all(), F_IN);
auto &C = cl.mat(CC, m, n, ldc, rg_all(), rw(be));
return cl.run([&] { CALL((side, uplo, m, n, al, A.p, A.ld, Bm.p, Bm.ld, be, C.p, C.ld)); });'''
B['Rsymm'] = SYMM.replace('T al', 'REAL al').replace('cl.S', 'cl.re').replace('CC', 'cc')
TRXM = '''const long na = lsm(side, "L") ? m : n;
fcall cl((double)m * n * na / 2);
const T al = cl.S(alpha);
auto &A = cl.mat(a, na, na, lda, rg_tri(lsm(uplo, "U"), lsm(diag, "U")), F_IN);
auto &Bm = cl.mat(b, m, n, ldb, rg_all(), rw(al));
return cl.run([&] { CALL((side, uplo, transa, diag, m, n, al, A.p, A.ld, Bm.p, Bm.ld)); });'''
B['Rtrmm'] = TRXM.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
B['Rtrsm'] = B['Rtrmm']
GEMMTR = '''const bool nota = lsm(transa, "N"), notb = lsm(transb, "N");
fcall cl((double)n * n * k / 2);
const T al = cl.S(alpha), be = cl.S(beta);
auto &A = cl.mat(a, nota ? n : k, nota ? k : n, lda, rg_all(), F_IN);
auto &Bm = cl.mat(b, notb ? k : n, notb ? n : k, ldb, rg_all(), F_IN);
auto &C = cl.mat(c, n, n, ldc, rg_tri(lsm(uplo, "U")), rw(be));
return cl.run([&] { CALL((uplo, transa, transb, n, k, al, A.p, A.ld, Bm.p, Bm.ld, be, C.p, C.ld)); });'''
B['Rgemmtr'] = GEMMTR.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
# ---- real Level 2
GEMV = '''const bool notr = lsm(trans, "N");
fcall cl((double)m * n);
const T al = cl.S(alpha), be = cl.S(beta);
auto &A = cl.mat(a, m, n, lda, rg_all(), F_IN);
auto &X = cl.vec(x, notr ? n : m, incx, F_IN);
auto &Y = cl.vec(y, notr ? m : n, incy, rw(be));
return cl.run([&] { CALL((trans, m, n, al, A.p, A.ld, X.p, X.ld, be, Y.p, Y.ld)); });'''
B['Rgemv'] = GEMV.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SYMV = '''fcall cl((double)n * n);
const T al = cl.S(alpha), be = cl.S(beta);
auto &A = cl.mat(a, n, n, lda, rg_tri(lsm(uplo, "U")), F_IN);
auto &X = cl.vec(x, n, incx, F_IN);
auto &Y = cl.vec(y, n, incy, rw(be));
return cl.run([&] { CALL((uplo, n, al, A.p, A.ld, X.p, X.ld, be, Y.p, Y.ld)); });'''
B['Rsymv'] = SYMV.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
TRMV = '''fcall cl((double)n * n / 2);
auto &A = cl.mat(a, n, n, lda, rg_tri(lsm(uplo, "U"), lsm(diag, "U")), F_IN);
auto &X = cl.vec(x, n, incx, F_INOUT);
return cl.run([&] { CALL((uplo, trans, diag, n, A.p, A.ld, X.p, X.ld)); });'''
B['Rtrmv'] = TRMV
B['Rtrsv'] = TRMV
GER = '''fcall cl((double)m * n);
const T al = cl.S(alpha);
auto &X = cl.vec(x, m, incx, F_IN);
auto &Y = cl.vec(y, n, incy, F_IN);
auto &A = cl.mat(a, m, n, lda, rg_all(), F_INOUT);
return cl.run([&] { CALL((m, n, al, X.p, X.ld, Y.p, Y.ld, A.p, A.ld)); });'''
B['Rger'] = GER.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SYR = '''fcall cl((double)n * n / 2);
const T al = cl.S(alpha);
auto &X = cl.vec(x, n, incx, F_IN);
auto &A = cl.mat(a, n, n, lda, rg_tri(lsm(uplo, "U")), F_INOUT);
return cl.run([&] { CALL((uplo, n, al, X.p, X.ld, A.p, A.ld)); });'''
B['Rsyr'] = SYR.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SYR2 = '''fcall cl((double)n * n);
const T al = cl.S(alpha);
auto &X = cl.vec(x, n, incx, F_IN);
auto &Y = cl.vec(y, n, incy, F_IN);
auto &A = cl.mat(a, n, n, lda, rg_tri(lsm(uplo, "U")), F_INOUT);
return cl.run([&] { CALL((uplo, n, al, X.p, X.ld, Y.p, Y.ld, A.p, A.ld)); });'''
B['Rsyr2'] = SYR2.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
GBMV = '''const bool notr = lsm(trans, "N");
fcall cl((double)n * (kl + ku + 1));
const T al = cl.S(alpha), be = cl.S(beta);
auto &A = cl.mat(a, kl + ku + 1, n, lda, rg_band(kl, ku, m), F_IN);
auto &X = cl.vec(x, notr ? n : m, incx, F_IN);
auto &Y = cl.vec(y, notr ? m : n, incy, rw(be));
return cl.run([&] { CALL((trans, m, n, kl, ku, al, A.p, A.ld, X.p, X.ld, be, Y.p, Y.ld)); });'''
B['Rgbmv'] = GBMV.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SBMV = '''const bool upper = lsm(uplo, "U");
fcall cl((double)n * (2 * k + 1));
const T al = cl.S(alpha), be = cl.S(beta);
auto &A = cl.mat(a, k + 1, n, lda, upper ? rg_band(0, k, n) : rg_band(k, 0, n), F_IN);
auto &X = cl.vec(x, n, incx, F_IN);
auto &Y = cl.vec(y, n, incy, rw(be));
return cl.run([&] { CALL((uplo, n, k, al, A.p, A.ld, X.p, X.ld, be, Y.p, Y.ld)); });'''
B['Rsbmv'] = SBMV.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SPMV = '''fcall cl((double)n * n);
const T al = cl.S(alpha), be = cl.S(beta);
auto &P = cl.packed(ap, n, F_IN);
auto &X = cl.vec(x, n, incx, F_IN);
auto &Y = cl.vec(y, n, incy, rw(be));
return cl.run([&] { CALL((uplo, n, al, P.p, X.p, X.ld, be, Y.p, Y.ld)); });'''
B['Rspmv'] = SPMV.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SPR = '''fcall cl((double)n * n / 2);
const T al = cl.S(alpha);
auto &X = cl.vec(x, n, incx, F_IN);
auto &P = cl.packed(ap, n, F_INOUT);
return cl.run([&] { CALL((uplo, n, al, X.p, X.ld, P.p)); });'''
B['Rspr'] = SPR.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
SPR2 = '''fcall cl((double)n * n);
const T al = cl.S(alpha);
auto &X = cl.vec(x, n, incx, F_IN);
auto &Y = cl.vec(y, n, incy, F_IN);
auto &P = cl.packed(ap, n, F_INOUT);
return cl.run([&] { CALL((uplo, n, al, X.p, X.ld, Y.p, Y.ld, P.p)); });'''
B['Rspr2'] = SPR2.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
TBMV = '''const bool upper = lsm(uplo, "U"), unit = lsm(diag, "U");
fcall cl((double)n * (k + 1));
auto &A = cl.mat(a, k + 1, n, lda, upper ? rg_band(0, k, n, unit) : rg_band(k, 0, n, unit), F_IN);
auto &X = cl.vec(x, n, incx, F_INOUT);
return cl.run([&] { CALL((uplo, trans, diag, n, k, A.p, A.ld, X.p, X.ld)); });'''
B['Rtbmv'] = TBMV
B['Rtbsv'] = TBMV
TPMV = '''fcall cl((double)n * n / 2);
auto &P = cl.packed(ap, n, F_IN);
auto &X = cl.vec(x, n, incx, F_INOUT);
return cl.run([&] { CALL((uplo, trans, diag, n, P.p, X.p, X.ld)); });'''
B['Rtpmv'] = TPMV
B['Rtpsv'] = TPMV
# ---- real Level 1
B['Rscal'] = '''if (incx <= 0) // nothing to do
    return false;
fcall cl((double)n);
const REAL a = cl.re(da);
auto &X = cl.vec(dx, n, incx, F_INOUT);
return cl.run([&] { CALL((n, a, X.p, X.ld)); });'''
B['Rrot'] = '''fcall cl(2.0 * n);
const REAL c = cl.re(c_), sn = cl.re(s);
auto &X = cl.vec(dx, n, incx, F_INOUT);
auto &Y = cl.vec(dy, n, incy, F_INOUT);
return cl.run([&] { CALL((n, X.p, X.ld, Y.p, Y.ld, c, sn)); });'''
B['Rasum'] = '''if (incx <= 0)
    return false;
fcall cl((double)n);
auto &X = cl.vec(dx, n, incx, F_IN);
cl.ok = cl.ok && prec_ok(result);
REAL r;
if (!cl.run([&] { if (!Rasum_ns::Rasum_omp(n, X.p, X.ld, r)) r = Rasum_ns::Rasum_reference_cpu(n, X.p, X.ld); }))
    return false;
store(result, r);
return true;'''
B['Rdot'] = '''fcall cl((double)n);
auto &X = cl.vec(dx, n, incx, F_IN);
auto &Y = cl.vec(dy, n, incy, F_IN);
cl.ok = cl.ok && prec_ok(result);
REAL r;
if (!cl.run([&] { if (!Rdot_ns::Rdot_par(n, X.p, X.ld, Y.p, Y.ld, r)) r = Rdot_ns::Rdot_reference_cpu(n, X.p, X.ld, Y.p, Y.ld); }))
    return false;
store(result, r);
return true;'''
B['Raxpy'] = '''fcall cl((double)n);
const REAL a = cl.re(da);
auto &X = cl.vec(dx, n, incx, F_IN);
auto &Y = cl.vec(dy, n, incy, F_INOUT);
return cl.run([&] {
    // the elements are independent: parts of the vectors in parallel
    if (X.ld != 1 || Y.ld != 1 || !mplapack_blas::blas_split(2.0 * n, n, [&](long p0, long p1) { Raxpy_ns::Raxpy_reference_cpu(p1 - p0, a, X.p + p0, 1, Y.p + p0, 1); }))
        Raxpy_ns::Raxpy_reference_cpu(n, a, X.p, X.ld, Y.p, Y.ld);
});'''
# ---- complex: the same shapes
def cx(t, name):
    return t.replace('T al', 'COMPLEX al').replace('cl.S', 'cl.cx').replace('REAL al', 'COMPLEX al').replace('cl.re(', 'cl.cx(')
B['Cgemm'] = '''const bool nota = lsm(transa, "N"), notb = lsm(transb, "N");
fcall cl(4.0 * m * n * k);
const COMPLEX al = cl.cx(alpha), be = cl.cx(beta);
auto &A = cl.mat(a, nota ? m : k, nota ? k : m, lda, rg_all(), F_IN);
auto &Bm = cl.mat(b, notb ? k : n, notb ? n : k, ldb, rg_all(), F_IN);
auto &C = cl.mat(c, m, n, ldc, rg_all(), rw(be));
return cl.run([&] { CALL((transa, transb, m, n, k, al, A.p, A.ld, Bm.p, Bm.ld, be, C.p, C.ld)); });'''
B['Csymm'] = cx(SYMM, 'Csymm').replace('CC', 'c')
B['Chemm'] = B['Csymm']
B['Ctrmm'] = cx(TRXM, 'Ctrmm')
B['Ctrsm'] = B['Ctrmm']
B['Csyrk'] = cx(B['Rsyrk'], 'Csyrk').replace('cc', 'c').replace('(double)n * n * k / 2', '2.0 * n * n * k')
B['Cherk'] = B['Rsyrk'].replace('cc', 'c').replace('(double)n * n * k / 2', '2.0 * n * n * k')
B['Csyr2k'] = cx(B['Rsyr2k'], 'Csyr2k').replace('cc', 'c').replace('(double)n * n * k', '4.0 * n * n * k')
B['Cher2k'] = B['Csyr2k'].replace('be = cl.cx(beta)', 'be2 = cl.cx(alpha)').replace('const COMPLEX al = cl.cx(alpha), be2 = cl.cx(alpha);',
                                                                              'const COMPLEX al = cl.cx(alpha);\nconst REAL be = cl.re(beta);')
B['Cgemmtr'] = cx(GEMMTR, 'Cgemmtr').replace('(double)n * n * k / 2', '2.0 * n * n * k')
B['Cgemv'] = cx(GEMV, 'Cgemv')
B['Cgbmv'] = cx(GBMV, 'Cgbmv')
B['Chemv'] = cx(SYMV, 'Chemv')
B['Chbmv'] = cx(SBMV, 'Chbmv')
B['Chpmv'] = cx(SPMV, 'Chpmv')
B['Cgeru'] = cx(GER, 'Cgeru')
B['Cgerc'] = B['Cgeru']
B['Cher'] = SYR.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
B['Cher2'] = cx(SYR2, 'Cher2')
B['Chpr'] = SPR.replace('T al', 'REAL al').replace('cl.S', 'cl.re')
B['Chpr2'] = cx(SPR2, 'Chpr2')
B['Ctrmv'] = TRMV
B['Ctbmv'] = TBMV
B['Ctpmv'] = TPMV
for s, m in [('Ctrsv', TRMV), ('Ctbsv', TBMV), ('Ctpsv', TPMV)]:
    B[s] = m
B['Caxpy'] = '''fcall cl(4.0 * n);
const COMPLEX a = cl.cx(za);
auto &X = cl.vec(zx, n, incx, F_IN);
auto &Y = cl.vec(zy, n, incy, F_INOUT);
return cl.run([&] { CALL((n, a, X.p, X.ld, Y.p, Y.ld)); });'''
B['Cscal'] = '''if (incx <= 0)
    return false;
fcall cl(4.0 * n);
const COMPLEX a = cl.cx(za);
auto &X = cl.vec(zx, n, incx, F_INOUT);
return cl.run([&] { CALL((n, a, X.p, X.ld)); });'''
B['CRscal'] = B['Cscal'].replace('const COMPLEX a = cl.cx(za);', 'const REAL a = cl.re(da);').replace('4.0 * n', '2.0 * n')
B['CRrot'] = '''fcall cl(8.0 * n);
const REAL cr = cl.re(c), sr = cl.re(s);
auto &X = cl.vec(zx, n, incx, F_INOUT);
auto &Y = cl.vec(zy, n, incy, F_INOUT);
return cl.run([&] { CALL((n, X.p, X.ld, Y.p, Y.ld, cr, sr)); });'''
DOT = '''fcall cl(4.0 * n);
auto &X = cl.vec(zx, n, incx, F_IN);
auto &Y = cl.vec(zy, n, incy, F_IN);
cl.ok = cl.ok && prec_ok(result);
COMPLEX r;
if (!cl.run([&] { if (!NAME_ns::NAME_omp(n, X.p, X.ld, Y.p, Y.ld, r)) r = NAME_ns::NAME_reference_cpu(n, X.p, X.ld, Y.p, Y.ld); }))
    return false;
store(result, r);
return true;'''
B['Cdotc'] = DOT.replace('NAME', 'Cdotc')
B['Cdotu'] = DOT.replace('NAME', 'Cdotu')
B['RCasum'] = '''if (incx <= 0)
    return false;
fcall cl(2.0 * n);
auto &X = cl.vec(zx, n, incx, F_IN);
cl.ok = cl.ok && prec_ok(result);
REAL r;
if (!cl.run([&] { if (!RCasum_ns::RCasum_omp(n, X.p, X.ld, r)) r = RCasum_ns::RCasum_reference_cpu(n, X.p, X.ld); }))
    return false;
store(result, r);
return true;'''

# Level 1 and 2 routines: one or two operations per element copied (see
# fcall in Rfreal_blas_inst.h); the others are Level 3
L1 = 'Rscal Rrot Rasum Rdot Raxpy Caxpy Cscal CRscal CRrot Cdotc Cdotu RCasum'.split()
L3 = 'Rsyrk Rsyr2k Rsymm Rtrmm Rtrsm Rgemmtr Cgemm Csymm Chemm Ctrmm Ctrsm Csyrk Cherk Csyr2k Cher2k Cgemmtr'.split()

# the LAPACK routines: the matrix is copied once, the routine of the
# instance runs on it (its BLAS calls on the copies too)
GETRF = '''if (WINO)
    return false;
fcall cl((double)m * n * std::min(m, n), 3);
auto &A = cl.mat(a, m, n, lda, rg_all(), F_INOUT);
INTEGER inf = 0;
if (!cl.run([&] { NS::NAME(m, n, A.p, A.ld, ipiv, inf); }))
    return false;
info = inf;
return true;'''
TRI = '''if (WINO)
    return false;
fcall cl((double)n * n * n / 3, 3);
auto &A = cl.mat(a, n, n, lda, rg_tri(lsm(uplo, "U")UNIT), F_INOUT);
INTEGER inf = 0;
if (!cl.run([&] { NS::NAME(ARGS); }))
    return false;
info = inf;
return true;'''
# a Winograd Rgemm (opt-in) does not give the reference result: the MPFR code then
WINO_R = 'mplapack_blas::blas_env_long("MPLAPACK_MPFR_GEMM_WINOGRAD_CUTOFF", 0) > 0'
LB = {}
for r, c in [('R', 'false'), ('C', 'false')]:
    w = WINO_R if r == 'R' else 'false'
    LB[r + 'getrf'] = GETRF.replace('NS', r + 'getrf_ns').replace('NAME', r + 'getrf').replace('WINO', w)
    LB[r + 'getrf2'] = GETRF.replace('NS', r + 'getrf_ns').replace('NAME', r + 'getrf2').replace('WINO', w)
    LB[r + 'potrf'] = TRI.replace('NS', r + 'potrf_ns').replace('NAME', r + 'potrf').replace('ARGS', 'uplo, n, A.p, A.ld, inf').replace('UNIT', '').replace('WINO', w)
    LB[r + 'trtri'] = TRI.replace('NS', r + 'trtri_ns').replace('NAME', r + 'trtri').replace('ARGS', 'uplo, diag, n, A.p, A.ld, inf').replace('UNIT', ', lsm(diag, "U")').replace('WINO', w)
    LB[r + 'lauum'] = TRI.replace('NS', r + 'lauum_ns').replace('NAME', r + 'lauum').replace('ARGS', 'uplo, n, A.p, A.ld, inf').replace('UNIT', '').replace('WINO', w)
for k in list(LB):
    LB[k] = LB[k].replace('if (false)\n    return false;\n', '')

table = open(D + 'Rfreal_blas.h').read()
rows = re.findall(r'^    X\((\w+), \((.*?)\), \((.*?)\)\)', table, re.M)
lrows = re.findall(r'^    X\((\w+), \((.*?)\), \((.*?)\)\)', table.split('MPLAPACK_FREAL_LAPACK_ROUTINES(X)')[1], re.M)
rows = [r for r in rows if r not in lrows]
assert len(lrows) == len(LB), (lrows, LB.keys())
assert len(rows) == len(B), (len(rows), len(B), set(B) ^ set(r[0] for r in rows))

licence = open(D + 'Rfreal_state.h').read().split('*/')[0] + '*/\n'
out = [licence, '// GENERATED by misc/gen_freal_blas.py; do not edit.\n//',
       '// The drivers of the fixed-precision MPFR BLAS (Rfreal_blas.h), one per',
       '// routine, inside Rfreal_blas_inst.h: the operands with their referenced',
       '// regions and accesses, and the call of <name>_omp (else the reference code)',
       '// of the instance on the copies.\n', 'namespace FRB_NS {\n']
for name, params, args in rows:
    body = B[name].replace('CALL(', 'FRB_CALL(%s, ' % name)
    level = 1 if name in L1 else 3 if name in L3 else 2
    body = re.sub(r'fcall cl\((.*?)\);', lambda m: 'fcall cl(%s, %d);' % (m.group(1), level), body)
    assert ', %d);' % level in body, name
    out.append('bool %s_freal(%s)\n{\n%s\n}\n' % (name, params, '\n'.join(('    ' + l) if l else l for l in body.split('\n'))))
for name, params, args in lrows:
    body = LB[name]
    out.append('bool %s_freal(%s)\n{\n%s\n}\n' % (name, params, '\n'.join(('    ' + l) if l else l for l in body.split('\n'))))
out.append('} // namespace FRB_NS\n')
text = '\n'.join(out)
text = text.replace('namespace FRB_NS {\n', '#define FRB_CALL(name, args)                \\\n    do {                                    \\\n        if (!name##_ns::name##_omp args)    \\\n            name##_ns::name##_reference_cpu args; \\\n    } while (0)\n\nnamespace FRB_NS {\n', 1)
open(D + 'Rfreal_blas_drivers.h', 'w').write(text)

# BLAS of fr<PB> / fc<PB> under their names (for the LAPACK code of the
# instance): <name>_omp, else the reference code
names = [licence, '// GENERATED by misc/gen_freal_blas.py; do not edit.\n//',
         '// The BLAS of REAL = fr<PB> and COMPLEX = fc<PB> under their own names, in',
         '// the namespace of Rfreal_blas_inst.h, for the LAPACK code compiled there.\n',
         'namespace FRB_NS {\n']
for name, params, args0 in rows:
    args = '(' + args0 + ')'
    fp = params.replace('mpfr_class', 'REAL').replace('mpc_class', 'COMPLEX').replace('mplapackint', 'INTEGER')
    if name == 'Raxpy':
        names.append('inline void Raxpy(%s) { Raxpy_ns::Raxpy_reference_cpu%s; }' % (fp, args))
        continue
    m = re.match(r'(.*), (REAL|COMPLEX) &result$', fp)
    if m:
        a = args.replace(', result)', ')')
        omp = 'Rdot_ns::Rdot_par' if name == 'Rdot' else '%s_ns::%s_omp' % (name, name)
        names.append('inline %s %s(%s)\n{\n    %s r;\n    if (!%s%s)\n        r = %s_ns::%s_reference_cpu%s;\n    return r;\n}' %
                     (m.group(2), name, m.group(1), m.group(2), omp, args.replace('result', 'r'), name, name, a))
        continue
    names.append('inline void %s(%s)\n{\n    if (!%s_ns::%s_omp%s)\n        %s_ns::%s_reference_cpu%s;\n}' % (name, fp, name, name, args, name, name, args))
names.append('} // namespace FRB_NS\n')
open(D + 'Rfreal_blas_names.h', 'w').write('\n'.join(names))

for p in PRECISIONS:
    open(D + 'Rfreal_blas_%d.cpp' % p, 'w').write(licence + '''
// GENERATED by misc/gen_freal_blas.py; do not edit.
// The fixed-precision MPFR BLAS (Rfreal_blas.h) for the precision %d.

#define FRB_PB %d
#include "Rfreal_blas_inst.h"
''' % (p, p))
