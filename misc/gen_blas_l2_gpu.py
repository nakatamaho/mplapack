#!/usr/bin/env python3
"""Generate the GPU side of the Level 2 range variants (misc/gen_blas_l2_ranges.py).

Writes the jobs (mpblas/optimized/common/Cblas_cuda_jobs_common.h: one
element of y or x, or one column of AP, per GPU thread), the host routines
<R>_gpu_t and the entry points <R>_<backend>_cuda
(mpblas/optimized/common/blas_cuda/Cblas_bridge_body.h) and the
dispatchers mpblas/optimized/common/blas_cuda/<R>.cpp of Chemv, Chbmv,
Rsbmv, Chpmv, Rspmv, Rgbmv, Cgbmv, Ctrmv, Rtbmv, Ctbmv, Rtpmv, Ctpmv, Rspr,
Rspr2, Chpr and Chpr2.  Run once on files without these routines (it
inserts, it does not replace).

usage: misc/gen_blas_l2_gpu.py   (from the top source directory)
"""
# spec: name, element type (C/R), chars, longs, scalars [(type,name)], operands [(name, host-arg, ld, rows, cols, download)], units, rows-call args, work, check
S = []
# alpha == 0: the reference only scales y (or returns), unlike the range variants
ZERO = {'C': 'COMPLEX(0.0, 0.0)', 'R': 'REAL(0.0)'}
def add(**k): S.append(k)
sym = lambda n, T, band, packed: add(
    n=n, T=T, chars=['uplo'], longs=['n'] + (['k', 'lda'] if band else ([] if packed else ['lda'])),
    scal=[(T, 'alpha'), (T, 'beta')],
    ops=[('A', 'ap' if packed else 'a', '1' if packed else 'lda', '1' if packed else ('k + 1' if band else 'n'), 'n * (n + 1) / 2' if packed else 'n', 'false'),
         ('X', 'x', '1', '1', 'n', 'false'), ('Y', 'y', '1', '1', 'n', 'true')],
    units='n',
    call='uplo, n, %salpha, A, %sX, 1, beta, Y, 1' % ('k, ' if band else '', '' if packed else 'lda, '),
    work='4.0 * n * %s' % ('(k + 1)' if band else 'n'),
    check='(!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == %s%s' % (ZERO[T], ' || k < 0 || lda < k + 1' if band else ('' if packed else ' || lda < std::max((INTEGER)1, n)')))
sym('Chemv', 'C', False, False); sym('Chbmv', 'C', True, False); sym('Rsbmv', 'R', True, False)
sym('Chpmv', 'C', False, True); sym('Rspmv', 'R', False, True)
for n, T in [('Rgbmv', 'R'), ('Cgbmv', 'C')]:
    add(n=n, T=T, chars=['tr'], longs=['m', 'n', 'kl', 'ku', 'lda'], scal=[(T, 'alpha'), (T, 'beta')],
        ops=[('A', 'a', 'lda', 'kl + ku + 1', 'n', 'false'), ('X', 'x', '1', '1', 'ls(trans, "N") ? n : m', 'false'), ('Y', 'y', '1', '1', 'ls(trans, "N") ? m : n', 'true')],
        units='lsame(tr, "N") ? m : n', call='tr, m, n, kl, ku, alpha, A, lda, X, 1, beta, Y, 1',
        work='4.0 * (ls(trans, "N") ? m : n) * (kl + ku + 1)',
        check='(!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || kl < 0 || ku < 0 || lda < kl + ku + 1 || incx != 1 || incy != 1 || alpha == ' + ZERO[T],
        trans_arg=True)
for n, T, kind in [('Ctrmv', 'C', 'full'), ('Rtbmv', 'R', 'band'), ('Ctbmv', 'C', 'band'), ('Rtpmv', 'R', 'packed'), ('Ctpmv', 'C', 'packed')]:
    band = kind == 'band'; packed = kind == 'packed'
    add(n=n, T=T, chars=['uplo', 'tr', 'diag'], longs=['n'] + (['k', 'lda'] if band else ([] if packed else ['lda'])), scal=[],
        ops=[('A', 'ap' if packed else 'a', '1' if packed else 'lda', '1' if packed else ('k + 1' if band else 'n'), 'n * (n + 1) / 2' if packed else 'n', 'false'),
             ('X0', 'x', '1', '1', 'n', 'false'), ('X', 'x', '1', '1', 'n', 'true')],
        units='n', call='uplo, tr, diag, n, %sA, %sX, 1, X0' % ('k, ' if band else '', '' if packed else 'lda, '),
        work='4.0 * n * %s' % ('(k + 1)' if band else 'n'),
        check='(!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1%s' % (' || k < 0 || lda < k + 1' if band else ('' if packed else ' || lda < std::max((INTEGER)1, n)')),
        trans_arg=True)
for n, T, two, at in [('Rspr', 'R', False, 'R'), ('Rspr2', 'R', True, 'R'), ('Chpr', 'C', False, 'R'), ('Chpr2', 'C', True, 'C')]:
    ops = [('X', 'x', '1', '1', 'n', 'false')] + ([('Y', 'y', '1', '1', 'n', 'false')] if two else []) + [('P', 'ap', '1', '1', 'n * (n + 1) / 2', 'true')]
    add(n=n, T=T, chars=['uplo'], longs=['n'], scal=[(at, 'alpha')], ops=ops, units='n',
        call='uplo, n, alpha, X, 1, %sP' % ('Y, 1, ' if two else ''), work='%s * n * n' % ('4.0' if two else '2.0'),
        check='(!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1%s || alpha == %s' % (' || incy != 1' if two else '', ZERO[at]))

def ty(t): return 'DR' if t == 'R' else 'DC'
jobs = []; bodies = []
for s in S:
    E = ty(s['T'])
    names = [o[0] for o in s['ops']]
    fields = ''
    if s['chars']: fields += '    char %s;\n' % ', '.join('%s[2]' % c for c in s['chars'])
    fields += '    long %s;\n' % ', '.join(s['longs'])
    for t, nm in s['scal']: fields += '    %s %s;\n' % (ty(t), nm)
    fields += '    E %s;\n' % ', '.join('*' + nm for nm in names)
    bind = ' '.join('%s = d[%d];' % (nm, i) for i, nm in enumerate(names))
    jobs.append('''template <class Ops> struct job_%(n)s {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef %(E)s E;
%(fields)s    MPLAPACK_DEV_HD long units() const { return %(units)s; }
    MPLAPACK_DEV_HD void bind(E **d) { %(bind)s }
    MPLAPACK_DEV_HD void operator()(long u) const { %(n)s_rows_dev<DR, DC>(%(call)s, u, u + 1); }
};
''' % dict(n=s['n'], E=E, fields=fields, units=s['units'], bind=bind, call=s['call']))
    # bridge body
    src = open('mpblas/optimized/common/blas/%s.cpp' % s['n']).read()
    import re
    sig = re.search(r'^void %s\((.*?)\)\n\{' % s['n'], src, re.M).group(1)
    fill = []
    for c in s['chars']:
        fill.append('ch(j.%s, %s);' % (c, 'trans' if c == 'tr' else c))
    fill.append(' '.join('j.%s = %s;' % (l, l) for l in s['longs']))
    for t, nm in s['scal']: fill.append('ctx.sc(j.%s, %s);' % (nm, nm))
    ops = ', '.join('ctx.template op<E>(%s, %s, %s, %s, %s)' % (h, ld, r, c, d) for nm, h, ld, r, c, d in s['ops'])
    bodies.append('''template <class Ops, class X> bool %(n)s_gpu_t(X &ctx, %(sig)s)
{
    if (%(check)s)
        return false;
    if (!ctx.eligible(false, %(work)s))
        return false;
    typedef job_%(n)s<Ops> J;
    typedef typename J::E E;
    J j;
    %(fill)s
    const gpu_operand<E> ops[] = {%(ops)s};
    return ctx.run(j, ops, %(no)d);
}
''' % dict(n=s['n'], sig=sig, check=s['check'], work=s['work'], fill=' '.join(fill), ops=ops, no=len(s['ops'])))
    s['sig'] = sig
    s['args'] = [a.strip().split()[-1].lstrip('*') for a in sig.split(',')]

J = 'mpblas/optimized/common/Cblas_cuda_jobs_common.h'
c = open(J).read()
anchor = '// Defined in' if '// Defined in' in c else 'template <class Job> int cgpu_run('
i = c.index('template <class Job> int cgpu_run(')
# put the jobs before the comment preceding cgpu_run
k = c.rfind('\n\n', 0, i) + 2
c = c[:k] + '// Level 2 range variants (blas_dev/l2_rows_dev.h): one element of y or x,\n// or one column of AP, per thread.\n\n' + ''.join(j + '\n' for j in jobs) + c[k:]
open(J, 'w').write(c)

B = 'mpblas/optimized/common/blas_cuda/Cblas_bridge_body.h'
c = open(B).read()
i = c.index('} // namespace mplapack_cblas_gpu')
c = c[:i] + '// Level 2 range variants (Level 2: eligible(false, ...))\n\n' + '\n'.join(bodies) + '\n' + c[i:]
lines = []
for s in S:
    lines.append('#define MPLAPACK_CBLAS_GPU_FUNCTION_%s(RUN) \\\n    bool MPLAPACK_CUDA_ENTRY(%s)(%s) { return RUN(%s_gpu_t, %s); }' % (s['n'], s['n'], s['sig'], s['n'], ', '.join(s['args'])))
i = c.index('#define MPLAPACK_CBLAS_GPU_FUNCTIONS(RUN)')
c = c[:i] + '\n'.join(lines) + '\n\n' + c[i:]
c = c.replace('    MPLAPACK_CBLAS_GPU_FUNCTION_Cher2(RUN)\n', '    MPLAPACK_CBLAS_GPU_FUNCTION_Cher2(RUN) \\\n' + ' \\\n'.join('    MPLAPACK_CBLAS_GPU_FUNCTION_%s(RUN)' % s['n'] for s in S) + '\n')
open(B, 'w').write(c)

# dispatchers
lic = open('mpblas/optimized/common/blas/Cher.cpp').read().split('*/', 1)[0] + '*/\n'
for s in S:
    a = ', '.join(s['args'])
    open('mpblas/optimized/common/blas_cuda/%s.cpp' % s['n'], 'w').write(lic + '''
// %(n)s of the CUDA libraries: eligible calls (Level 2: only with
// <P>_BLAS_CUDA_L12=1 / MPLAPACK_MPFR_CUDA_L12=1) run on the GPU
// (Cblas_bridge_body.h, ../blas_dev/l2_rows_dev.h), the others on the CPU
// (../blas/%(n)s_omp.cpp, then the reference code).  All of them give the
// reference result bit for bit.

#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"

bool %(n)s_omp(%(sig)s);
void %(n)s_reference_cpu(%(sig)s);

void %(n)s(%(sig)s)
{
    if (MPLAPACK_CUDA_ENTRY(%(n)s)(%(a)s))
        return;
    if (%(n)s_omp(%(a)s))
        return;
    %(n)s_reference_cpu(%(a)s);
}
''' % dict(n=s['n'], sig=s['sig'], a=a))
print(' '.join(s['n'] for s in S))
