//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpfr_class rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpfr_class rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpfr_class *a, int len) {
    mpfr_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpfr_class *a, int lda) {
    mpfr_class mtmp;
    printf("[ ");
    for (int i = 0; i < n; i++) {
        printf("[ ");
        for (int j = 0; j < m; j++) {
            mtmp = a[i + j * lda];
            printnum(mtmp);
            if (j < m - 1)
                printf(", ");
        }
        if (i < n - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}

mpfr_class maxabs(mpfr_class a, mpfr_class b) {
    mpfr_class d = abs(a - b);
    return d;
}

mpfr_class max_solution_error(mplapackint n, mplapackint nrhs, mpfr_class *x, mplapackint ldx, mpfr_class *xexact, mplapackint ldxexact) {
    mpfr_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpfr_class d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpfr_class max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpfr_class *a, mplapackint lda, mpfr_class *x, mplapackint ldx, mpfr_class *b, mplapackint ldb) {
    mpfr_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpfr_class s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpfr_class d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint m = 3, n = 2, p = 1, lda = m, ldb = p, info, lwork = -1;
    mpfr_class *a = new mpfr_class[lda * n];
    mpfr_class *bmat = new mpfr_class[ldb * n];
    mpfr_class *c = new mpfr_class[m];
    mpfr_class *d = new mpfr_class[p];
    mpfr_class *x = new mpfr_class[n];
    mpfr_class *xexact = new mpfr_class[n];
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 1.0;
    a[0 + lda] = 0.0;
    a[1 + lda] = 1.0;
    a[2 + lda] = 1.0;
    bmat[0] = 1.0;
    bmat[0 + ldb] = 1.0;
    xexact[0] = 1.0;
    xexact[1] = 2.0;
    for (mplapackint i = 0; i < m; i++)
        c[i] = a[i] * xexact[0] + a[i + lda] * xexact[1];
    d[0] = 3.0;
    mpfr_class wk;
    Rgglse(m, n, p, a, lda, bmat, ldb, c, d, x, &wk, lwork, info);
    lwork = castINTEGER_mpfr(wk);
    mpfr_class *work = new mpfr_class[lwork];
    Rgglse(m, n, p, a, lda, bmat, ldb, c, d, x, work, lwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("constraint B*x-d = "); printnum(x[0] + x[1] - d[0]); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, (mplapackint)1, x, n, xexact, n)); printf("\n");
    delete[] work;
    delete[] xexact;
    delete[] x;
    delete[] d;
    delete[] c;
    delete[] bmat;
    delete[] a;
    return info != 0 ? 1 : 0;
}
