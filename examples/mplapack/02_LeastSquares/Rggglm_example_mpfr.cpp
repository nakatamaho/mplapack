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
    mplapackint n = 3, m = 2, p = 3, lda = n, ldb = n, info, lwork = -1;
    mpfr_class *a = new mpfr_class[lda * m];
    mpfr_class *b = new mpfr_class[ldb * p];
    mpfr_class *d = new mpfr_class[n];
    mpfr_class *x = new mpfr_class[m];
    mpfr_class *y = new mpfr_class[p];
    mpfr_class *xexact = new mpfr_class[m];
    mpfr_class *yexact = new mpfr_class[p];
    for (mplapackint i = 0; i < lda * m; i++)
        a[i] = 0.0;
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 1.0;
    a[0 + lda] = 0.0;
    a[1 + lda] = 1.0;
    a[2 + lda] = 1.0;
    for (mplapackint i = 0; i < ldb * p; i++)
        b[i] = 0.0;
    for (mplapackint i = 0; i < n; i++)
        b[i + i * ldb] = 1.0;
    xexact[0] = 1.0;
    xexact[1] = 2.0;
    yexact[0] = mpfr_class(0.5);
    yexact[1] = mpfr_class(-0.5);
    yexact[2] = 1.0;
    for (mplapackint i = 0; i < n; i++)
        d[i] = a[i] * xexact[0] + a[i + lda] * xexact[1] + yexact[i];
    mpfr_class wk;
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, &wk, lwork, info);
    lwork = castINTEGER_mpfr(wk);
    mpfr_class *work = new mpfr_class[lwork];
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
    printf("x = "); printvec(x, m); printf("\n");
    printf("y = "); printvec(y, p); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(m, (mplapackint)1, x, m, xexact, m)); printf("\n");
    printf("max |y-y_exact| = "); printnum(max_solution_error(p, (mplapackint)1, y, p, yexact, p)); printf("\n");
    delete[] work;
    delete[] yexact;
    delete[] xexact;
    delete[] y;
    delete[] x;
    delete[] d;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
