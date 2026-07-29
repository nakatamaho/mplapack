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

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpreal *a, int len) {
    mpreal tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpreal *a, int lda) {
    mpreal mtmp;
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

mpreal maxabs(mpreal a, mpreal b) {
    mpreal d = abs(a - b);
    return d;
}

mpreal max_solution_error(mplapackint n, mplapackint nrhs, mpreal *x, mplapackint ldx, mpreal *xexact, mplapackint ldxexact) {
    mpreal err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpreal d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpreal max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpreal *a, mplapackint lda, mpreal *x, mplapackint ldx, mpreal *b, mplapackint ldb) {
    mpreal err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpreal s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpreal d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 3, m = 2, p = 3, lda = n, ldb = n, info, lwork = -1;
    mpreal *a = new mpreal[lda * m];
    mpreal *b = new mpreal[ldb * p];
    mpreal *d = new mpreal[n];
    mpreal *x = new mpreal[m];
    mpreal *y = new mpreal[p];
    mpreal *xexact = new mpreal[m];
    mpreal *yexact = new mpreal[p];
    for (mplapackint i = 0; i < lda * m; i++)
        a[i] = 0;
    a[0] = 1;
    a[1] = 0;
    a[2] = 1;
    a[0 + lda] = 0;
    a[1 + lda] = 1;
    a[2 + lda] = 1;
    for (mplapackint i = 0; i < ldb * p; i++)
        b[i] = 0;
    for (mplapackint i = 0; i < n; i++)
        b[i + i * ldb] = 1;
    xexact[0] = 1;
    xexact[1] = 2;
    yexact[0] = mpreal(0.5);
    yexact[1] = mpreal(-0.5);
    yexact[2] = 1;
    for (mplapackint i = 0; i < n; i++)
        d[i] = a[i] * xexact[0] + a[i + lda] * xexact[1] + yexact[i];
    mpreal wk;
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, &wk, lwork, info);
    lwork = castINTEGER_mpfr(wk);
    mpreal *work = new mpreal[lwork];
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
