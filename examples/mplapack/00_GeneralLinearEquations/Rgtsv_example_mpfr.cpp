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
    mplapackint n = 6, nrhs = 1, ldb = n, info;
    mpreal *dl = new mpreal[n - 1];
    mpreal *d = new mpreal[n];
    mpreal *du = new mpreal[n - 1];
    mpreal *a = new mpreal[n * n];
    mpreal *b = new mpreal[n];
    mpreal *borg = new mpreal[n];
    mpreal *xexact = new mpreal[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0;
    for (mplapackint i = 0; i < n; i++) {
        d[i] = 2;
        xexact[i] = i + 1;
        a[i + i * n] = 2;
    }
    for (mplapackint i = 0; i < n - 1; i++) {
        dl[i] = -1;
        du[i] = -1;
        a[i + 1 + i * n] = -1;
        a[i + (i + 1) * n] = -1;
    }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint k = 0; k < n; k++) b[i] = b[i] + a[i + k * n] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n, n, a, n); printf("\n");
    Rgtsv(n, nrhs, dl, d, du, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, a, n, b, ldb, borg, ldb)); printf("\n");
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] a;
    delete[] du;
    delete[] d;
    delete[] dl;
    return info != 0 ? 1 : 0;
}
