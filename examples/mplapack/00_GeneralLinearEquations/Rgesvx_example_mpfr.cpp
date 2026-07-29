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
    mplapackint n = 5, nrhs = 1, lda = n, ldb = n, info;
    mpreal *a = new mpreal[n * n];
    mpreal *af = new mpreal[n * n];
    mpreal *b = new mpreal[n];
    mpreal *x = new mpreal[n];
    mpreal *xexact = new mpreal[n];
    mpreal *r = new mpreal[n];
    mpreal *c = new mpreal[n];
    mpreal *ferr = new mpreal[nrhs];
    mpreal *berr = new mpreal[nrhs];
    mpreal *work = new mpreal[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    char equed[2]; equed[0] = 'N'; equed[1] = '\0';
    mpreal rcond;
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) a[i + j * lda] = mpreal(1.0) / mpreal(i + j + 1);
    for (mplapackint i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? mpreal(1.0) : mpreal(-1.0);
    for (mplapackint i = 0; i < n; i++) { b[i] = 0; for (mplapackint k = 0; k < n; k++) b[i] = b[i] + a[i + k * lda] * xexact[k]; }
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Rgesvx("N", "N", n, nrhs, a, lda, af, lda, ipiv, equed, r, c, b, ldb, x, ldb, rcond, ferr, berr, work, iwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("rcond = "); printnum(rcond); printf("\n");
    printf("ferr = "); printvec(ferr, nrhs); printf("\n");
    printf("berr = "); printvec(berr, nrhs); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, x, ldb, xexact, n)); printf("\n");
    delete[] ipiv; delete[] iwork; delete[] work; delete[] berr; delete[] ferr; delete[] c; delete[] r; delete[] xexact; delete[] x; delete[] b; delete[] af; delete[] a;
    return info != 0 ? 1 : 0;
}
