//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(mplapack_binary128_t rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(mplapack_binary128_t *a, int len) {
    mplapack_binary128_t tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mplapack_binary128_t *a, int lda)
{
    mplapack_binary128_t mtmp;

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
mplapack_binary128_t maxabs(mplapack_binary128_t a, mplapack_binary128_t b) {
    mplapack_binary128_t d = abs(a - b);
    return d;
}

mplapack_binary128_t max_solution_error(mplapackint n, mplapackint nrhs, mplapack_binary128_t *x, mplapackint ldx, mplapack_binary128_t *xexact, mplapackint ldxexact) {
    mplapack_binary128_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mplapack_binary128_t d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mplapack_binary128_t max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mplapack_binary128_t *a, mplapackint lda, mplapack_binary128_t *x, mplapackint ldx, mplapack_binary128_t *b, mplapackint ldb) {
    mplapack_binary128_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mplapack_binary128_t s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mplapack_binary128_t d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 5, nrhs = 1, lda = n, ldb = n, info;
    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *af = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *b = new mplapack_binary128_t[n];
    mplapack_binary128_t *x = new mplapack_binary128_t[n];
    mplapack_binary128_t *xexact = new mplapack_binary128_t[n];
    mplapack_binary128_t *r = new mplapack_binary128_t[n];
    mplapack_binary128_t *c = new mplapack_binary128_t[n];
    mplapack_binary128_t *ferr = new mplapack_binary128_t[nrhs];
    mplapack_binary128_t *berr = new mplapack_binary128_t[nrhs];
    mplapack_binary128_t *work = new mplapack_binary128_t[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    char equed[2];
    equed[0] = 'N';
    equed[1] = '\0';
    mplapack_binary128_t rcond;
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) a[i + j * lda] = mplapack_binary128_t(1.0) / mplapack_binary128_t(i + j + 1);
    for (mplapackint i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? mplapack_binary128_t(1.0) : mplapack_binary128_t(-1.0);
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0.0;
        for (mplapackint k = 0; k < n; k++)
            b[i] = b[i] + a[i + k * lda] * xexact[k];
    }
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Rgesvx("N", "N", n, nrhs, a, lda, af, lda, ipiv, equed, r, c, b, ldb, x, ldb, rcond, ferr, berr, work, iwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("rcond = "); printnum(rcond); printf("\n");
    printf("ferr = "); printvec(ferr, nrhs); printf("\n");
    printf("berr = "); printvec(berr, nrhs); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, x, ldb, xexact, n)); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] berr;
    delete[] ferr;
    delete[] c;
    delete[] r;
    delete[] xexact;
    delete[] x;
    delete[] b;
    delete[] af;
    delete[] a;
    return info != 0 ? 1 : 0;
}
