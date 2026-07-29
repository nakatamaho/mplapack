//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
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
    mplapackint m = 3, n = 2, p = 1, lda = m, ldb = p, info, lwork = -1;
    mplapack_binary128_t *a = new mplapack_binary128_t[lda * n];
    mplapack_binary128_t *bmat = new mplapack_binary128_t[ldb * n];
    mplapack_binary128_t *c = new mplapack_binary128_t[m];
    mplapack_binary128_t *d = new mplapack_binary128_t[p];
    mplapack_binary128_t *x = new mplapack_binary128_t[n];
    mplapack_binary128_t *xexact = new mplapack_binary128_t[n];
    a[0] = 1;
    a[1] = 0;
    a[2] = 1;
    a[0 + lda] = 0;
    a[1 + lda] = 1;
    a[2 + lda] = 1;
    bmat[0] = 1;
    bmat[0 + ldb] = 1;
    xexact[0] = 1;
    xexact[1] = 2;
    for (mplapackint i = 0; i < m; i++)
        c[i] = a[i] * xexact[0] + a[i + lda] * xexact[1];
    d[0] = 3;
    mplapack_binary128_t wk;
    Rgglse(m, n, p, a, lda, bmat, ldb, c, d, x, &wk, lwork, info);
    lwork = castINTEGER_binary128(wk);
    mplapack_binary128_t *work = new mplapack_binary128_t[lwork];
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
