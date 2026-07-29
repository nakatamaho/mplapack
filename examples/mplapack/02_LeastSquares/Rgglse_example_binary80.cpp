//public domain
#include <mpblas_binary80.h>
#include <mplapack_binary80.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#define BUFLEN 1024
void printnum(mplapack_binary80_t rtmp)
{
    int width = 25;
    char buf[BUFLEN];
#if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_STRFROMF64X
    strfromf64x(buf, sizeof(buf), "%.21e", rtmp);
#elif MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL
    snprintf(buf, sizeof(buf), "%*.21Le", width, rtmp);
#else
    #error "unsupported binary80 type"
#endif
    if (rtmp >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(mplapack_binary80_t *a, int len) {
    mplapack_binary80_t tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mplapack_binary80_t *a, int lda)
{
    mplapack_binary80_t mtmp;

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
mplapack_binary80_t maxabs(mplapack_binary80_t a, mplapack_binary80_t b) {
    mplapack_binary80_t d = abs(a - b);
    return d;
}

mplapack_binary80_t max_solution_error(mplapackint n, mplapackint nrhs, mplapack_binary80_t *x, mplapackint ldx, mplapack_binary80_t *xexact, mplapackint ldxexact) {
    mplapack_binary80_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mplapack_binary80_t d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mplapack_binary80_t max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mplapack_binary80_t *a, mplapackint lda, mplapack_binary80_t *x, mplapackint ldx, mplapack_binary80_t *b, mplapackint ldb) {
    mplapack_binary80_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mplapack_binary80_t s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mplapack_binary80_t d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint m = 3, n = 2, p = 1, lda = m, ldb = p, info, lwork = -1;
    mplapack_binary80_t *a = new mplapack_binary80_t[lda * n];
    mplapack_binary80_t *bmat = new mplapack_binary80_t[ldb * n];
    mplapack_binary80_t *c = new mplapack_binary80_t[m];
    mplapack_binary80_t *d = new mplapack_binary80_t[p];
    mplapack_binary80_t *x = new mplapack_binary80_t[n];
    mplapack_binary80_t *xexact = new mplapack_binary80_t[n];
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
    mplapack_binary80_t wk;
    Rgglse(m, n, p, a, lda, bmat, ldb, c, d, x, &wk, lwork, info);
    lwork = castINTEGER_binary80(wk);
    mplapack_binary80_t *work = new mplapack_binary80_t[lwork];
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
