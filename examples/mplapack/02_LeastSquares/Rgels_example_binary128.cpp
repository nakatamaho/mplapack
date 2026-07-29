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

void set_problem(mplapackint m, mplapackint n, mplapack_binary128_t *a, mplapackint lda, mplapack_binary128_t *b, mplapackint ldb, mplapack_binary128_t *xexact) {
    for (mplapackint i = 0; i < m; i++) {
        a[i + 0 * lda] = 1;
        a[i + 1 * lda] = i;
    }
    xexact[0] = 1;
    xexact[1] = 2;
    for (mplapackint i = 0; i < m; i++)
        b[i] = a[i + 0 * lda] * xexact[0] + a[i + 1 * lda] * xexact[1];
}

int main() {
    mplapackint m = 4, n = 2, nrhs = 1, lda = m, ldb = m, info, lwork = -1;
    mplapack_binary128_t *a = new mplapack_binary128_t[lda * n];
    mplapack_binary128_t *aorg = new mplapack_binary128_t[lda * n];
    mplapack_binary128_t *b = new mplapack_binary128_t[ldb];
    mplapack_binary128_t *borg = new mplapack_binary128_t[ldb];
    mplapack_binary128_t *xexact = new mplapack_binary128_t[n];
    set_problem(m, n, a, lda, b, ldb, xexact);
    for (mplapackint i = 0; i < lda * n; i++)
        aorg[i] = a[i];
    for (mplapackint i = 0; i < ldb; i++)
        borg[i] = b[i];
    mplapack_binary128_t wk;
    Rgels("N", m, n, nrhs, a, lda, b, ldb, &wk, lwork, info);
    lwork = castINTEGER_binary128(wk);
    mplapack_binary128_t *work = new mplapack_binary128_t[lwork];
    Rgels("N", m, n, nrhs, a, lda, b, ldb, work, lwork, info);
    printf("A = "); printmat(m, n, aorg, lda); printf("\n");
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(m, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");
    delete[] work;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
