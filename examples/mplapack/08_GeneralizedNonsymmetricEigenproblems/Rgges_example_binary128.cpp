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
bool select_none(mplapack_binary128_t ar, mplapack_binary128_t ai, mplapack_binary128_t beta) {
    return false;
}
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, sdim, info, lwork = -1;
    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *b = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *alphar = new mplapack_binary128_t[n];
    mplapack_binary128_t *alphai = new mplapack_binary128_t[n];
    mplapack_binary128_t *beta = new mplapack_binary128_t[n];
    mplapack_binary128_t *vsl = new mplapack_binary128_t[1];
    mplapack_binary128_t *vsr = new mplapack_binary128_t[1];
    bool *bwork = new bool[n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0.0;
        b[i] = 0.0;
    }
    a[0] = 1.0;
    a[4] = 2.0;
    a[8] = 3.0;
    b[0] = 1.0;
    b[4] = 1.0;
    b[8] = 0.0;
    mplapack_binary128_t wk;
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, &wk, lwork, bwork, info);
    lwork = castINTEGER_binary128(wk);
    mplapack_binary128_t *work = new mplapack_binary128_t[lwork];
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, work, lwork, bwork, info);
    printf("S = "); printmat(n, n, a, lda); printf("\n");
    printf("T = "); printmat(n, n, b, ldb); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        printf("lambda[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_binary128("E"))
            printf("Inf\n");
        else {
            printnum(alphar[i] / beta[i]);
            printf(" + "); printnum(alphai[i] / beta[i]); printf("i\n");
        }
    }
    delete[] work;
    delete[] bwork;
    delete[] vsr;
    delete[] vsl;
    delete[] beta;
    delete[] alphai;
    delete[] alphar;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
