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
int main() {
    mplapackint m = 2, n = 3, p = 2, k, l, lda = m, ldb = p, ldu = m, ldv = p, ldq = n, info, lwork = -1;
    mplapack_binary80_t *a = new mplapack_binary80_t[lda * n];
    mplapack_binary80_t *b = new mplapack_binary80_t[ldb * n];
    mplapack_binary80_t *alpha = new mplapack_binary80_t[n];
    mplapack_binary80_t *beta = new mplapack_binary80_t[n];
    mplapack_binary80_t *u = new mplapack_binary80_t[ldu * m];
    mplapack_binary80_t *v = new mplapack_binary80_t[ldv * p];
    mplapack_binary80_t *q = new mplapack_binary80_t[ldq * n];
    mplapackint *iwork = new mplapackint[n];
    for (mplapackint i = 0; i < lda * n; i++)
        a[i] = 0;
    for (mplapackint i = 0; i < ldb * n; i++)
        b[i] = 0;
    a[0] = 1;
    a[1 + lda] = 2;
    a[0 + 2 * lda] = 1;
    b[0] = 1;
    b[1 + ldb] = 3;
    b[0 + 2 * ldb] = 1;
    mplapack_binary80_t wk;
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, &wk, lwork, iwork, info);
    lwork = castINTEGER_binary80(wk);
    mplapack_binary80_t *work = new mplapack_binary80_t[lwork];
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
    printf("k = %ld, l = %ld\n", (long)k, (long)l);
    printf("alpha = "); printvec(alpha, n); printf("\n");
    printf("beta = "); printvec(beta, n); printf("\n");
    for (mplapackint i = k; i < k + l; i++) {
        printf("gsv[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_binary80("E"))
            printf("Inf\n");
        else {
            printnum(alpha[i] / beta[i]);
            printf("\n");
        }
    }
    delete[] work;
    delete[] iwork;
    delete[] q;
    delete[] v;
    delete[] u;
    delete[] beta;
    delete[] alpha;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
