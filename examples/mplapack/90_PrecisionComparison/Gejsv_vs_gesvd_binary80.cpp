//public domain
#include <mpblas_binary80.h>
#include <mplapack_binary80.h>
#include <iostream>
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
    mplapackint m = 8, n = 4, lda = m, info, lwork = -1;
    mplapack_binary80_t *a = new mplapack_binary80_t[m * n];
    mplapack_binary80_t *b = new mplapack_binary80_t[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            mplapack_binary80_t scale = pow(mplapack_binary80_t(10.0), mplapack_binary80_t(i - 4));
            a[i + j * lda] = scale / (mplapack_binary80_t(i + j + 1));
            b[i + j * lda] = a[i + j * lda];
        }
    mplapack_binary80_t *s1 = new mplapack_binary80_t[n];
    mplapack_binary80_t *s2 = new mplapack_binary80_t[n];
    mplapack_binary80_t *u = new mplapack_binary80_t[1];
    mplapack_binary80_t *vt = new mplapack_binary80_t[1];
    mplapack_binary80_t wk;
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_binary80(wk);
    mplapack_binary80_t *work = new mplapack_binary80_t[lwork];
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, info);
    delete[] work;
    mplapackint *iwork = new mplapackint[m + 3 * n + 10];
    lwork = 5 * n * n + 9 * n + m;
    if (lwork < 2 * m + n)
        lwork = 2 * m + n;
    if (lwork < 4 * n + n * n)
        lwork = 4 * n + n * n;
    if (lwork < 7)
        lwork = 7;
    work = new mplapack_binary80_t[lwork];
    Rgejsv("G", "N", "N", "N", "N", "N", m, n, b, lda, s2, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, iwork, info);
    printf("Rgesvd singular values = "); printvec(s1, n); printf("\n");
    printf("Rgejsv singular values = "); printvec(s2, n); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        mplapack_binary80_t rel = abs(s1[i] - s2[i]) / (abs(s2[i]) + Rlamch_binary80("S"));
        printf("relative difference[%ld] = ", (long)i); printnum(rel); printf("\n");
    }
    delete[] work;
    delete[] iwork;
    delete[] vt;
    delete[] u;
    delete[] s2;
    delete[] s1;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
