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
void sort_real(mplapackint n, mplapack_binary80_t *x) {
    for (mplapackint i = 0; i < n; i++)
        for (mplapackint j = i + 1; j < n; j++)
            if (x[j] < x[i]) {
                mplapack_binary80_t t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
}
int main() {
    mplapackint n = 20, lda = n, ldv = 1, info, lwork = -1;
    mplapack_binary80_t *coef = new mplapack_binary80_t[n + 1];
    for (mplapackint i = 0; i <= n; i++)
        coef[i] = 0;
    coef[0] = 1;
    for (mplapackint k = 1; k <= n; k++) {
        for (mplapackint j = k; j >= 1; j--)
            coef[j] = coef[j] - mplapack_binary80_t(k) * coef[j - 1];
    }
    mplapack_binary80_t *a = new mplapack_binary80_t[n * n];
    for (mplapackint i = 0; i < n * n; i++)
        a[i] = 0;
    for (mplapackint i = 1; i < n; i++)
        a[i + (i - 1) * lda] = 1;
    for (mplapackint j = 0; j < n; j++)
        a[j + (n - 1) * lda] = -coef[n - j] / coef[0];
    mplapack_binary80_t *wr = new mplapack_binary80_t[n];
    mplapack_binary80_t *wi = new mplapack_binary80_t[n];
    mplapack_binary80_t *vl = new mplapack_binary80_t[1];
    mplapack_binary80_t *vr = new mplapack_binary80_t[1];
    mplapack_binary80_t wk;
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_binary80(wk);
    mplapack_binary80_t *work = new mplapack_binary80_t[lwork];
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, work, lwork, info);
    sort_real(n, wr);
    mplapack_binary80_t maxerr = 0;
    for (mplapackint i = 0; i < n; i++) {
        mplapack_binary80_t err = abs(wr[i] - mplapack_binary80_t(i + 1));
        if (maxerr < err)
            maxerr = err;
        printf("root[%ld] = ", (long)i); printnum(wr[i]); printf(", error = "); printnum(err); printf("\n");
    }
    printf("max root error = "); printnum(maxerr); printf("\n");
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] wi;
    delete[] wr;
    delete[] a;
    delete[] coef;
    return info != 0 ? 1 : 0;
}
