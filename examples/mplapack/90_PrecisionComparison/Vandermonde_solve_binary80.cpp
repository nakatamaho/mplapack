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
mplapack_binary80_t max_solution_error(mplapackint n, mplapack_binary80_t *x, mplapack_binary80_t *xexact) {
    mplapack_binary80_t err = 0;
    for (mplapackint i = 0; i < n; i++) {
        mplapack_binary80_t d = abs(x[i] - xexact[i]);
        if (err < d)
            err = d;
    }
    return err;
}
int main() {
    mplapackint n = 15, lda = n, ldb = n, info;
    mplapack_binary80_t *a = new mplapack_binary80_t[n * n];
    mplapack_binary80_t *b = new mplapack_binary80_t[n];
    mplapack_binary80_t *xexact = new mplapack_binary80_t[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        xexact[j] = mplapack_binary80_t(j % 3 - 1);
    for (mplapackint i = 0; i < n; i++) {
        mplapack_binary80_t node = mplapack_binary80_t(i + 1);
        mplapack_binary80_t p = 1;
        for (mplapackint j = 0; j < n; j++) {
            a[i + j * lda] = p;
            p = p * node;
        }
    }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint j = 0; j < n; j++)
            b[i] = b[i] + a[i + j * lda] * xexact[j];
    }
    Rgesv(n, (mplapackint)1, a, lda, ipiv, b, ldb, info);
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, b, xexact)); printf("\n");
    delete[] ipiv;
    delete[] xexact;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
