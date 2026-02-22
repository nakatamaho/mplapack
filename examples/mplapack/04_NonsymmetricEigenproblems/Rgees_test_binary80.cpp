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
bool rselect(mplapack_binary80_t ar, mplapack_binary80_t ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    mplapack_binary80_t *a = new mplapack_binary80_t[n * n];
    mplapack_binary80_t *vs = new mplapack_binary80_t[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 3 * n;
    mplapack_binary80_t *wr = new mplapack_binary80_t[n];
    mplapack_binary80_t *wi = new mplapack_binary80_t[n];
    mplapack_binary80_t *work = new mplapack_binary80_t[lwork];
    bool bwork[n];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = -2.0;     a[0 + 1 * n] = 2.0;   a[0 + 2 * n] = 2.0;    a[0 + 3 * n] = 2.0;
    a[1 + 0 * n] = -3.0;     a[1 + 1 * n] = 3.0;   a[1 + 2 * n] = 2.0;    a[1 + 3 * n] = 2.0;
    a[2 + 0 * n] = -2.0;     a[2 + 1 * n] = 0.0;   a[2 + 2 * n] = 4.0;    a[2 + 3 * n] = 2.0;
    a[3 + 0 * n] = -1.0;     a[3 + 1 * n] = 0.0;   a[3 + 2 * n] = 0.0;    a[3 + 3 * n] = 5.0;

    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "S", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
    printf("eig(a)\n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }
    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vs;
    delete[] a;
}
