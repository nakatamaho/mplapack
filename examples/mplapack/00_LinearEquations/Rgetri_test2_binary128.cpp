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
int main()
{
    mplapackint n = 4;
    mplapackint lwork, info;

    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
//https://www.tuhh.de/ti3/paper/rump/NiRuOi11.pdf
//Nonlinear Theory and Its Applications, IEICE, vol. 2, no. 2, pp. 226-245
//DOI: 10.1588/nolta.2.226
    a[0 + 0 * n] = 17;   a[0 + 1 * n] = -864; a[0 + 2 * n] = 716;    a[0 + 3 * n] = -799;
    a[1 + 0 * n] = 1;    a[1 + 1 * n] = -50;  a[1 + 2 * n] = 0.0;    a[1 + 3 * n] = 0.0;
    a[2 + 0 * n] = 0.0;  a[2 + 1 * n] = 1;    a[2 + 2 * n] = -50;    a[2 + 3 * n] = 0.0;
    a[3 + 0 * n] = 0.0;  a[3 + 1 * n] = 0.0;  a[3 + 2 * n] = 1;      a[3 + 3 * n] = -50;

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    mplapack_binary128_t *work = new mplapack_binary128_t[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_binary128 (work[0]);
    delete[]work;
    work = new mplapack_binary128_t[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("a * ainv - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
