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
int main()
{
    mplapackint n = 4;
    mplapackint info;
    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *b = new mplapack_binary128_t[n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
///Collection of Matrices for Testing Computational Algorithms 1969/1/21 Robert T. Gregory,  David L. Karney
    a[0 + 0 * n] = 1.0;   a[0 + 1 * n] = -2.0;  a[0 + 2 * n] = 3.0;    a[0 + 3 * n] = 1.0;
    a[1 + 0 * n] = -2.0;  a[1 + 1 * n] = 1.0;   a[1 + 2 * n] = -2.0;   a[1 + 3 * n] = -1.0;
    a[2 + 0 * n] = 3.0;   a[2 + 1 * n] = -2.0;  a[2 + 2 * n] = 1.0;    a[2 + 3 * n] = 5.0;
    a[3 + 0 * n] = 1.0;   a[3 + 1 * n] = -1.0;  a[3 + 2 * n] = 5.0;    a[3 + 3 * n] = 3.0;

    b[0] = 3.0;
    b[1] = -4.0;
    b[2] = 7.0;
    b[3] = 8.0;

//answer is [1, 1, 1, 1]
    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printvec(b, n); printf("\n");

//Solve linear equation
    Rgesv(n, (mplapackint)1, a, n, ipiv, b, n, info);

    printf("x ="); printvec(b, n); printf("\n");
    printf("a*x'\n");
    delete[]ipiv;
    delete[]b;
    delete[]a;
}
