//public domain
#include <mpblas_binary128.h>
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
    mplapackint n = 3;

    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *b = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *c = new mplapack_binary128_t[n * n];
    mplapack_binary128_t alpha, beta;

//setting A matrix
    a[0 + 0 * n] = 1.0;    a[0 + 1 * n] = 8.0;    a[0 + 2 * n] = 3.0;
    a[1 + 0 * n] = 2.0;    a[1 + 1 * n] = 10.0;   a[1 + 2 * n] = 8.0;
    a[2 + 0 * n] = 9.0;    a[2 + 1 * n] = -5.0;   a[2 + 2 * n] = -1.0;

    b[0 + 0 * n] = 9.0;    b[0 + 1 * n] = 8.0;    b[0 + 2 * n] = 3.0;
    b[1 + 0 * n] = 3.0;    b[1 + 1 * n] = -11.0;  b[1 + 2 * n] = 8.0;
    b[2 + 0 * n] = -8.0;   b[2 + 1 * n] = 6.0;    b[2 + 2 * n] = 1.0;

    c[0 + 0 * n] = 3.0;    c[0 + 1 * n] = 3.0;    c[0 + 2 * n] = -9.0;
    c[1 + 0 * n] = 8.0;    c[1 + 1 * n] = 4.0;    c[1 + 2 * n] = 8.0;
    c[2 + 0 * n] = 6.0;    c[2 + 1 * n] = 1.0;    c[2 + 2 * n] = -2.0;

    printf("# Rgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = 3.0;
    beta = -2.0;
    Rgemm("n", "n", n, n, n, alpha, a, n, b, n, beta, c, n);

    printf("alpha = "); printnum(alpha); printf("\n");
    printf("beta = "); printnum(beta); printf("\n");
    printf("ans ="); printmat(n, n, c, n); printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * a * b + beta * c \n");
    delete[]c;
    delete[]b;
    delete[]a;
}
