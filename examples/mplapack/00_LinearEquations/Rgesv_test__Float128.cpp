//public domain
#include <mpblas__Float128.h>
#include <mplapack__Float128.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(_Float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp);
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(_Float128 *a, int len) {
    _Float128 tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, _Float128 *a, int lda)
{
    _Float128 mtmp;

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
    _Float128 *a = new _Float128[n * n];
    _Float128 *b = new _Float128[n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
///Collection of Matrices for Testing Computational Algorithms 1969/1/21 Robert T. Gregory,  David L. Karney
    a[0 + 0 * n] = 1;   a[0 + 1 * n] = -2;  a[0 + 2 * n] = 3;    a[0 + 3 * n] = 1;
    a[1 + 0 * n] = -2;  a[1 + 1 * n] = 1;   a[1 + 2 * n] = -2;   a[1 + 3 * n] = -1;
    a[2 + 0 * n] = 3;   a[2 + 1 * n] = -2;  a[2 + 2 * n] = 1;    a[2 + 3 * n] = 5;
    a[3 + 0 * n] = 1;   a[3 + 1 * n] = -1;  a[3 + 2 * n] = 5;    a[3 + 3 * n] = 3;

    b[0] = 3;
    b[1] = -4;
    b[2] = 7;
    b[3] = 8;

//answer is [1, 1, 1, 1]
    printf("a ="); printmat(n, n, a, n); printf("\n");

//Solve linear equation
    Rgesv(n, (mplapackint)1, a, n, ipiv, b, n, info);

    printf("x ="); printvec(b, n); printf("\n");
    delete[]ipiv;
    delete[]b;
    delete[]a;
}
