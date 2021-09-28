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
    mplapackint lwork, info;

    _Float128 *a = new _Float128[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
//https://academic.oup.com/qjmam/article/1/1/253/1883468
//The Quarterly Journal of Mechanics and Applied Mathematics, Volume 1, Issue 1, 1948, Pages 253280, https://doi.org/10.1093/qjmam/1.1.253
//https://babel.hathitrust.org/cgi/pt?id=mdp.39015023899019&view=1up&seq=30
//Marcus, M. (1960). Basic theorems in matrix theory. Washington: U.S. Govt. Print. Off.

    a[0 + 0 * n] = 5;    a[0 + 1 * n] = 7;    a[0 + 2 * n] = 6;      a[0 + 3 * n] = 5;
    a[1 + 0 * n] = 7;    a[1 + 1 * n] = 10;   a[1 + 2 * n] = 8;      a[1 + 3 * n] = 7;
    a[2 + 0 * n] = 6;    a[2 + 1 * n] = 8;    a[2 + 2 * n] = 10;     a[2 + 3 * n] = 9;
    a[3 + 0 * n] = 5;    a[3 + 1 * n] = 7;    a[3 + 2 * n] = 9;      a[3 + 3 * n] = 10;

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    _Float128 *work = new _Float128[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER__Float128 (work[0]);
    delete[]work;
    work = new _Float128[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("a * ainv - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
