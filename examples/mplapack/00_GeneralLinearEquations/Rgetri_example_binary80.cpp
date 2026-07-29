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
int main()
{
    mplapackint n = 4;
    mplapackint lwork, info;

    mplapack_binary80_t *a = new mplapack_binary80_t[n * n];
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
    mplapack_binary80_t *work = new mplapack_binary80_t[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_binary80 (work[0]);
    delete[]work;
    work = new mplapack_binary80_t[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("a * ainv - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
