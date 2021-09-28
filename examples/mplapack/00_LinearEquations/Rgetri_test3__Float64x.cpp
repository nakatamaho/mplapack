//public domain
#include <mpblas__Float64x.h>
#include <mplapack__Float64x.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define FLOAT64X_FORMAT "%+25.21Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

void printnum(_Float64x rtmp)
{
    printf(FLOAT64X_FORMAT, rtmp);
    return;
}

//Matlab/Octave format
void printvec(_Float64x *a, int len) {
    _Float64x tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, _Float64x *a, int lda)
{
    _Float64x mtmp;

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

    _Float64x *a = new _Float64x[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
//https://www.tuhh.de/ti3/paper/rump/NiRuOi11.pdf
//Nonlinear Theory and Its Applications, IEICE, vol. 2, no. 2, pp. 226-245
//DOI: 10.1588/nolta.2.226
    a[0 + 0 * n] = 300;  a[0 + 1 * n] = -590; a[0 + 2 * n] = 850;    a[0 + 3 * n] = -561;
    a[1 + 0 * n] = 1;    a[1 + 1 * n] = -1;   a[1 + 2 * n] = 0.0;    a[1 + 3 * n] = 0.0;
    a[2 + 0 * n] = 0.0;  a[2 + 1 * n] = 1;    a[2 + 2 * n] = -1;     a[2 + 3 * n] = 0.0;
    a[3 + 0 * n] = 0.0;  a[3 + 1 * n] = 0.0;  a[3 + 2 * n] = 1;      a[3 + 3 * n] = -1;

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    _Float64x *work = new _Float64x[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER__Float64x (work[0]);
    delete[]work;
    work = new _Float64x[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("a * ainv - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
