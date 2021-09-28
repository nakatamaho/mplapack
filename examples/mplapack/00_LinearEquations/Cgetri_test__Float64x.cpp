//public domain
#include <mpblas__Float64x.h>
#include <mplapack__Float64x.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define FLOAT64X_FORMAT "%+25.21Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

void printnum(_Float64x rtmp) { printf(FLOAT64X_FORMAT, rtmp); return;}
void printnum(std::complex<_Float64x> ctmp) { printf(FLOAT64X_FORMAT FLOAT64X_FORMAT "i", ctmp.real(), ctmp.imag()); }

//Matlab/Octave format
template <class X> void printvec(X *a, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

template <class X> void printmat(int n, int m, X *a, int lda)
{
    X mtmp;

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
//taking from Collection of Matrices for Testing Computational Algorithms 1969 Robert T. Gregory, David L. Karney pp.30
int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    std::complex<_Float64x> *a = new std::complex<_Float64x>[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix


    a[0 + 0 * n] = std::complex<_Float64x>(1.0, 0.0);   a[0 + 1 * n] = std::complex<_Float64x>(1.0, 2.0);    a[0 + 2 * n] = std::complex<_Float64x>(2.0, 10.0);
    a[1 + 0 * n] = std::complex<_Float64x>(1.0, 1.0);   a[1 + 1 * n] = std::complex<_Float64x>(0.0, 3.0);    a[1 + 2 * n] = std::complex<_Float64x>(-5.0, 14.0);
    a[2 + 0 * n] = std::complex<_Float64x>(1.0, 1.0);   a[2 + 1 * n] = std::complex<_Float64x>(0.0, 5.0);    a[2 + 2 * n] = std::complex<_Float64x>(-8.0, 20.0);

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    std::complex<_Float64x> *work = new std::complex<_Float64x>[1];

    Cgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER__Float64x (work[0].real());
    delete[]work;
    work = new std::complex<_Float64x>[std::max(1, (int) lwork)];

//inverse matrix
    Cgetrf(n, n, a, n, ipiv, info);
    Cgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("ainv * a - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
