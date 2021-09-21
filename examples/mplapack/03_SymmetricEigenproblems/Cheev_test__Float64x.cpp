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
int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    std::complex<_Float64x> *A = new std::complex<_Float64x>[n * n];
    _Float64x *w = new _Float64x[n];
    _Float64x *rwork = new _Float64x[3 * n - 1];

//setting A matrix
    A[0 + 0 * n] = 2.0;               A[0 + 1 * n] = std::complex<_Float64x>(0.0, -1.0);    A[0 + 2 * n] = 0.0;
    A[1 + 0 * n] = std::complex<_Float64x>(0.0, 1.0); A[1 + 1 * n] = 2.0;                   A[1 + 2 * n] = 0.0;
    A[2 + 0 * n] = 0.0;               A[2 + 1 * n] = 0.0;                   A[2 + 2 * n] = 3.0;

    printf("A ="); printmat(n, n, A, n); printf("\n");
//work space query
    lwork = -1;
    std::complex<_Float64x> *work = new std::complex<_Float64x>[1];

    Cheev("V", "U", n, A, n, w, work, lwork, rwork, info);
    lwork = (int) cast2double (work[0].real());
    delete[]work;
    work = new std::complex<_Float64x>[std::max((mplapackint) 1, lwork)];
//inverse matrix
    Cheev("V", "U", n, A, n, w, work, lwork, rwork, info);
//print out some results.
    printf("#eigenvalues \n");
    printf("w ="); printmat(n, 1, w, 1); printf("\n");

    printf("#eigenvecs \n");
    printf("U ="); printmat(n, n, A, n); printf("\n");
    printf("#you can check eigenvalues using octave/Matlab by:\n");
    printf("eig(A)\n");
    printf("#you can check eigenvectors using octave/Matlab by:\n");
    printf("U'*A*U\n");

    delete[]work;
    delete[]w;
    delete[]A;
}
