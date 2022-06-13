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
bool cselect(std::complex<_Float64x> a) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    std::complex<_Float64x> *a = new std::complex<_Float64x>[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 2 * n;
    std::complex<_Float64x> *w = new std::complex<_Float64x>[n];
    std::complex<_Float64x> *vs = new std::complex<_Float64x>[n * n];
    std::complex<_Float64x> *work = new std::complex<_Float64x>[lwork];
    _Float64x *rwork = new _Float64x[n];
    bool bwork[n];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = std::complex<_Float64x>(3.0,  0.0); a[0 + 1 * n] = std::complex<_Float64x>(1.0, 0.0);   a[0 + 2 * n] = std::complex<_Float64x>(0.0, 0.0);  a[0 + 3 * n] = std::complex<_Float64x>(0.0, 2.0);
    a[1 + 0 * n] = std::complex<_Float64x>(1.0,  0.0); a[1 + 1 * n] = std::complex<_Float64x>(3.0, 0.0);   a[1 + 2 * n] = std::complex<_Float64x>(0.0, -2.0); a[1 + 3 * n] = std::complex<_Float64x>(0.0, 0.0);
    a[2 + 0 * n] = std::complex<_Float64x>(0.0,  0.0); a[2 + 1 * n] = std::complex<_Float64x>(0.0, 2.0);   a[2 + 2 * n] = std::complex<_Float64x>(1.0, 0.0);  a[2 + 3 * n] = std::complex<_Float64x>(1.0, 0.0);
    a[3 + 0 * n] = std::complex<_Float64x>(0.0, -2.0); a[3 + 1 * n] = std::complex<_Float64x>(0.0, 0.0);   a[3 + 2 * n] = std::complex<_Float64x>(1.0, 0.0);  a[3 + 3 * n] = std::complex<_Float64x>(1.0, 0.0); 

    printf("# Ex. 6.6 p. 116, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgees("V", "S", cselect, n, a, n, sdim, w, vs, n, work, lwork, rwork, bwork, info);
    printf("w ="); printvec(w, n); printf("\n");
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
    printf("eig(a)\n");

    delete[] rwork;
    delete[] work;
    delete[] vs;
    delete[] w;
    delete[] a;
}
