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
#include <mplapack_utils__Float64x.h>

bool rselect(_Float64x ar, _Float64x ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 10;
    std::complex<_Float64x> *a = new std::complex<_Float64x>[n * n];
    std::complex<_Float64x> *w = new std::complex<_Float64x>[n];
    std::complex<_Float64x> *vl = new std::complex<_Float64x>[n * n];
    std::complex<_Float64x> *vr = new std::complex<_Float64x>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<_Float64x> *work = new std::complex<_Float64x>[lwork];    
    _Float64x *rwork = new _Float64x[lwork];
    mplapackint info;
    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a [ (i - 1) + (j - 1) * n ] = 0.0;
        }
    }
    //Tridiagonal Toeplitz matrices: properties and novel applications 
    //https://doi.org/10.1002/nla.1811
    //http://www.math.kent.edu/~reichel/publications/toep3.pdf

    std::complex<_Float64x> sigma = std::complex<_Float64x>(4.0, 3.0) / _Float64x(8.0);
    std::complex<_Float64x> delta = std::complex<_Float64x>(16.0, -3.0);
    std::complex<_Float64x> tau   = std::complex<_Float64x>(0.0, -5.0);

    for (int i = 1; i <= n; i++) {
        a [ (i - 1) + (i - 1) * n ] = delta;
    }

    for (int i = 1; i <= n - 1; i++) {
        a [ (i - 1) + i * n ] = sigma;
        a [ i + (i - 1) * n ] = tau;
    }

    printf("# Tridiagonal Toeplitz matrices: properties and novel applications, https://doi.org/10.1002/nla.1811 http://www.math.kent.edu/~reichel/publications/toep3.pdf\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgeev("V", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
    printf("lambda ="); printvec(w,n); printf("\n");

    std::complex<_Float64x> _pi = pi(_Float64x(0.0));
    std::complex<_Float64x> *lambda = new std::complex<_Float64x>[n];
    for (int h = 1; h <= n; h++) {
        lambda [h - 1] = delta + std::complex<_Float64x>(2.0, 0.0) * sqrt (sigma * tau) * cos( (_Float64x(h) * _pi) / _Float64x((int)n + 1) );
    }
    printf("lambda_true = "); printvec(lambda, n); printf("\n");
    printf("vr ="); printmat(n,n,vr,n); printf("\n");    

    delete[] lambda;
    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
}
