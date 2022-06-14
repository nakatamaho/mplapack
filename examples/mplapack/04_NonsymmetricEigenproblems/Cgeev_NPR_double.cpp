//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }
inline void printnum(std::complex<double> ctmp) { printf(DOUBLE_FORMAT DOUBLE_FORMAT "i", ctmp.real(), ctmp.imag()); }

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
#include <mplapack_utils_double.h>

bool rselect(double ar, double ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 10;
    std::complex<double> *a = new std::complex<double>[n * n];
    std::complex<double> *w = new std::complex<double>[n];
    std::complex<double> *vl = new std::complex<double>[n * n];
    std::complex<double> *vr = new std::complex<double>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<double> *work = new std::complex<double>[lwork];    
    double *rwork = new double[lwork];
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

    std::complex<double> sigma = std::complex<double>(4.0, 3.0) / double(8.0);
    std::complex<double> delta = std::complex<double>(16.0, -3.0);
    std::complex<double> tau   = std::complex<double>(0.0, -5.0);

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

    std::complex<double> _pi = pi(double(0.0));
    std::complex<double> *lambda = new std::complex<double>[n];
    for (int h = 1; h <= n; h++) {
        lambda [h - 1] = delta + std::complex<double>(2.0, 0.0) * sqrt (sigma * tau) * cos( (double(h) * _pi) / double((int)n + 1) );
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
