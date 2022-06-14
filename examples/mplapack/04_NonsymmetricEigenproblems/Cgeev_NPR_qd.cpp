//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum(qd_complex rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}

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
#include <mplapack_utils_qd.h>

bool rselect(qd_real ar, qd_real ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 10;
    qd_complex *a = new qd_complex[n * n];
    qd_complex *w = new qd_complex[n];
    qd_complex *vl = new qd_complex[n * n];
    qd_complex *vr = new qd_complex[n * n];
    mplapackint lwork = 4 * n;
    qd_complex *work = new qd_complex[lwork];    
    qd_real *rwork = new qd_real[lwork];
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

    qd_complex sigma = qd_complex(4.0, 3.0) / qd_real(8.0);
    qd_complex delta = qd_complex(16.0, -3.0);
    qd_complex tau   = qd_complex(0.0, -5.0);

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

    qd_complex _pi = pi(qd_real(0.0));
    qd_complex *lambda = new qd_complex[n];
    for (int h = 1; h <= n; h++) {
        lambda [h - 1] = delta + qd_complex(2.0, 0.0) * sqrt (sigma * tau) * cos( (qd_real(h) * _pi) / qd_real((int)n + 1) );
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
