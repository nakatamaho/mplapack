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
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, info, lwork = -1;
    std::complex<double> *a = new std::complex<double>[n * n];
    std::complex<double> *b = new std::complex<double>[n * n];
    std::complex<double> *alpha = new std::complex<double>[n];
    std::complex<double> *beta = new std::complex<double>[n];
    std::complex<double> *vl = new std::complex<double>[1];
    std::complex<double> *vr = new std::complex<double>[1];
    double *rwork = new double[8 * n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = std::complex<double>(0.0, 0.0);
        b[i] = std::complex<double>(0.0, 0.0);
    }
    a[0] = std::complex<double>(1.0, 1.0);
    a[4] = std::complex<double>(2.0, 0.0);
    a[8] = std::complex<double>(3.0, -1.0);
    b[0] = std::complex<double>(1.0, 0.0);
    b[4] = std::complex<double>(1.0, 0.0);
    b[8] = std::complex<double>(0.0, 0.0);
    std::complex<double> wk;
    Cggev("N", "N", n, a, lda, b, ldb, alpha, beta, vl, ldv, vr, ldv, &wk, lwork, rwork, info);
    lwork = castINTEGER_double(wk.real());
    std::complex<double> *work = new std::complex<double>[lwork];
    Cggev("N", "N", n, a, lda, b, ldb, alpha, beta, vl, ldv, vr, ldv, work, lwork, rwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alpha[i]); printf(", beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_double("E"))
            printf(", lambda = Inf\n");
        else {
            printf(", lambda = "); printnum(alpha[i] / beta[i]); printf("\n");
        }
    }
    delete[] work;
    delete[] rwork;
    delete[] vr;
    delete[] vl;
    delete[] beta;
    delete[] alpha;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
