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
bool cselect(std::complex<double> a) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    std::complex<double> *a = new std::complex<double>[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 2 * n;
    std::complex<double> *w = new std::complex<double>[n];
    std::complex<double> *vs = new std::complex<double>[n * n];
    std::complex<double> *work = new std::complex<double>[lwork];
    double *rwork = new double[n];
    bool bwork[n];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = std::complex<double>(3.0,  0.0); a[0 + 1 * n] = std::complex<double>(1.0, 0.0);   a[0 + 2 * n] = std::complex<double>(0.0, 0.0);  a[0 + 3 * n] = std::complex<double>(0.0, 2.0);
    a[1 + 0 * n] = std::complex<double>(1.0,  0.0); a[1 + 1 * n] = std::complex<double>(3.0, 0.0);   a[1 + 2 * n] = std::complex<double>(0.0, -2.0); a[1 + 3 * n] = std::complex<double>(0.0, 0.0);
    a[2 + 0 * n] = std::complex<double>(0.0,  0.0); a[2 + 1 * n] = std::complex<double>(0.0, 2.0);   a[2 + 2 * n] = std::complex<double>(1.0, 0.0);  a[2 + 3 * n] = std::complex<double>(1.0, 0.0);
    a[3 + 0 * n] = std::complex<double>(0.0, -2.0); a[3 + 1 * n] = std::complex<double>(0.0, 0.0);   a[3 + 2 * n] = std::complex<double>(1.0, 0.0);  a[3 + 3 * n] = std::complex<double>(1.0, 0.0); 

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
