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
bool rselect(double ar, double ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    std::complex<double> *a = new std::complex<double>[n * n];
    std::complex<double> *w = new std::complex<double>[n];
    std::complex<double> *vl = new std::complex<double>[n * n];
    std::complex<double> *vr = new std::complex<double>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<double> *work = new std::complex<double>[lwork];    
    double *rwork = new double[lwork];
    mplapackint info;
    // setting A matrix
    //# Example 6.5 "Collection of Matrices for Testing Computational Algorithms", Robert T. Gregory, David L. Karney    
    a[0 + 0 * n] = std::complex<double>(7.0, 0.0);   a[0 + 1 * n] = std::complex<double>(3.0, 0.0);  a[0 + 2 * n] = std::complex<double>(1.0, 2.0);   a[0 + 3 * n] = std::complex<double>(-1.0, 2.0);
    a[1 + 0 * n] = std::complex<double>(3.0, 0.0);   a[1 + 1 * n] = std::complex<double>(7.0, 0.0);  a[1 + 2 * n] = std::complex<double>(1.0, -2.0);  a[1 + 3 * n] = std::complex<double>(-1.0, -2.0);
    a[2 + 0 * n] = std::complex<double>(1.0, -2.0);  a[2 + 1 * n] = std::complex<double>(1.0, 2.0);  a[2 + 2 * n] = std::complex<double>(7.0, 0.0);   a[2 + 3 * n] = std::complex<double>(-3.0, 0.0);
    a[3 + 0 * n] = std::complex<double>(-1.0, -2.0); a[3 + 1 * n] = std::complex<double>(-1.0, 2.0); a[3 + 2 * n] = std::complex<double>(-3.0, 0.0);  a[3 + 3 * n] = std::complex<double>(7.0, 0.0); 

    printf("# Ex. 6.7 p. 117, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgeev("V", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
    printf("lambda ="); printvec(w,n); printf("\n");    
    printf("vr ="); printmat(n,n,vr,n); printf("\n");    

    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
}
