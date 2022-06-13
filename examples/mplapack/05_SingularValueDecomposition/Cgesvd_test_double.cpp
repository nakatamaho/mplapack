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
    mplapackint n = 4;
    mplapackint m = 4;

    std::complex<double> *a = new std::complex<double>[m * n];
    double *s = new double[std::min(m, n)];
    std::complex<double> *u = new std::complex<double>[m * m];
    std::complex<double> *vt = new std::complex<double>[n * n];
    mplapackint lwork = std::max((mplapackint)1, 2 * std::min(m, n) + std::max(m, n));
    std::complex<double> *work = new std::complex<double>[lwork];
    double *rwork = new double[5 * std::min(m, n)];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = std::complex<double>(0.25, 0.0); a[0 + 1 * n] = std::complex<double>(0.0, -2.25); a[0 + 2 * n] = std::complex<double>(-1.75,0.0);   a[0 + 3 * n] = std::complex<double>(0.0, 0.25);
    a[1 + 0 * n] = std::complex<double>(0.0,-2.25); a[1 + 1 * n] = std::complex<double>(-0.25, 0.0); a[1 + 2 * n] = std::complex<double>(0.0, -0.25);  a[1 + 3 * n] = std::complex<double>(-1.75, 0.0);
    a[2 + 0 * n] = std::complex<double>(-1.75,0.0); a[2 + 1 * n] = std::complex<double>(0.0, -0.25); a[2 + 2 * n] = std::complex<double>(0.25, 0.0);   a[2 + 3 * n] = std::complex<double>(0.0, 2.25);
    a[3 + 0 * n] = std::complex<double>(0.0, 0.25); a[3 + 1 * n] = std::complex<double>(-1.75, 0.0); a[3 + 2 * n] = std::complex<double>(0.0, 2.25);   a[3 + 3 * n] = std::complex<double>(-0.25, 0.0);
    
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Cgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, rwork, info);
    printf("s="); printvec(s, std::min(m, n)); printf("\n");
    if (m < n)
        printf("padding=zeros(%d, %d-%d)\n", (int)m, (int)n, (int)m);
    if (n < m)
        printf("padding=zeros(%d-%d,%d)\n", (int)m, (int)n, (int)n);
    printf("u ="); printmat(m, m, u, m); printf("\n");
    printf("vt ="); printmat(n, n, vt, n); printf("\n");
    printf("svd(a)\n");
    if (m < n)
        printf("sigma=[diag(s) padding] \n");
    if (n < m)
        printf("sigma=[diag(s); padding] \n");
    if (n == m)
        printf("sigma=[diag(s)] \n");
    printf("sigma \n");
    printf("u * sigma  * vt\n");
    delete[] rwork;
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
