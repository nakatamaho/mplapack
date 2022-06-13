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
int main()
{
    mplapackint n = 3;

    std::complex<double> *a = new std::complex<double>[n * n];
    std::complex<double> *b = new std::complex<double>[n * n];
    std::complex<double> *c = new std::complex<double>[n * n];
    std::complex<double> alpha, beta;

//setting A matrix
    a[0 + 0 * n] = std::complex<double>(1.0,-1.0);    a[0 + 1 * n] = std::complex<double>(8.0, 2.2);    a[0 + 2 * n] = std::complex<double>(0.0, -10.0);
    a[1 + 0 * n] = std::complex<double>(2.0, 0.0);    a[1 + 1 * n] = std::complex<double>(10.0,0.0);    a[1 + 2 * n] = std::complex<double>(8.1, 2.2);
    a[2 + 0 * n] = std::complex<double>(-9.0,3.0);    a[2 + 1 * n] = std::complex<double>(-5.0,3.0);    a[2 + 2 * n] = std::complex<double>(-1.0, 0.0);

    b[0 + 0 * n] = std::complex<double>(9.0, 0.0);    b[0 + 1 * n] = std::complex<double>(8.0, -0.01);  b[0 + 2 * n] = std::complex<double>(3.0, 1.001);
    b[1 + 0 * n] = std::complex<double>(3.0, -8.0);   b[1 + 1 * n] = std::complex<double>(-11.0, 0.1);  b[1 + 2 * n] = std::complex<double>(8.0, 0.00001);
    b[2 + 0 * n] = std::complex<double>(-8.0, 1.0);   b[2 + 1 * n] = std::complex<double>(6.0, 0.0);    b[2 + 2 * n] = std::complex<double>(1.1, 1.0);

    c[0 + 0 * n] = std::complex<double>(3.0, 1.0);   c[0 + 1 * n] = std::complex<double>(-3.0, 9.99);   c[0 + 2 * n] = std::complex<double>(-9.0, -11.0);
    c[1 + 0 * n] = std::complex<double>(8.0, -1.0);  c[1 + 1 * n] = std::complex<double>(4.0, 4.44);    c[1 + 2 * n] = std::complex<double>(8.0, 9.0);
    c[2 + 0 * n] = std::complex<double>(6.0, 0.0);   c[2 + 1 * n] = std::complex<double>(-1.0, 0.0);    c[2 + 2 * n] = std::complex<double>(-2.0, 1.0);

    printf("# Cgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = std::complex<double>(3.0,-1.2);
    beta = std::complex<double>(-2.0, -2.0);
    Cgemm("n", "n", n, n, n, alpha, a, n, b, n, beta, c, n);

    printf("alpha = "); printnum(alpha); printf("\n");
    printf("beta = "); printnum(beta); printf("\n");
    printf("ans ="); printmat(n, n, c, n); printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * a * b + beta * c \n");
    delete[]c;
    delete[]b;
    delete[]a;
}
