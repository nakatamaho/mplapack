//public domain
#include <mpblas__Float128.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(_Float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

void printnum(std::complex<_Float128> rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp.real());
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp.real());
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.real());
#endif
    if (rtmp.real() >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp.imag());
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp.imag());
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.imag());
#endif
    if (rtmp.imag() >= 0.0)
        printf ("+%si", buf);
    else
        printf ("%si", buf);
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
int main()
{
    mplapackint n = 3;

    std::complex<_Float128> *a = new std::complex<_Float128>[n * n];
    std::complex<_Float128> *b = new std::complex<_Float128>[n * n];
    std::complex<_Float128> *c = new std::complex<_Float128>[n * n];
    std::complex<_Float128> alpha, beta;

//setting A matrix
    a[0 + 0 * n] = std::complex<_Float128>(1.0,-1.0);    a[0 + 1 * n] = std::complex<_Float128>(8.0, 2.2);    a[0 + 2 * n] = std::complex<_Float128>(0.0, -10.0);
    a[1 + 0 * n] = std::complex<_Float128>(2.0, 0.0);    a[1 + 1 * n] = std::complex<_Float128>(10.0,0.0);    a[1 + 2 * n] = std::complex<_Float128>(8.1, 2.2);
    a[2 + 0 * n] = std::complex<_Float128>(-9.0,3.0);    a[2 + 1 * n] = std::complex<_Float128>(-5.0,3.0);    a[2 + 2 * n] = std::complex<_Float128>(-1.0, 0.0);

    b[0 + 0 * n] = std::complex<_Float128>(9.0, 0.0);    b[0 + 1 * n] = std::complex<_Float128>(8.0, -0.01);  b[0 + 2 * n] = std::complex<_Float128>(3.0, 1.001);
    b[1 + 0 * n] = std::complex<_Float128>(3.0, -8.0);   b[1 + 1 * n] = std::complex<_Float128>(-11.0, 0.1);  b[1 + 2 * n] = std::complex<_Float128>(8.0, 0.00001);
    b[2 + 0 * n] = std::complex<_Float128>(-8.0, 1.0);   b[2 + 1 * n] = std::complex<_Float128>(6.0, 0.0);    b[2 + 2 * n] = std::complex<_Float128>(1.1, 1.0);

    c[0 + 0 * n] = std::complex<_Float128>(3.0, 1.0);   c[0 + 1 * n] = std::complex<_Float128>(-3.0, 9.99);   c[0 + 2 * n] = std::complex<_Float128>(-9.0, -11.0);
    c[1 + 0 * n] = std::complex<_Float128>(8.0, -1.0);  c[1 + 1 * n] = std::complex<_Float128>(4.0, 4.44);    c[1 + 2 * n] = std::complex<_Float128>(8.0, 9.0);
    c[2 + 0 * n] = std::complex<_Float128>(6.0, 0.0);   c[2 + 1 * n] = std::complex<_Float128>(-1.0, 0.0);    c[2 + 2 * n] = std::complex<_Float128>(-2.0, 1.0);

    printf("# Cgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = std::complex<_Float128>(3.0,-1.2);
    beta = std::complex<_Float128>(-2.0, -2.0);
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
