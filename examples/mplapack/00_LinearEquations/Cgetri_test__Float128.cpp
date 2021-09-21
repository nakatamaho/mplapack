//public domain
#include <mpblas__Float128.h>
#include <mplapack__Float128.h>
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
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp);
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
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp.real());
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
    n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp.imag());
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp.imag());
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.imag());
#endif
    if (rtmp.imag() >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
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
    mplapackint lwork, info;

    std::complex<_Float128> *a = new std::complex<_Float128>[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix


    a[0 + 0 * n] = std::complex<_Float128>(1.0, 0.0);   a[0 + 1 * n] = std::complex<_Float128>(1.0, 2.0);    a[0 + 2 * n] = std::complex<_Float128>(2.0, 10.0);
    a[1 + 0 * n] = std::complex<_Float128>(1.0, 1.0);   a[1 + 1 * n] = std::complex<_Float128>(0.0, 3.0);    a[1 + 2 * n] = std::complex<_Float128>(-5.0, 14.0);
    a[2 + 0 * n] = std::complex<_Float128>(1.0, 1.0);   a[2 + 1 * n] = std::complex<_Float128>(0.0, 5.0);    a[2 + 2 * n] = std::complex<_Float128>(-8.0, 20.0);

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    std::complex<_Float128> *work = new std::complex<_Float128>[1];

    Cgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER__Float128 (work[0].real());
    delete[]work;
    work = new std::complex<_Float128>[std::max(1, (int) lwork)];

//inverse matrix
    Cgetrf(n, n, a, n, ipiv, info);
    Cgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
