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
int main() {
    mplapackint n = 4;
    mplapackint m = 4;

    std::complex<_Float128> *a = new std::complex<_Float128>[m * n];
    _Float128 *s = new _Float128[std::min(m, n)];
    std::complex<_Float128> *u = new std::complex<_Float128>[m * m];
    std::complex<_Float128> *vt = new std::complex<_Float128>[n * n];
    mplapackint lwork = std::max((mplapackint)1, 2 * std::min(m, n) + std::max(m, n));
    std::complex<_Float128> *work = new std::complex<_Float128>[lwork];
    _Float128 *rwork = new _Float128[5 * std::min(m, n)];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = std::complex<_Float128>(0.25, 0.0); a[0 + 1 * n] = std::complex<_Float128>(0.0, -2.25); a[0 + 2 * n] = std::complex<_Float128>(-1.75,0.0);   a[0 + 3 * n] = std::complex<_Float128>(0.0, 0.25);
    a[1 + 0 * n] = std::complex<_Float128>(0.0,-2.25); a[1 + 1 * n] = std::complex<_Float128>(-0.25, 0.0); a[1 + 2 * n] = std::complex<_Float128>(0.0, -0.25);  a[1 + 3 * n] = std::complex<_Float128>(-1.75, 0.0);
    a[2 + 0 * n] = std::complex<_Float128>(-1.75,0.0); a[2 + 1 * n] = std::complex<_Float128>(0.0, -0.25); a[2 + 2 * n] = std::complex<_Float128>(0.25, 0.0);   a[2 + 3 * n] = std::complex<_Float128>(0.0, 2.25);
    a[3 + 0 * n] = std::complex<_Float128>(0.0, 0.25); a[3 + 1 * n] = std::complex<_Float128>(-1.75, 0.0); a[3 + 2 * n] = std::complex<_Float128>(0.0, 2.25);   a[3 + 3 * n] = std::complex<_Float128>(-0.25, 0.0);
    
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
