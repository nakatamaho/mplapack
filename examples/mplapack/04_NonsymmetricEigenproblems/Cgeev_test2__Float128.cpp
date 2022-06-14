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
bool rselect(_Float128 ar, _Float128 ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    std::complex<_Float128> *a = new std::complex<_Float128>[n * n];
    std::complex<_Float128> *w = new std::complex<_Float128>[n];
    std::complex<_Float128> *vl = new std::complex<_Float128>[n * n];
    std::complex<_Float128> *vr = new std::complex<_Float128>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<_Float128> *work = new std::complex<_Float128>[lwork];    
    _Float128 *rwork = new _Float128[lwork];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = std::complex<_Float128>(7.0, 0.0);   a[0 + 1 * n] = std::complex<_Float128>(3.0, 0.0);  a[0 + 2 * n] = std::complex<_Float128>(1.0, 2.0);   a[0 + 3 * n] = std::complex<_Float128>(-1.0, 2.0);
    a[1 + 0 * n] = std::complex<_Float128>(3.0, 0.0);   a[1 + 1 * n] = std::complex<_Float128>(7.0, 0.0);  a[1 + 2 * n] = std::complex<_Float128>(1.0, -2.0);  a[1 + 3 * n] = std::complex<_Float128>(-1.0, -2.0);
    a[2 + 0 * n] = std::complex<_Float128>(1.0, -2.0);  a[2 + 1 * n] = std::complex<_Float128>(1.0, 2.0);  a[2 + 2 * n] = std::complex<_Float128>(7.0, 0.0);   a[2 + 3 * n] = std::complex<_Float128>(-3.0, 0.0);
    a[3 + 0 * n] = std::complex<_Float128>(-1.0, -2.0); a[3 + 1 * n] = std::complex<_Float128>(-1.0, 2.0); a[3 + 2 * n] = std::complex<_Float128>(-3.0, 0.0);  a[3 + 3 * n] = std::complex<_Float128>(7.0, 0.0); 

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
