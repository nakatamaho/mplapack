//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(mplapack_binary128_t rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

void printnum(std::complex<mplapack_binary128_t> rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.real());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp.real());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp.real());
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp.real() >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.imag());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp.imag());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp.imag());
#else
    #error "unsupported binary128 type"
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
bool rselect(mplapack_binary128_t ar, mplapack_binary128_t ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    std::complex<mplapack_binary128_t> *a = new std::complex<mplapack_binary128_t>[n * n];
    std::complex<mplapack_binary128_t> *w = new std::complex<mplapack_binary128_t>[n];
    std::complex<mplapack_binary128_t> *vl = new std::complex<mplapack_binary128_t>[n * n];
    std::complex<mplapack_binary128_t> *vr = new std::complex<mplapack_binary128_t>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<mplapack_binary128_t> *work = new std::complex<mplapack_binary128_t>[lwork];    
    mplapack_binary128_t *rwork = new mplapack_binary128_t[lwork];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = std::complex<mplapack_binary128_t>(5.0, 9.0); a[0 + 1 * n] = std::complex<mplapack_binary128_t>(5.0, 5.0);   a[0 + 2 * n] = std::complex<mplapack_binary128_t>(-6.0, -6.0); a[0 + 3 * n] = std::complex<mplapack_binary128_t>(-7.0, -7.0);
    a[1 + 0 * n] = std::complex<mplapack_binary128_t>(3.0, 3.0); a[1 + 1 * n] = std::complex<mplapack_binary128_t>(6.0, 10.0);  a[1 + 2 * n] = std::complex<mplapack_binary128_t>(-5.0, -5.0); a[1 + 3 * n] = std::complex<mplapack_binary128_t>(-6.0, -6.0);
    a[2 + 0 * n] = std::complex<mplapack_binary128_t>(2.0, 2.0); a[2 + 1 * n] = std::complex<mplapack_binary128_t>(3.0, 3.0);   a[2 + 2 * n] = std::complex<mplapack_binary128_t>(-1.0,  3.0); a[2 + 3 * n] = std::complex<mplapack_binary128_t>(-5.0, -5.0);
    a[3 + 0 * n] = std::complex<mplapack_binary128_t>(1.0, 1.0); a[3 + 1 * n] = std::complex<mplapack_binary128_t>(2.0, 2.0);   a[3 + 2 * n] = std::complex<mplapack_binary128_t>(-3.0, -3.0); a[3 + 3 * n] = std::complex<mplapack_binary128_t>(0.0, 4.0); 

    printf("# Ex. 6.5 p. 116, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
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
