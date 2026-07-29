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
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, info, lwork = -1;
    std::complex<mplapack_binary128_t> *a = new std::complex<mplapack_binary128_t>[n * n];
    std::complex<mplapack_binary128_t> *b = new std::complex<mplapack_binary128_t>[n * n];
    std::complex<mplapack_binary128_t> *alpha = new std::complex<mplapack_binary128_t>[n];
    std::complex<mplapack_binary128_t> *beta = new std::complex<mplapack_binary128_t>[n];
    std::complex<mplapack_binary128_t> *vl = new std::complex<mplapack_binary128_t>[1];
    std::complex<mplapack_binary128_t> *vr = new std::complex<mplapack_binary128_t>[1];
    mplapack_binary128_t *rwork = new mplapack_binary128_t[8 * n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = std::complex<mplapack_binary128_t>(0.0, 0.0);
        b[i] = std::complex<mplapack_binary128_t>(0.0, 0.0);
    }
    a[0] = std::complex<mplapack_binary128_t>(1.0, 1.0);
    a[4] = std::complex<mplapack_binary128_t>(2.0, 0.0);
    a[8] = std::complex<mplapack_binary128_t>(3.0, -1.0);
    b[0] = std::complex<mplapack_binary128_t>(1.0, 0.0);
    b[4] = std::complex<mplapack_binary128_t>(1.0, 0.0);
    b[8] = std::complex<mplapack_binary128_t>(0.0, 0.0);
    std::complex<mplapack_binary128_t> wk;
    Cggev("N", "N", n, a, lda, b, ldb, alpha, beta, vl, ldv, vr, ldv, &wk, lwork, rwork, info);
    lwork = castINTEGER_binary128(wk.real());
    std::complex<mplapack_binary128_t> *work = new std::complex<mplapack_binary128_t>[lwork];
    Cggev("N", "N", n, a, lda, b, ldb, alpha, beta, vl, ldv, vr, ldv, work, lwork, rwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alpha[i]); printf(", beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_binary128("E"))
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
