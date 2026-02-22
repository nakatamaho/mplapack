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
    mplapackint n = 4;
    mplapackint m = 4;

    std::complex<mplapack_binary128_t> *a = new std::complex<mplapack_binary128_t>[m * n];
    mplapack_binary128_t *s = new mplapack_binary128_t[std::min(m, n)];
    std::complex<mplapack_binary128_t> *u = new std::complex<mplapack_binary128_t>[m * m];
    std::complex<mplapack_binary128_t> *vt = new std::complex<mplapack_binary128_t>[n * n];
    mplapackint lwork = std::max((mplapackint)1, 2 * std::min(m, n) + std::max(m, n));
    std::complex<mplapack_binary128_t> *work = new std::complex<mplapack_binary128_t>[lwork];
    mplapack_binary128_t *rwork = new mplapack_binary128_t[5 * std::min(m, n)];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = std::complex<mplapack_binary128_t>(0.9, -1.0); a[0 + 1 * n] = std::complex<mplapack_binary128_t>(20.0, -2.25);  a[0 + 2 * n] = std::complex<mplapack_binary128_t>(1.75, -0.5);  a[0 + 3 * n] = std::complex<mplapack_binary128_t>(0.0, 0.5);
    a[1 + 0 * n] = std::complex<mplapack_binary128_t>(8.0,-2.25); a[1 + 1 * n] = std::complex<mplapack_binary128_t>(-0.25, 0.0);   a[1 + 2 * n] = std::complex<mplapack_binary128_t>(1.25, -0.25); a[1 + 3 * n] = std::complex<mplapack_binary128_t>(-3.75, 0.0);
    a[2 + 0 * n] = std::complex<mplapack_binary128_t>(-1.75,0.0); a[2 + 1 * n] = std::complex<mplapack_binary128_t>(-80.0,  1.25); a[2 + 2 * n] = std::complex<mplapack_binary128_t>(1.5, 0.0);    a[2 + 3 * n] = std::complex<mplapack_binary128_t>(30.0, 2.25);
    a[3 + 0 * n] = std::complex<mplapack_binary128_t>(3.0, 0.25); a[3 + 1 * n] = std::complex<mplapack_binary128_t>(1.75, 0.0);    a[3 + 2 * n] = std::complex<mplapack_binary128_t>(0.0, 2.25);   a[3 + 3 * n] = std::complex<mplapack_binary128_t>(-0.25, -80.0);
    
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
