//public domain
#include <mpblas_binary80.h>
#include <mplapack_binary80.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#define BUFLEN 1024
void printnum(mplapack_binary80_t rtmp)
{
    int width = 25;
    char buf[BUFLEN];
#if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_STRFROMF64X
    strfromf64x(buf, sizeof(buf), "%.21e", rtmp);
#elif MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL
    snprintf(buf, sizeof(buf), "%*.21Le", width, rtmp);
#else
    #error "unsupported binary80 type"
#endif
    if (rtmp >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    return;
}
void printnum(std::complex<mplapack_binary80_t> ctmp)
{
    int width = 25;
    char buf[BUFLEN];
#if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_STRFROMF64X
    strfromf64x(buf, sizeof(buf), "%.21e", ctmp.real());
    if (ctmp.real() >= 0.0) printf("+%s", buf); else printf("%s", buf);
    strfromf64x(buf, sizeof(buf), "%.21e", ctmp.imag());
    if (ctmp.imag() >= 0.0) printf("+%s", buf); else printf("%s", buf);
#elif MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL
    snprintf(buf, sizeof(buf), "%*.21Le", width, ctmp.real());
    if (ctmp.real() >= 0.0) printf("+%s", buf); else printf("%s", buf);
    snprintf(buf, sizeof(buf), "%*.21Le", width, ctmp.imag());
    if (ctmp.imag() >= 0.0) printf("+%s", buf); else printf("%s", buf);
#else
    #error "unsupported binary80 type"
#endif
    printf("i");
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
//taking from Collection of Matrices for Testing Computational Algorithms 1969 Robert T. Gregory, David L. Karney pp.30
int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    std::complex<mplapack_binary80_t> *a = new std::complex<mplapack_binary80_t>[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix


    a[0 + 0 * n] = std::complex<mplapack_binary80_t>(1.0, 0.0);   a[0 + 1 * n] = std::complex<mplapack_binary80_t>(1.0, 2.0);    a[0 + 2 * n] = std::complex<mplapack_binary80_t>(2.0, 10.0);
    a[1 + 0 * n] = std::complex<mplapack_binary80_t>(1.0, 1.0);   a[1 + 1 * n] = std::complex<mplapack_binary80_t>(0.0, 3.0);    a[1 + 2 * n] = std::complex<mplapack_binary80_t>(-5.0, 14.0);
    a[2 + 0 * n] = std::complex<mplapack_binary80_t>(1.0, 1.0);   a[2 + 1 * n] = std::complex<mplapack_binary80_t>(0.0, 5.0);    a[2 + 2 * n] = std::complex<mplapack_binary80_t>(-8.0, 20.0);

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    std::complex<mplapack_binary80_t> *work = new std::complex<mplapack_binary80_t>[1];

    Cgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_binary80 (work[0].real());
    delete[]work;
    work = new std::complex<mplapack_binary80_t>[std::max(1, (int) lwork)];

//inverse matrix
    Cgetrf(n, n, a, n, ipiv, info);
    Cgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("ainv * a - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
