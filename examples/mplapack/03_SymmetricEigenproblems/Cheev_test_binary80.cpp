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
int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    std::complex<mplapack_binary80_t> *A = new std::complex<mplapack_binary80_t>[n * n];
    mplapack_binary80_t *w = new mplapack_binary80_t[n];
    mplapack_binary80_t *rwork = new mplapack_binary80_t[3 * n - 1];

//setting A matrix
    A[0 + 0 * n] = 2.0;               A[0 + 1 * n] = std::complex<mplapack_binary80_t>(0.0, -1.0);    A[0 + 2 * n] = 0.0;
    A[1 + 0 * n] = std::complex<mplapack_binary80_t>(0.0, 1.0); A[1 + 1 * n] = 2.0;                   A[1 + 2 * n] = 0.0;
    A[2 + 0 * n] = 0.0;               A[2 + 1 * n] = 0.0;                   A[2 + 2 * n] = 3.0;

    printf("A ="); printmat(n, n, A, n); printf("\n");
//work space query
    lwork = -1;
    std::complex<mplapack_binary80_t> *work = new std::complex<mplapack_binary80_t>[1];

    Cheev("V", "U", n, A, n, w, work, lwork, rwork, info);
    lwork = (int) cast2double (work[0].real());
    delete[]work;
    work = new std::complex<mplapack_binary80_t>[std::max((mplapackint) 1, lwork)];
//inverse matrix
    Cheev("V", "U", n, A, n, w, work, lwork, rwork, info);
//print out some results.
    printf("#eigenvalues \n");
    printf("w ="); printmat(n, 1, w, 1); printf("\n");

    printf("#eigenvecs \n");
    printf("U ="); printmat(n, n, A, n); printf("\n");
    printf("#you can check eigenvalues using octave/Matlab by:\n");
    printf("eig(A)\n");
    printf("#you can check eigenvectors using octave/Matlab by:\n");
    printf("U'*A*U\n");

    delete[]work;
    delete[]w;
    delete[]A;
}
