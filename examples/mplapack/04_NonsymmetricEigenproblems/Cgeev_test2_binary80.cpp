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
bool rselect(mplapack_binary80_t ar, mplapack_binary80_t ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    std::complex<mplapack_binary80_t> *a = new std::complex<mplapack_binary80_t>[n * n];
    std::complex<mplapack_binary80_t> *w = new std::complex<mplapack_binary80_t>[n];
    std::complex<mplapack_binary80_t> *vl = new std::complex<mplapack_binary80_t>[n * n];
    std::complex<mplapack_binary80_t> *vr = new std::complex<mplapack_binary80_t>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<mplapack_binary80_t> *work = new std::complex<mplapack_binary80_t>[lwork];    
    mplapack_binary80_t *rwork = new mplapack_binary80_t[lwork];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = std::complex<mplapack_binary80_t>(7.0, 0.0);   a[0 + 1 * n] = std::complex<mplapack_binary80_t>(3.0, 0.0);  a[0 + 2 * n] = std::complex<mplapack_binary80_t>(1.0, 2.0);   a[0 + 3 * n] = std::complex<mplapack_binary80_t>(-1.0, 2.0);
    a[1 + 0 * n] = std::complex<mplapack_binary80_t>(3.0, 0.0);   a[1 + 1 * n] = std::complex<mplapack_binary80_t>(7.0, 0.0);  a[1 + 2 * n] = std::complex<mplapack_binary80_t>(1.0, -2.0);  a[1 + 3 * n] = std::complex<mplapack_binary80_t>(-1.0, -2.0);
    a[2 + 0 * n] = std::complex<mplapack_binary80_t>(1.0, -2.0);  a[2 + 1 * n] = std::complex<mplapack_binary80_t>(1.0, 2.0);  a[2 + 2 * n] = std::complex<mplapack_binary80_t>(7.0, 0.0);   a[2 + 3 * n] = std::complex<mplapack_binary80_t>(-3.0, 0.0);
    a[3 + 0 * n] = std::complex<mplapack_binary80_t>(-1.0, -2.0); a[3 + 1 * n] = std::complex<mplapack_binary80_t>(-1.0, 2.0); a[3 + 2 * n] = std::complex<mplapack_binary80_t>(-3.0, 0.0);  a[3 + 3 * n] = std::complex<mplapack_binary80_t>(7.0, 0.0); 

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
