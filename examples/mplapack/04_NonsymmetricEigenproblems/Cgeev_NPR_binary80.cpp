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
#include <mplapack_utils_binary80.h>

bool rselect(mplapack_binary80_t ar, mplapack_binary80_t ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 10;
    std::complex<mplapack_binary80_t> *a = new std::complex<mplapack_binary80_t>[n * n];
    std::complex<mplapack_binary80_t> *w = new std::complex<mplapack_binary80_t>[n];
    std::complex<mplapack_binary80_t> *vl = new std::complex<mplapack_binary80_t>[n * n];
    std::complex<mplapack_binary80_t> *vr = new std::complex<mplapack_binary80_t>[n * n];
    mplapackint lwork = 4 * n;
    std::complex<mplapack_binary80_t> *work = new std::complex<mplapack_binary80_t>[lwork];    
    mplapack_binary80_t *rwork = new mplapack_binary80_t[lwork];
    mplapackint info;
    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a [ (i - 1) + (j - 1) * n ] = 0.0;
        }
    }
    //Tridiagonal Toeplitz matrices: properties and novel applications 
    //https://doi.org/10.1002/nla.1811
    //http://www.math.kent.edu/~reichel/publications/toep3.pdf

    std::complex<mplapack_binary80_t> sigma = std::complex<mplapack_binary80_t>(4.0, 3.0) / mplapack_binary80_t(8.0);
    std::complex<mplapack_binary80_t> delta = std::complex<mplapack_binary80_t>(16.0, -3.0);
    std::complex<mplapack_binary80_t> tau   = std::complex<mplapack_binary80_t>(0.0, -5.0);

    for (int i = 1; i <= n; i++) {
        a [ (i - 1) + (i - 1) * n ] = delta;
    }

    for (int i = 1; i <= n - 1; i++) {
        a [ (i - 1) + i * n ] = sigma;
        a [ i + (i - 1) * n ] = tau;
    }

    printf("# Tridiagonal Toeplitz matrices: properties and novel applications, https://doi.org/10.1002/nla.1811 http://www.math.kent.edu/~reichel/publications/toep3.pdf\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgeev("V", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
    printf("lambda ="); printvec(w,n); printf("\n");

    std::complex<mplapack_binary80_t> _pi = pi(mplapack_binary80_t(0.0));
    std::complex<mplapack_binary80_t> *lambda = new std::complex<mplapack_binary80_t>[n];
    for (int h = 1; h <= n; h++) {
        lambda [h - 1] = delta + std::complex<mplapack_binary80_t>(2.0, 0.0) * sqrt (sigma * tau) * cos( (mplapack_binary80_t((double)h) * _pi) / mplapack_binary80_t((double)((int)n + 1)) );
    }
    printf("lambda_true = "); printvec(lambda, n); printf("\n");
    printf("vr ="); printmat(n,n,vr,n); printf("\n");    

    delete[] lambda;
    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
}
