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
mplapack_binary80_t max_solution_error(mplapackint n, mplapackint nrhs, std::complex<mplapack_binary80_t> *x, mplapackint ldx, std::complex<mplapack_binary80_t> *xexact, mplapackint ldxexact) {
    mplapack_binary80_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mplapack_binary80_t d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mplapack_binary80_t max_residual(mplapackint m, mplapackint n, mplapackint nrhs, std::complex<mplapack_binary80_t> *a, mplapackint lda, std::complex<mplapack_binary80_t> *x, mplapackint ldx, std::complex<mplapack_binary80_t> *b, mplapackint ldb) {
    mplapack_binary80_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            std::complex<mplapack_binary80_t> s = std::complex<mplapack_binary80_t>(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mplapack_binary80_t d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 2, nrhs = 2, lda = n, ldb = n, info;
    std::complex<mplapack_binary80_t> *a = new std::complex<mplapack_binary80_t>[n*n]; std::complex<mplapack_binary80_t> *aorg = new std::complex<mplapack_binary80_t>[n*n]; std::complex<mplapack_binary80_t> *b = new std::complex<mplapack_binary80_t>[n*nrhs]; std::complex<mplapack_binary80_t> *borg = new std::complex<mplapack_binary80_t>[n*nrhs]; std::complex<mplapack_binary80_t> *xexact = new std::complex<mplapack_binary80_t>[n*nrhs];
    a[0]=std::complex<mplapack_binary80_t>(5.0,0.0); a[1]=std::complex<mplapack_binary80_t>(1.0,-1.0); a[2]=std::complex<mplapack_binary80_t>(1.0,1.0); a[3]=std::complex<mplapack_binary80_t>(4.0,0.0);
    xexact[0]=std::complex<mplapack_binary80_t>(1.0,1.0); xexact[1]=std::complex<mplapack_binary80_t>(2.0,-1.0); xexact[0+n]=std::complex<mplapack_binary80_t>(-1.0,0.0); xexact[1+n]=std::complex<mplapack_binary80_t>(0.0,2.0);
    for(mplapackint i=0;i<n*n;i++) aorg[i]=a[i];
    for(mplapackint j=0;j<nrhs;j++) for(mplapackint i=0;i<n;i++){ b[i+j*ldb]=std::complex<mplapack_binary80_t>(0.0,0.0); for(mplapackint k=0;k<n;k++) b[i+j*ldb]=b[i+j*ldb]+aorg[i+k*lda]*xexact[k+j*n]; borg[i+j*ldb]=b[i+j*ldb]; }
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    Cpotrf("L", n, a, lda, info);
    if (info == 0) Cpotrs("L", n, nrhs, a, lda, b, ldb, info);
    printf("x = "); printmat(n,nrhs,b,ldb); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n");
    delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a; return info != 0 ? 1 : 0;
}
