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
mplapack_binary128_t max_eigen_residual(mplapackint n, std::complex<mplapack_binary128_t> *a, std::complex<mplapack_binary128_t> *b, mplapack_binary128_t lambda, std::complex<mplapack_binary128_t> *z) {
    mplapack_binary128_t err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        std::complex<mplapack_binary128_t> s = std::complex<mplapack_binary128_t>(0.0, 0.0), t = std::complex<mplapack_binary128_t>(0.0, 0.0);
        for (mplapackint j = 0; j < n; j++) { s = s + a[i + j * n] * z[j]; t = t + b[i + j * n] * z[j]; }
        mplapack_binary128_t d = abs(s - lambda * t);
        if (err < d) err = d;
    }
    return err;
}

int main(){ mplapackint n=2,lda=n,ldb=n,info,lwork=-1; std::complex<mplapack_binary128_t> *a=new std::complex<mplapack_binary128_t>[n*n]; std::complex<mplapack_binary128_t> *b=new std::complex<mplapack_binary128_t>[n*n]; std::complex<mplapack_binary128_t> *aorg=new std::complex<mplapack_binary128_t>[n*n]; std::complex<mplapack_binary128_t> *borg=new std::complex<mplapack_binary128_t>[n*n]; mplapack_binary128_t *w=new mplapack_binary128_t[n]; mplapack_binary128_t *rwork=new mplapack_binary128_t[3*n]; for(mplapackint i=0;i<n*n;i++){a[i]=std::complex<mplapack_binary128_t>(0.0,0.0);b[i]=std::complex<mplapack_binary128_t>(0.0,0.0);} a[0]=std::complex<mplapack_binary128_t>(2.0,0.0); a[3]=std::complex<mplapack_binary128_t>(6.0,0.0); b[0]=std::complex<mplapack_binary128_t>(1.0,0.0); b[3]=std::complex<mplapack_binary128_t>(2.0,0.0); for(mplapackint i=0;i<n*n;i++){aorg[i]=a[i]; borg[i]=b[i];} std::complex<mplapack_binary128_t> wk; Chegv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,&wk,lwork,rwork,info); lwork=castINTEGER_binary128(wk.real()); std::complex<mplapack_binary128_t> *work=new std::complex<mplapack_binary128_t>[lwork]; Chegv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,work,lwork,rwork,info); printf("eigenvalues = "); printvec(w,n); printf("\n"); printf("eigenvectors = "); printmat(n,n,a,lda); printf("\n"); for(mplapackint j=0;j<n;j++){ printf("residual[%ld] = ",(long)j); printnum(max_eigen_residual(n,aorg,borg,w[j],&a[j*lda])); printf("\n"); } delete[] work; delete[] rwork; delete[] w; delete[] borg; delete[] aorg; delete[] b; delete[] a; return info!=0?1:0; }
