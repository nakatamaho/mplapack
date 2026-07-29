//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
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
mpreal max_eigen_residual(mplapackint n, mpcomplex *a, mpcomplex *b, mpreal lambda, mpcomplex *z) {
    mpreal err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        mpcomplex s = mpcomplex(0.0, 0.0), t = mpcomplex(0.0, 0.0);
        for (mplapackint j = 0; j < n; j++) { s = s + a[i + j * n] * z[j]; t = t + b[i + j * n] * z[j]; }
        mpreal d = abs(s - lambda * t);
        if (err < d) err = d;
    }
    return err;
}

int main(){ mplapackint n=2,lda=n,ldb=n,info,lwork=-1; mpcomplex *a=new mpcomplex[n*n]; mpcomplex *b=new mpcomplex[n*n]; mpcomplex *aorg=new mpcomplex[n*n]; mpcomplex *borg=new mpcomplex[n*n]; mpreal *w=new mpreal[n]; mpreal *rwork=new mpreal[3*n]; for(mplapackint i=0;i<n*n;i++){a[i]=mpcomplex(0.0,0.0);b[i]=mpcomplex(0.0,0.0);} a[0]=mpcomplex(2.0,0.0); a[3]=mpcomplex(6.0,0.0); b[0]=mpcomplex(1.0,0.0); b[3]=mpcomplex(2.0,0.0); for(mplapackint i=0;i<n*n;i++){aorg[i]=a[i]; borg[i]=b[i];} mpcomplex wk; Chegv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,&wk,lwork,rwork,info); lwork=castINTEGER_mpfr(wk.real()); mpcomplex *work=new mpcomplex[lwork]; Chegv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,work,lwork,rwork,info); printf("eigenvalues = "); printvec(w,n); printf("\n"); printf("eigenvectors = "); printmat(n,n,a,lda); printf("\n"); for(mplapackint j=0;j<n;j++){ printf("residual[%ld] = ",(long)j); printnum(max_eigen_residual(n,aorg,borg,w[j],&a[j*lda])); printf("\n"); } delete[] work; delete[] rwork; delete[] w; delete[] borg; delete[] aorg; delete[] b; delete[] a; return info!=0?1:0; }
