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

// Matlab/Octave format
void printvec(mpreal *a, int len) {
    mpreal tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpreal *a, int lda) {
    mpreal mtmp;
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

mpreal max_eigen_residual(mplapackint n, mpreal *a, mpreal *b, mpreal lambda, mpreal *z) {
    mpreal err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        mpreal s = 0.0, t = 0.0;
        for (mplapackint j = 0; j < n; j++) { s = s + a[i + j * n] * z[j]; t = t + b[i + j * n] * z[j]; }
        mpreal d = abs(s - lambda * t);
        if (err < d) err = d;
    }
    return err;
}

int main(){ mplapackint n=2,lda=n,ldb=n,info,lwork=-1; mpreal *a=new mpreal[n*n]; mpreal *b=new mpreal[n*n]; mpreal *aorg=new mpreal[n*n]; mpreal *borg=new mpreal[n*n]; mpreal *w=new mpreal[n]; for(mplapackint i=0;i<n*n;i++){a[i]=0;b[i]=0;} a[0]=2; a[3]=6; b[0]=1; b[3]=2; for(mplapackint i=0;i<n*n;i++){aorg[i]=a[i]; borg[i]=b[i];} mpreal wk; Rsygv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,&wk,lwork,info); lwork=castINTEGER_mpfr(wk); mpreal *work=new mpreal[lwork]; Rsygv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,work,lwork,info); printf("eigenvalues = "); printvec(w,n); printf("\n"); printf("eigenvectors = "); printmat(n,n,a,lda); printf("\n"); for(mplapackint j=0;j<n;j++){ printf("residual[%ld] = ",(long)j); printnum(max_eigen_residual(n,aorg,borg,w[j],&a[j*lda])); printf("\n"); } delete[] work; delete[] w; delete[] borg; delete[] aorg; delete[] b; delete[] a; return info!=0?1:0; }
