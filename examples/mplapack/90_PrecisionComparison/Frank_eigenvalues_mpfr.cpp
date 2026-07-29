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

void sort_real(mplapackint n, mpreal *x){ for(mplapackint i=0;i<n;i++) for(mplapackint j=i+1;j<n;j++) if(x[j]<x[i]){ mpreal t=x[i]; x[i]=x[j]; x[j]=t; } }
int main(){ mplapackint n=20,lda=n,ldv=1,info,lwork=-1; mpreal *a=new mpreal[n*n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<n;i++) a[i+j*lda]=(i<=j)?mpreal(n-j):mpreal(n-i); mpreal *wr=new mpreal[n]; mpreal *wi=new mpreal[n]; mpreal *vl=new mpreal[1]; mpreal *vr=new mpreal[1]; mpreal wk; Rgeev("N","N",n,a,lda,wr,wi,vl,ldv,vr,ldv,&wk,lwork,info); lwork=castINTEGER_mpfr(wk); mpreal *work=new mpreal[lwork]; Rgeev("N","N",n,a,lda,wr,wi,vl,ldv,vr,ldv,work,lwork,info); sort_real(n,wr); printf("Frank eigenvalues = "); printvec(wr,n); printf("\n"); delete[] work; delete[] vr; delete[] vl; delete[] wi; delete[] wr; delete[] a; return info!=0?1:0; }
