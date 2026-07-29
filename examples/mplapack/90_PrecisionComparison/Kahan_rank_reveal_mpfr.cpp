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

int main(){ mplapackint n=12,m=n,lda=m,info,lwork=-1; mpreal theta=mpreal(0.1), c=cos(theta), s=sin(theta); mpreal *a=new mpreal[m*n]; mpreal *asvd=new mpreal[m*n]; mpreal *aqr=new mpreal[m*n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<m;i++){ mpreal scale=pow(s,mpreal(i)); mpreal val=(i==j)?mpreal(1.0):((i<j)?-c:mpreal(0.0)); a[i+j*lda]=scale*val; asvd[i+j*lda]=a[i+j*lda]; aqr[i+j*lda]=a[i+j*lda]; } mpreal *sigma=new mpreal[n]; mpreal *u=new mpreal[1]; mpreal *vt=new mpreal[1]; mpreal wk; Rgesvd("N","N",m,n,asvd,lda,sigma,u,(mplapackint)1,vt,(mplapackint)1,&wk,lwork,info); lwork=castINTEGER_mpfr(wk); mpreal *work=new mpreal[lwork]; Rgesvd("N","N",m,n,asvd,lda,sigma,u,(mplapackint)1,vt,(mplapackint)1,work,lwork,info); delete[] work; mplapackint *jpvt=new mplapackint[n]; mpreal *tau=new mpreal[n]; for(mplapackint i=0;i<n;i++) jpvt[i]=0; lwork=-1; Rgeqp3(m,n,aqr,lda,jpvt,tau,&wk,lwork,info); lwork=castINTEGER_mpfr(wk); work=new mpreal[lwork]; Rgeqp3(m,n,aqr,lda,jpvt,tau,work,lwork,info); printf("smallest singular value = "); printnum(sigma[n-1]); printf("\n"); printf("|R(n,n)| from Rgeqp3 = "); printnum(abs(aqr[n-1+(n-1)*lda])); printf("\n"); delete[] work; delete[] tau; delete[] jpvt; delete[] vt; delete[] u; delete[] sigma; delete[] aqr; delete[] asvd; delete[] a; return info!=0?1:0; }
