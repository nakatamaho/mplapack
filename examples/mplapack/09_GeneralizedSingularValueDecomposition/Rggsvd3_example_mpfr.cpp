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

int main(){ mplapackint m=2,n=3,p=2,k,l,lda=m,ldb=p,ldu=m,ldv=p,ldq=n,info,lwork=-1; mpreal *a=new mpreal[lda*n]; mpreal *b=new mpreal[ldb*n]; mpreal *alpha=new mpreal[n]; mpreal *beta=new mpreal[n]; mpreal *u=new mpreal[ldu*m]; mpreal *v=new mpreal[ldv*p]; mpreal *q=new mpreal[ldq*n]; mplapackint *iwork=new mplapackint[n]; for(mplapackint i=0;i<lda*n;i++) a[i]=0; for(mplapackint i=0;i<ldb*n;i++) b[i]=0; a[0]=1; a[1+lda]=2; a[0+2*lda]=1; b[0]=1; b[1+ldb]=3; b[0+2*ldb]=1; mpreal wk; Rggsvd3("U","V","Q",m,n,p,k,l,a,lda,b,ldb,alpha,beta,u,ldu,v,ldv,q,ldq,&wk,lwork,iwork,info); lwork=castINTEGER_mpfr(wk); mpreal *work=new mpreal[lwork]; Rggsvd3("U","V","Q",m,n,p,k,l,a,lda,b,ldb,alpha,beta,u,ldu,v,ldv,q,ldq,work,lwork,iwork,info); printf("k = %ld, l = %ld\n",(long)k,(long)l); printf("alpha = "); printvec(alpha,n); printf("\n"); printf("beta = "); printvec(beta,n); printf("\n"); for(mplapackint i=k;i<k+l;i++){ printf("gsv[%ld] = ",(long)i); if(abs(beta[i])<=Rlamch_mpfr("E")) printf("Inf\n"); else { printnum(alpha[i]/beta[i]); printf("\n"); }} delete[] work; delete[] iwork; delete[] q; delete[] v; delete[] u; delete[] beta; delete[] alpha; delete[] b; delete[] a; return info!=0?1:0; }
