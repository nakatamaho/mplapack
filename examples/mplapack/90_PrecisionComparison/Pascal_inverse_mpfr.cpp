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

mpreal binom(mplapackint n, mplapackint k){ mpreal r=1; for(mplapackint i=1;i<=k;i++) r=r*mpreal(n-k+i)/mpreal(i); return r; }
mpreal nearest_integer_error(mpreal x){ mpreal f=floor(x); mpreal c=f+1; mpreal df=abs(x-f); mpreal dc=abs(x-c); return df<dc?df:dc; }
int main(){ mplapackint n=8,lda=n,info,lwork=-1; mpreal *a=new mpreal[n*n]; mplapackint *ipiv=new mplapackint[n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<n;i++) a[i+j*lda]=binom(i+j,i); Rgetrf(n,n,a,lda,ipiv,info); mpreal wk; if(info==0) Rgetri(n,a,lda,ipiv,&wk,lwork,info); lwork=castINTEGER_mpfr(wk); mpreal *work=new mpreal[lwork]; if(info==0) Rgetri(n,a,lda,ipiv,work,lwork,info); mpreal err=0; for(mplapackint i=0;i<n*n;i++){ mpreal d=nearest_integer_error(a[i]); if(err<d) err=d; } printf("P inverse = "); printmat(n,n,a,lda); printf("\n"); printf("max distance to integer = "); printnum(err); printf("\n"); delete[] work; delete[] ipiv; delete[] a; return info!=0?1:0; }
