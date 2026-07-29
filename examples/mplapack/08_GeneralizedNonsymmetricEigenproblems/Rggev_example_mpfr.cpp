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

int main(){ mplapackint n=3,lda=n,ldb=n,ldv=1,info,lwork=-1; mpreal *a=new mpreal[n*n]; mpreal *b=new mpreal[n*n]; mpreal *alphar=new mpreal[n]; mpreal *alphai=new mpreal[n]; mpreal *beta=new mpreal[n]; mpreal *vl=new mpreal[1]; mpreal *vr=new mpreal[1]; for(mplapackint i=0;i<n*n;i++){a[i]=0;b[i]=0;} a[0]=1; a[4]=2; a[8]=3; b[0]=1; b[4]=1; b[8]=0; mpreal wk; Rggev("N","N",n,a,lda,b,ldb,alphar,alphai,beta,vl,ldv,vr,ldv,&wk,lwork,info); lwork=castINTEGER_mpfr(wk); mpreal *work=new mpreal[lwork]; Rggev("N","N",n,a,lda,b,ldb,alphar,alphai,beta,vl,ldv,vr,ldv,work,lwork,info); for(mplapackint i=0;i<n;i++){ printf("alpha[%ld] = ",(long)i); printnum(alphar[i]); printf(" + "); printnum(alphai[i]); printf("i, beta = "); printnum(beta[i]); if(abs(beta[i])<=Rlamch_mpfr("E")){ printf(", lambda = Inf\n"); } else { printf(", lambda = "); printnum(alphar[i]/beta[i]); printf(" + "); printnum(alphai[i]/beta[i]); printf("i\n"); }} delete[] work; delete[] vr; delete[] vl; delete[] beta; delete[] alphai; delete[] alphar; delete[] b; delete[] a; return info!=0?1:0; }
