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

//Matlab/Octave format
void printvec(mplapack_binary80_t *a, int len) {
    mplapack_binary80_t tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mplapack_binary80_t *a, int lda)
{
    mplapack_binary80_t mtmp;

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
mplapack_binary80_t binom(mplapackint n, mplapackint k){ mplapack_binary80_t r=1; for(mplapackint i=1;i<=k;i++) r=r*mplapack_binary80_t(n-k+i)/mplapack_binary80_t(i); return r; }
mplapack_binary80_t nearest_integer_error(mplapack_binary80_t x){ mplapack_binary80_t f=floor(x); mplapack_binary80_t c=f+1; mplapack_binary80_t df=abs(x-f); mplapack_binary80_t dc=abs(x-c); return df<dc?df:dc; }
int main(){ mplapackint n=8,lda=n,info,lwork=-1; mplapack_binary80_t *a=new mplapack_binary80_t[n*n]; mplapackint *ipiv=new mplapackint[n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<n;i++) a[i+j*lda]=binom(i+j,i); Rgetrf(n,n,a,lda,ipiv,info); mplapack_binary80_t wk; if(info==0) Rgetri(n,a,lda,ipiv,&wk,lwork,info); lwork=castINTEGER_binary80(wk); mplapack_binary80_t *work=new mplapack_binary80_t[lwork]; if(info==0) Rgetri(n,a,lda,ipiv,work,lwork,info); mplapack_binary80_t err=0; for(mplapackint i=0;i<n*n;i++){ mplapack_binary80_t d=nearest_integer_error(a[i]); if(err<d) err=d; } printf("P inverse = "); printmat(n,n,a,lda); printf("\n"); printf("max distance to integer = "); printnum(err); printf("\n"); delete[] work; delete[] ipiv; delete[] a; return info!=0?1:0; }
