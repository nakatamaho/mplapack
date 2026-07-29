//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_dd.h>
#include <mplapack_dd.h>

#define DD_PRECISION_SHORT 16

inline void printnum(dd_real rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(dd_real *a, int len) {
    dd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, dd_real * a, int lda)
{
    dd_real mtmp;
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
void sort_real(mplapackint n, dd_real *x){ for(mplapackint i=0;i<n;i++) for(mplapackint j=i+1;j<n;j++) if(x[j]<x[i]){ dd_real t=x[i]; x[i]=x[j]; x[j]=t; } }
int main(){ mplapackint n=20,lda=n,ldv=1,info,lwork=-1; dd_real *coef=new dd_real[n+1]; for(mplapackint i=0;i<=n;i++) coef[i]=0; coef[0]=1; for(mplapackint k=1;k<=n;k++){ for(mplapackint j=k;j>=1;j--) coef[j]=coef[j]-dd_real(k)*coef[j-1]; } dd_real *a=new dd_real[n*n]; for(mplapackint i=0;i<n*n;i++) a[i]=0; for(mplapackint i=1;i<n;i++) a[i+(i-1)*lda]=1; for(mplapackint j=0;j<n;j++) a[j+(n-1)*lda]=-coef[n-j]/coef[0]; dd_real *wr=new dd_real[n]; dd_real *wi=new dd_real[n]; dd_real *vl=new dd_real[1]; dd_real *vr=new dd_real[1]; dd_real wk; Rgeev("N","N",n,a,lda,wr,wi,vl,ldv,vr,ldv,&wk,lwork,info); lwork=castINTEGER_dd(wk); dd_real *work=new dd_real[lwork]; Rgeev("N","N",n,a,lda,wr,wi,vl,ldv,vr,ldv,work,lwork,info); sort_real(n,wr); dd_real maxerr=0; for(mplapackint i=0;i<n;i++){ dd_real err=abs(wr[i]-dd_real(i+1)); if(maxerr<err) maxerr=err; printf("root[%ld] = ",(long)i); printnum(wr[i]); printf(", error = "); printnum(err); printf("\n"); } printf("max root error = "); printnum(maxerr); printf("\n"); delete[] work; delete[] vr; delete[] vl; delete[] wi; delete[] wr; delete[] a; delete[] coef; return info!=0?1:0; }
