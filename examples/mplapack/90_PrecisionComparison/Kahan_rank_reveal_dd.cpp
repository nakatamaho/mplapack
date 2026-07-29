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
int main(){ mplapackint n=12,m=n,lda=m,info,lwork=-1; dd_real theta=dd_real(0.1), c=cos(theta), s=sin(theta); dd_real *a=new dd_real[m*n]; dd_real *asvd=new dd_real[m*n]; dd_real *aqr=new dd_real[m*n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<m;i++){ dd_real scale=pow(s,dd_real(i)); dd_real val=(i==j)?dd_real(1.0):((i<j)?-c:dd_real(0.0)); a[i+j*lda]=scale*val; asvd[i+j*lda]=a[i+j*lda]; aqr[i+j*lda]=a[i+j*lda]; } dd_real *sigma=new dd_real[n]; dd_real *u=new dd_real[1]; dd_real *vt=new dd_real[1]; dd_real wk; Rgesvd("N","N",m,n,asvd,lda,sigma,u,(mplapackint)1,vt,(mplapackint)1,&wk,lwork,info); lwork=castINTEGER_dd(wk); dd_real *work=new dd_real[lwork]; Rgesvd("N","N",m,n,asvd,lda,sigma,u,(mplapackint)1,vt,(mplapackint)1,work,lwork,info); delete[] work; mplapackint *jpvt=new mplapackint[n]; dd_real *tau=new dd_real[n]; for(mplapackint i=0;i<n;i++) jpvt[i]=0; lwork=-1; Rgeqp3(m,n,aqr,lda,jpvt,tau,&wk,lwork,info); lwork=castINTEGER_dd(wk); work=new dd_real[lwork]; Rgeqp3(m,n,aqr,lda,jpvt,tau,work,lwork,info); printf("smallest singular value = "); printnum(sigma[n-1]); printf("\n"); printf("|R(n,n)| from Rgeqp3 = "); printnum(abs(aqr[n-1+(n-1)*lda])); printf("\n"); delete[] work; delete[] tau; delete[] jpvt; delete[] vt; delete[] u; delete[] sigma; delete[] aqr; delete[] asvd; delete[] a; return info!=0?1:0; }
