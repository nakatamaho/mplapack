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
int main(){ mplapackint m=8,n=4,lda=m,info,lwork=-1; dd_real *a=new dd_real[m*n]; dd_real *b=new dd_real[m*n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<m;i++){ dd_real scale=pow(dd_real(10.0),dd_real(i-4)); a[i+j*lda]=scale/(dd_real(i+j+1)); b[i+j*lda]=a[i+j*lda]; } dd_real *s1=new dd_real[n]; dd_real *s2=new dd_real[n]; dd_real *u=new dd_real[1]; dd_real *vt=new dd_real[1]; dd_real wk; Rgesvd("N","N",m,n,a,lda,s1,u,(mplapackint)1,vt,(mplapackint)1,&wk,lwork,info); lwork=castINTEGER_dd(wk); dd_real *work=new dd_real[lwork]; Rgesvd("N","N",m,n,a,lda,s1,u,(mplapackint)1,vt,(mplapackint)1,work,lwork,info); delete[] work; mplapackint *iwork=new mplapackint[m+n]; lwork=-1; Rgejsv("G","N","N","N","N","N",m,n,b,lda,s2,u,(mplapackint)1,vt,(mplapackint)1,&wk,lwork,iwork,info); lwork=castINTEGER_dd(wk); work=new dd_real[lwork]; Rgejsv("G","N","N","N","N","N",m,n,b,lda,s2,u,(mplapackint)1,vt,(mplapackint)1,work,lwork,iwork,info); printf("Rgesvd singular values = "); printvec(s1,n); printf("\n"); printf("Rgejsv singular values = "); printvec(s2,n); printf("\n"); for(mplapackint i=0;i<n;i++){ dd_real rel=abs(s1[i]-s2[i])/(abs(s2[i])+Rlamch_dd("S")); printf("relative difference[%ld] = ",(long)i); printnum(rel); printf("\n"); } delete[] work; delete[] iwork; delete[] vt; delete[] u; delete[] s2; delete[] s1; delete[] b; delete[] a; return info!=0?1:0; }
