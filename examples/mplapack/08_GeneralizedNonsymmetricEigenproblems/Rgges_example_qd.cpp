//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(qd_real *a, int len) {
    qd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, qd_real * a, int lda)
{
    qd_real mtmp;
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
bool select_none(qd_real ar, qd_real ai, qd_real beta){ return false; }
int main(){ mplapackint n=3,lda=n,ldb=n,ldv=1,sdim,info,lwork=-1; qd_real *a=new qd_real[n*n]; qd_real *b=new qd_real[n*n]; qd_real *alphar=new qd_real[n]; qd_real *alphai=new qd_real[n]; qd_real *beta=new qd_real[n]; qd_real *vsl=new qd_real[1]; qd_real *vsr=new qd_real[1]; bool *bwork=new bool[n]; for(mplapackint i=0;i<n*n;i++){a[i]=0;b[i]=0;} a[0]=1; a[4]=2; a[8]=3; b[0]=1; b[4]=1; b[8]=0; qd_real wk; Rgges("N","N","N",select_none,n,a,lda,b,ldb,sdim,alphar,alphai,beta,vsl,ldv,vsr,ldv,&wk,lwork,bwork,info); lwork=castINTEGER_qd(wk); qd_real *work=new qd_real[lwork]; Rgges("N","N","N",select_none,n,a,lda,b,ldb,sdim,alphar,alphai,beta,vsl,ldv,vsr,ldv,work,lwork,bwork,info); printf("S = "); printmat(n,n,a,lda); printf("\n"); printf("T = "); printmat(n,n,b,ldb); printf("\n"); for(mplapackint i=0;i<n;i++){ printf("lambda[%ld] = ",(long)i); if(abs(beta[i])<=Rlamch_qd("E")) printf("Inf\n"); else { printnum(alphar[i]/beta[i]); printf(" + "); printnum(alphai[i]/beta[i]); printf("i\n"); }} delete[] work; delete[] bwork; delete[] vsr; delete[] vsl; delete[] beta; delete[] alphai; delete[] alphar; delete[] b; delete[] a; return info!=0?1:0; }
