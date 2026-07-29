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
qd_real maxabs(qd_real a, qd_real b) {
    qd_real d = abs(a - b);
    return d;
}

qd_real max_solution_error(mplapackint n, mplapackint nrhs, qd_real *x, mplapackint ldx, qd_real *xexact, mplapackint ldxexact) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            qd_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

qd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, qd_real *a, mplapackint lda, qd_real *x, mplapackint ldx, qd_real *b, mplapackint ldb) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            qd_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            qd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main(){ mplapackint m=3,n=2,p=1,lda=m,ldb=p,info,lwork=-1; qd_real *a=new qd_real[lda*n]; qd_real *bmat=new qd_real[ldb*n]; qd_real *c=new qd_real[m]; qd_real *d=new qd_real[p]; qd_real *x=new qd_real[n]; qd_real *xexact=new qd_real[n]; a[0]=1; a[1]=0; a[2]=1; a[0+lda]=0; a[1+lda]=1; a[2+lda]=1; bmat[0]=1; bmat[0+ldb]=1; xexact[0]=1; xexact[1]=2; for(mplapackint i=0;i<m;i++) c[i]=a[i]*xexact[0]+a[i+lda]*xexact[1]; d[0]=3; qd_real wk; Rgglse(m,n,p,a,lda,bmat,ldb,c,d,x,&wk,lwork,info); lwork=castINTEGER_qd(wk); qd_real *work=new qd_real[lwork]; Rgglse(m,n,p,a,lda,bmat,ldb,c,d,x,work,lwork,info); printf("x = "); printvec(x,n); printf("\n"); printf("constraint B*x-d = "); printnum(x[0]+x[1]-d[0]); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(n,(mplapackint)1,x,n,xexact,n)); printf("\n"); delete[] work; delete[] xexact; delete[] x; delete[] d; delete[] c; delete[] bmat; delete[] a; return info!=0?1:0; }
