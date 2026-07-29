//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }

// Matlab/Octave format
void printvec(double *a, int len) {
    double tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, double *a, int lda)
{
    double mtmp;

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
double maxabs(double a, double b) {
    double d = abs(a - b);
    return d;
}

double max_solution_error(mplapackint n, mplapackint nrhs, double *x, mplapackint ldx, double *xexact, mplapackint ldxexact) {
    double err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            double d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

double max_residual(mplapackint m, mplapackint n, mplapackint nrhs, double *a, mplapackint lda, double *x, mplapackint ldx, double *b, mplapackint ldb) {
    double err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            double s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            double d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main(){ mplapackint m=3,n=2,p=1,lda=m,ldb=p,info,lwork=-1; double *a=new double[lda*n]; double *bmat=new double[ldb*n]; double *c=new double[m]; double *d=new double[p]; double *x=new double[n]; double *xexact=new double[n]; a[0]=1; a[1]=0; a[2]=1; a[0+lda]=0; a[1+lda]=1; a[2+lda]=1; bmat[0]=1; bmat[0+ldb]=1; xexact[0]=1; xexact[1]=2; for(mplapackint i=0;i<m;i++) c[i]=a[i]*xexact[0]+a[i+lda]*xexact[1]; d[0]=3; double wk; Rgglse(m,n,p,a,lda,bmat,ldb,c,d,x,&wk,lwork,info); lwork=castINTEGER_double(wk); double *work=new double[lwork]; Rgglse(m,n,p,a,lda,bmat,ldb,c,d,x,work,lwork,info); printf("x = "); printvec(x,n); printf("\n"); printf("constraint B*x-d = "); printnum(x[0]+x[1]-d[0]); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(n,(mplapackint)1,x,n,xexact,n)); printf("\n"); delete[] work; delete[] xexact; delete[] x; delete[] d; delete[] c; delete[] bmat; delete[] a; return info!=0?1:0; }
