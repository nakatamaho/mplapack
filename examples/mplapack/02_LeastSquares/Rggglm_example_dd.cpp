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
dd_real maxabs(dd_real a, dd_real b) {
    dd_real d = abs(a - b);
    return d;
}

dd_real max_solution_error(mplapackint n, mplapackint nrhs, dd_real *x, mplapackint ldx, dd_real *xexact, mplapackint ldxexact) {
    dd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            dd_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

dd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, dd_real *a, mplapackint lda, dd_real *x, mplapackint ldx, dd_real *b, mplapackint ldb) {
    dd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            dd_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            dd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main(){ mplapackint n=3,m=2,p=3,lda=n,ldb=n,info,lwork=-1; dd_real *a=new dd_real[lda*m]; dd_real *b=new dd_real[ldb*p]; dd_real *d=new dd_real[n]; dd_real *x=new dd_real[m]; dd_real *y=new dd_real[p]; dd_real *xexact=new dd_real[m]; dd_real *yexact=new dd_real[p]; for(mplapackint i=0;i<lda*m;i++) a[i]=0; a[0]=1; a[1]=0; a[2]=1; a[0+lda]=0; a[1+lda]=1; a[2+lda]=1; for(mplapackint i=0;i<ldb*p;i++) b[i]=0; for(mplapackint i=0;i<n;i++) b[i+i*ldb]=1; xexact[0]=1; xexact[1]=2; yexact[0]=dd_real(0.5); yexact[1]=dd_real(-0.5); yexact[2]=1; for(mplapackint i=0;i<n;i++) d[i]=a[i]*xexact[0]+a[i+lda]*xexact[1]+yexact[i]; dd_real wk; Rggglm(n,m,p,a,lda,b,ldb,d,x,y,&wk,lwork,info); lwork=castINTEGER_dd(wk); dd_real *work=new dd_real[lwork]; Rggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info); printf("x = "); printvec(x,m); printf("\n"); printf("y = "); printvec(y,p); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(m,(mplapackint)1,x,m,xexact,m)); printf("\n"); printf("max |y-y_exact| = "); printnum(max_solution_error(p,(mplapackint)1,y,p,yexact,p)); printf("\n"); delete[] work; delete[] yexact; delete[] xexact; delete[] y; delete[] x; delete[] d; delete[] b; delete[] a; return info!=0?1:0; }
