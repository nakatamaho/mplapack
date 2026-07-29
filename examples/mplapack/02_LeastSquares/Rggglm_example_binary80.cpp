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
mplapack_binary80_t maxabs(mplapack_binary80_t a, mplapack_binary80_t b) {
    mplapack_binary80_t d = abs(a - b);
    return d;
}

mplapack_binary80_t max_solution_error(mplapackint n, mplapackint nrhs, mplapack_binary80_t *x, mplapackint ldx, mplapack_binary80_t *xexact, mplapackint ldxexact) {
    mplapack_binary80_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mplapack_binary80_t d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mplapack_binary80_t max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mplapack_binary80_t *a, mplapackint lda, mplapack_binary80_t *x, mplapackint ldx, mplapack_binary80_t *b, mplapackint ldb) {
    mplapack_binary80_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mplapack_binary80_t s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mplapack_binary80_t d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main(){ mplapackint n=3,m=2,p=3,lda=n,ldb=n,info,lwork=-1; mplapack_binary80_t *a=new mplapack_binary80_t[lda*m]; mplapack_binary80_t *b=new mplapack_binary80_t[ldb*p]; mplapack_binary80_t *d=new mplapack_binary80_t[n]; mplapack_binary80_t *x=new mplapack_binary80_t[m]; mplapack_binary80_t *y=new mplapack_binary80_t[p]; mplapack_binary80_t *xexact=new mplapack_binary80_t[m]; mplapack_binary80_t *yexact=new mplapack_binary80_t[p]; for(mplapackint i=0;i<lda*m;i++) a[i]=0; a[0]=1; a[1]=0; a[2]=1; a[0+lda]=0; a[1+lda]=1; a[2+lda]=1; for(mplapackint i=0;i<ldb*p;i++) b[i]=0; for(mplapackint i=0;i<n;i++) b[i+i*ldb]=1; xexact[0]=1; xexact[1]=2; yexact[0]=mplapack_binary80_t(0.5); yexact[1]=mplapack_binary80_t(-0.5); yexact[2]=1; for(mplapackint i=0;i<n;i++) d[i]=a[i]*xexact[0]+a[i+lda]*xexact[1]+yexact[i]; mplapack_binary80_t wk; Rggglm(n,m,p,a,lda,b,ldb,d,x,y,&wk,lwork,info); lwork=castINTEGER_binary80(wk); mplapack_binary80_t *work=new mplapack_binary80_t[lwork]; Rggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info); printf("x = "); printvec(x,m); printf("\n"); printf("y = "); printvec(y,p); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(m,(mplapackint)1,x,m,xexact,m)); printf("\n"); printf("max |y-y_exact| = "); printnum(max_solution_error(p,(mplapackint)1,y,p,yexact,p)); printf("\n"); delete[] work; delete[] yexact; delete[] xexact; delete[] y; delete[] x; delete[] d; delete[] b; delete[] a; return info!=0?1:0; }
