REAL maxabs(REAL a, REAL b) {
    REAL d = abs(a - b);
    return d;
}

REAL max_solution_error(INTEGER n, INTEGER nrhs, REAL *x, INTEGER ldx, REAL *xexact, INTEGER ldxexact) {
    REAL err = 0.0;
    for (INTEGER j = 0; j < nrhs; j++) {
        for (INTEGER i = 0; i < n; i++) {
            REAL d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

REAL max_residual(INTEGER m, INTEGER n, INTEGER nrhs, REAL *a, INTEGER lda, REAL *x, INTEGER ldx, REAL *b, INTEGER ldb) {
    REAL err = 0.0;
    for (INTEGER j = 0; j < nrhs; j++) {
        for (INTEGER i = 0; i < m; i++) {
            REAL s = 0.0;
            for (INTEGER k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            REAL d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main(){ INTEGER n=3,m=2,p=3,lda=n,ldb=n,info,lwork=-1; REAL *a=new REAL[lda*m]; REAL *b=new REAL[ldb*p]; REAL *d=new REAL[n]; REAL *x=new REAL[m]; REAL *y=new REAL[p]; REAL *xexact=new REAL[m]; REAL *yexact=new REAL[p]; for(INTEGER i=0;i<lda*m;i++) a[i]=0; a[0]=1; a[1]=0; a[2]=1; a[0+lda]=0; a[1+lda]=1; a[2+lda]=1; for(INTEGER i=0;i<ldb*p;i++) b[i]=0; for(INTEGER i=0;i<n;i++) b[i+i*ldb]=1; xexact[0]=1; xexact[1]=2; yexact[0]=REAL(0.5); yexact[1]=REAL(-0.5); yexact[2]=1; for(INTEGER i=0;i<n;i++) d[i]=a[i]*xexact[0]+a[i+lda]*xexact[1]+yexact[i]; REAL wk; Rggglm(n,m,p,a,lda,b,ldb,d,x,y,&wk,lwork,info); lwork=castInTEGER(wk); REAL *work=new REAL[lwork]; Rggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info); printf("x = "); printvec(x,m); printf("\n"); printf("y = "); printvec(y,p); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(m,(INTEGER)1,x,m,xexact,m)); printf("\n"); printf("max |y-y_exact| = "); printnum(max_solution_error(p,(INTEGER)1,y,p,yexact,p)); printf("\n"); delete[] work; delete[] yexact; delete[] xexact; delete[] y; delete[] x; delete[] d; delete[] b; delete[] a; return info!=0?1:0; }
