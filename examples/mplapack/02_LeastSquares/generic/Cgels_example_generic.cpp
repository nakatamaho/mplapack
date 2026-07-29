REAL max_solution_error(INTEGER n, INTEGER nrhs, COMPLEX *x, INTEGER ldx, COMPLEX *xexact, INTEGER ldxexact) {
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

REAL max_residual(INTEGER m, INTEGER n, INTEGER nrhs, COMPLEX *a, INTEGER lda, COMPLEX *x, INTEGER ldx, COMPLEX *b, INTEGER ldb) {
    REAL err = 0.0;
    for (INTEGER j = 0; j < nrhs; j++) {
        for (INTEGER i = 0; i < m; i++) {
            COMPLEX s = COMPLEX(0.0, 0.0);
            for (INTEGER k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            REAL d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

void set_problem(INTEGER m, INTEGER n, COMPLEX *a, INTEGER lda, COMPLEX *b, INTEGER ldb, COMPLEX *xexact) {
    for (INTEGER i = 0; i < m; i++) { a[i + 0 * lda] = COMPLEX(1.0, 0.0); a[i + 1 * lda] = COMPLEX(i, 1.0); }
    xexact[0] = COMPLEX(1.0, -1.0); xexact[1] = COMPLEX(2.0, 1.0);
    for (INTEGER i = 0; i < m; i++) b[i] = a[i + 0 * lda] * xexact[0] + a[i + 1 * lda] * xexact[1];
}

int main(){ INTEGER m=4,n=2,nrhs=1,lda=m,ldb=m,info,lwork=-1; COMPLEX *a=new COMPLEX[lda*n]; COMPLEX *aorg=new COMPLEX[lda*n]; COMPLEX *b=new COMPLEX[ldb]; COMPLEX *borg=new COMPLEX[ldb]; COMPLEX *xexact=new COMPLEX[n]; set_problem(m,n,a,lda,b,ldb,xexact); for(INTEGER i=0;i<lda*n;i++) aorg[i]=a[i]; for(INTEGER i=0;i<ldb;i++) borg[i]=b[i]; COMPLEX wk; Cgels("N",m,n,nrhs,a,lda,b,ldb,&wk,lwork,info); lwork=castInTEGER(wk.real()); COMPLEX *work=new COMPLEX[lwork]; Cgels("N",m,n,nrhs,a,lda,b,ldb,work,lwork,info); printf("A = "); printmat(m,n,aorg,lda); printf("\n"); printf("x = "); printvec(b,n); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n"); printf("max residual = "); printnum(max_residual(m,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n"); delete[] work; delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a; return info!=0?1:0; }
