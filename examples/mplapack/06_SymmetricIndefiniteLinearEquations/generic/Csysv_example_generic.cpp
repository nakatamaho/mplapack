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

int main(){ INTEGER n=2,nrhs=1,lda=n,ldb=n,info,lwork=-1; COMPLEX *a=new COMPLEX[n*n]; COMPLEX *aorg=new COMPLEX[n*n]; COMPLEX *b=new COMPLEX[n]; COMPLEX *borg=new COMPLEX[n]; COMPLEX *xexact=new COMPLEX[n]; INTEGER *ipiv=new INTEGER[n]; a[0]=COMPLEX(2.0,0.0); a[1]=COMPLEX(1.0,1.0); a[2]=COMPLEX(1.0,1.0); a[3]=COMPLEX(-3.0,0.0); xexact[0]=COMPLEX(1.0,-1.0); xexact[1]=COMPLEX(2.0,1.0); for(INTEGER i=0;i<n*n;i++) aorg[i]=a[i]; for(INTEGER i=0;i<n;i++){ b[i]=COMPLEX(0.0,0.0); for(INTEGER k=0;k<n;k++) b[i]=b[i]+aorg[i+k*lda]*xexact[k]; borg[i]=b[i]; } COMPLEX wk; Csysv("U",n,nrhs,a,lda,ipiv,b,ldb,&wk,lwork,info); lwork=castInTEGER(wk.real()); COMPLEX *work=new COMPLEX[lwork]; Csysv("U",n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info); printf("A = "); printmat(n,n,aorg,lda); printf("\n"); printf("x = "); printvec(b,n); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n"); printf("max residual = "); printnum(max_residual(n,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n"); delete[] work; delete[] ipiv; delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a; return info!=0?1:0; }
