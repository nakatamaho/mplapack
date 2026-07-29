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

int main() {
    INTEGER n = 2, nrhs = 2, lda = n, ldb = n, info;
    COMPLEX *a = new COMPLEX[n*n]; COMPLEX *aorg = new COMPLEX[n*n]; COMPLEX *b = new COMPLEX[n*nrhs]; COMPLEX *borg = new COMPLEX[n*nrhs]; COMPLEX *xexact = new COMPLEX[n*nrhs];
    a[0]=COMPLEX(5.0,0.0); a[1]=COMPLEX(1.0,-1.0); a[2]=COMPLEX(1.0,1.0); a[3]=COMPLEX(4.0,0.0);
    xexact[0]=COMPLEX(1.0,1.0); xexact[1]=COMPLEX(2.0,-1.0); xexact[0+n]=COMPLEX(-1.0,0.0); xexact[1+n]=COMPLEX(0.0,2.0);
    for(INTEGER i=0;i<n*n;i++) aorg[i]=a[i];
    for(INTEGER j=0;j<nrhs;j++) for(INTEGER i=0;i<n;i++){ b[i+j*ldb]=COMPLEX(0.0,0.0); for(INTEGER k=0;k<n;k++) b[i+j*ldb]=b[i+j*ldb]+aorg[i+k*lda]*xexact[k+j*n]; borg[i+j*ldb]=b[i+j*ldb]; }
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    Cpotrf("L", n, a, lda, info);
    if (info == 0) Cpotrs("L", n, nrhs, a, lda, b, ldb, info);
    printf("x = "); printmat(n,nrhs,b,ldb); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n");
    delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a; return info != 0 ? 1 : 0;
}
