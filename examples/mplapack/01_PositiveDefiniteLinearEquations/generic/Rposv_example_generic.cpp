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

int main() {
    INTEGER n = 2, nrhs = 1, lda = n, ldb = n, info;
    REAL *a = new REAL[n * n];
    REAL *aorg = new REAL[n * n];
    REAL *b = new REAL[n];
    REAL *borg = new REAL[n];
    REAL *xexact = new REAL[n];
    a[0]=4; a[1]=2; a[2]=2; a[3]=5;
    xexact[0]=1; xexact[1]=2;
    for (INTEGER i=0;i<n*n;i++) aorg[i]=a[i];
    for (INTEGER i = 0; i < n; i++) {
        b[i] = 0;
        for (INTEGER k = 0; k < n; k++)
            b[i] += aorg[i + k * lda] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    Rposv("L", n, nrhs, a, lda, b, ldb, info);
    printf("x = "); printvec(b,n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n");
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
