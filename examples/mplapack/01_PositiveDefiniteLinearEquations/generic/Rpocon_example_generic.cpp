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

REAL one_norm(INTEGER n, REAL *a, INTEGER lda) { REAL anorm=0.0; for(INTEGER j=0;j<n;j++){ REAL s=0.0; for(INTEGER i=0;i<n;i++) s=s+abs(a[i+j*lda]); if(anorm<s) anorm=s; } return anorm; }
int main() {
    INTEGER n = 2, lda = n, info;
    REAL *a = new REAL[n*n]; REAL *aorg = new REAL[n*n]; REAL *work = new REAL[3*n]; INTEGER *iwork = new INTEGER[n]; REAL rcond=0.0;
    a[0]=4; a[1]=2; a[2]=2; a[3]=5; for(INTEGER i=0;i<n*n;i++) aorg[i]=a[i];
    Rpotrf("L", n, a, lda, info);
    if (info == 0) Rpocon("L", n, a, lda, one_norm(n,aorg,lda), rcond, work, iwork, info);
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    printf("rcond_1 = "); printnum(rcond); printf("\n");
    delete[] iwork; delete[] work; delete[] aorg; delete[] a; return info != 0 ? 1 : 0;
}
