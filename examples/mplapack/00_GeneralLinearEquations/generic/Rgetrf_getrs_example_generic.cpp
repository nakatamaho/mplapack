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
    INTEGER n = 3, nrhs = 2, lda = n, ldb = n, info;
    REAL *a = new REAL[n * n];
    REAL *aorg = new REAL[n * n];
    REAL *b = new REAL[nrhs * ldb];
    REAL *borg = new REAL[nrhs * ldb];
    REAL *xexact = new REAL[nrhs * n];
    INTEGER *ipiv = new INTEGER[n];

    a[0 + 0 * lda] = 4; a[0 + 1 * lda] = 1; a[0 + 2 * lda] = 2;
    a[1 + 0 * lda] = 0; a[1 + 1 * lda] = 3; a[1 + 2 * lda] = -1;
    a[2 + 0 * lda] = 2; a[2 + 1 * lda] = 1; a[2 + 2 * lda] = 5;
    xexact[0] = 1; xexact[1] = -1; xexact[2] = 2;
    xexact[0 + 1 * n] = 2; xexact[1 + 1 * n] = 0; xexact[2 + 1 * n] = -1;
    for (INTEGER i = 0; i < n * n; i++) aorg[i] = a[i];
    for (INTEGER j = 0; j < nrhs; j++) for (INTEGER i = 0; i < n; i++) {
        b[i + j * ldb] = 0;
        for (INTEGER k = 0; k < n; k++) b[i + j * ldb] = b[i + j * ldb] + aorg[i + k * lda] * xexact[k + j * n];
        borg[i + j * ldb] = b[i + j * ldb];
    }

    printf("A = "); printmat(n, n, aorg, lda); printf("\n");
    printf("B = "); printmat(n, nrhs, b, ldb); printf("\n");
    Rgetrf(n, n, a, lda, ipiv, info);
    if (info == 0) Rgetrs("N", n, nrhs, a, lda, ipiv, b, ldb, info);
    printf("X = "); printmat(n, nrhs, b, ldb); printf("\n");
    printf("max |X-X_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");

    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
