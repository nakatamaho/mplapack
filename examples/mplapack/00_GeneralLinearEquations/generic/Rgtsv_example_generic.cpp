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
    INTEGER n = 6, nrhs = 1, ldb = n, info;
    REAL *dl = new REAL[n - 1];
    REAL *d = new REAL[n];
    REAL *du = new REAL[n - 1];
    REAL *a = new REAL[n * n];
    REAL *b = new REAL[n];
    REAL *borg = new REAL[n];
    REAL *xexact = new REAL[n];
    for (INTEGER i = 0; i < n * n; i++) a[i] = 0;
    for (INTEGER i = 0; i < n; i++) { d[i] = 2; xexact[i] = i + 1; a[i + i * n] = 2; }
    for (INTEGER i = 0; i < n - 1; i++) { dl[i] = -1; du[i] = -1; a[i + 1 + i * n] = -1; a[i + (i + 1) * n] = -1; }
    for (INTEGER i = 0; i < n; i++) {
        b[i] = 0;
        for (INTEGER k = 0; k < n; k++) b[i] = b[i] + a[i + k * n] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n, n, a, n); printf("\n");
    Rgtsv(n, nrhs, dl, d, du, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, a, n, b, ldb, borg, ldb)); printf("\n");
    delete[] xexact; delete[] borg; delete[] b; delete[] a; delete[] du; delete[] d; delete[] dl;
    return info != 0 ? 1 : 0;
}
