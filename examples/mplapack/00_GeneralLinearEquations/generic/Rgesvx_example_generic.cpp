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
    INTEGER n = 5, nrhs = 1, lda = n, ldb = n, info;
    REAL *a = new REAL[n * n];
    REAL *af = new REAL[n * n];
    REAL *b = new REAL[n];
    REAL *x = new REAL[n];
    REAL *xexact = new REAL[n];
    REAL *r = new REAL[n];
    REAL *c = new REAL[n];
    REAL *ferr = new REAL[nrhs];
    REAL *berr = new REAL[nrhs];
    REAL *work = new REAL[4 * n];
    INTEGER *iwork = new INTEGER[n];
    INTEGER *ipiv = new INTEGER[n];
    char equed[2];
    equed[0] = 'N';
    equed[1] = '\0';
    REAL rcond;
    for (INTEGER j = 0; j < n; j++) for (INTEGER i = 0; i < n; i++) a[i + j * lda] = REAL(1.0) / REAL((double)(i + j + 1));
    for (INTEGER i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? REAL(1.0) : REAL(-1.0);
    for (INTEGER i = 0; i < n; i++) {
        b[i] = 0.0;
        for (INTEGER k = 0; k < n; k++)
            b[i] = b[i] + a[i + k * lda] * xexact[k];
    }
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Rgesvx("N", "N", n, nrhs, a, lda, af, lda, ipiv, equed, r, c, b, ldb, x, ldb, rcond, ferr, berr, work, iwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("rcond = "); printnum(rcond); printf("\n");
    printf("ferr = "); printvec(ferr, nrhs); printf("\n");
    printf("berr = "); printvec(berr, nrhs); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, x, ldb, xexact, n)); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] berr;
    delete[] ferr;
    delete[] c;
    delete[] r;
    delete[] xexact;
    delete[] x;
    delete[] b;
    delete[] af;
    delete[] a;
    return info != 0 ? 1 : 0;
}
