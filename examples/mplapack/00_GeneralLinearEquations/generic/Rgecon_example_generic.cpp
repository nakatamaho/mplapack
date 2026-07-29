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

REAL one_norm(INTEGER n, REAL *a, INTEGER lda) {
    REAL anorm = 0.0;
    for (INTEGER j = 0; j < n; j++) { REAL s = 0.0; for (INTEGER i = 0; i < n; i++) s = s + abs(a[i + j * lda]); if (anorm < s) anorm = s; }
    return anorm;
}
REAL inf_norm(INTEGER n, REAL *a, INTEGER lda) {
    REAL anorm = 0.0;
    for (INTEGER i = 0; i < n; i++) { REAL s = 0.0; for (INTEGER j = 0; j < n; j++) s = s + abs(a[i + j * lda]); if (anorm < s) anorm = s; }
    return anorm;
}
int main() {
    INTEGER n = 3, lda = n, info;
    REAL *a = new REAL[n * n];
    REAL *lu = new REAL[n * n];
    REAL *work = new REAL[4 * n];
    INTEGER *iwork = new INTEGER[n];
    INTEGER *ipiv = new INTEGER[n];
    for (INTEGER i = 0; i < n * n; i++) a[i] = 0;
    a[0 + 0 * lda] = 1; a[1 + 1 * lda] = REAL(1.0e-3); a[2 + 2 * lda] = REAL(1.0e-6);
    for (INTEGER i = 0; i < n * n; i++) lu[i] = a[i];
    Rgetrf(n, n, lu, lda, ipiv, info);
    REAL rcond1 = 0.0, rcondi = 0.0;
    if (info == 0) Rgecon("1", n, lu, lda, one_norm(n, a, lda), rcond1, work, iwork, info);
    if (info == 0) Rgecon("I", n, lu, lda, inf_norm(n, a, lda), rcondi, work, iwork, info);
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("true cond_1 = "); printnum(REAL(1.0e6)); printf("\n");
    printf("rcond_1 = "); printnum(rcond1); printf("\n");
    printf("rcond_inf = "); printnum(rcondi); printf("\n");
    delete[] ipiv; delete[] iwork; delete[] work; delete[] lu; delete[] a;
    return info != 0 ? 1 : 0;
}
