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

void set_problem(INTEGER m, INTEGER n, REAL *a, INTEGER lda, REAL *b, INTEGER ldb, REAL *xexact) {
    for (INTEGER i = 0; i < m; i++) {
        a[i + 0 * lda] = 1.0;
        a[i + 1 * lda] = i;
    }
    xexact[0] = 1.0;
    xexact[1] = 2.0;
    for (INTEGER i = 0; i < m; i++)
        b[i] = a[i + 0 * lda] * xexact[0] + a[i + 1 * lda] * xexact[1];
}

int main() {
    INTEGER m = 4, n = 2, nrhs = 1, lda = m, ldb = m, info, lwork = -1, rank;
    REAL *a = new REAL[lda * n];
    REAL *aorg = new REAL[lda * n];
    REAL *b = new REAL[ldb];
    REAL *borg = new REAL[ldb];
    REAL *xexact = new REAL[n];
    INTEGER *jpvt = new INTEGER[n];
    for (INTEGER i = 0; i < n; i++)
        jpvt[i] = 0;
    set_problem(m, n, a, lda, b, ldb, xexact);
    for (INTEGER i = 0; i < lda * n; i++)
        aorg[i] = a[i];
    for (INTEGER i = 0; i < ldb; i++)
        borg[i] = b[i];
    REAL wk;
    Rgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, REAL(-1.0), rank, &wk, lwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, REAL(-1.0), rank, work, lwork, info);
    printf("A = "); printmat(m, n, aorg, lda); printf("\n");
    printf("x = "); printvec(b, n); printf("\n");
    printf("rank = %ld\n", (long)rank);
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(m, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");
    delete[] work;
    delete[] jpvt;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
