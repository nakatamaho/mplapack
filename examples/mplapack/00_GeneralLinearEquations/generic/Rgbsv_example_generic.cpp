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
    INTEGER n = 6, kl = 2, ku = 1, nrhs = 1, ldab = 2 * kl + ku + 1, ldb = n, info;
    REAL *a = new REAL[n * n];
    REAL *ab = new REAL[ldab * n];
    REAL *b = new REAL[n];
    REAL *borg = new REAL[n];
    REAL *xexact = new REAL[n];
    INTEGER *ipiv = new INTEGER[n];
    for (INTEGER i = 0; i < n * n; i++) a[i] = 0.0;
    for (INTEGER i = 0; i < ldab * n; i++) ab[i] = 0.0;
    for (INTEGER i = 0; i < n; i++) {
        xexact[i] = i + 1.0;
        for (INTEGER j = 0; j < n; j++) {
            if (i == j) a[i + j * n] = 6.0;
            else if (i > j && i - j <= kl) a[i + j * n] = -1.0;
            else if (j > i && j - i <= ku) a[i + j * n] = 2.0;
        }
    }
    for (INTEGER j = 0; j < n; j++) for (INTEGER i = 0; i < n; i++) if (a[i + j * n] != 0.0) {
        INTEGER row = kl + ku + i - j;
        ab[row + j * ldab] = a[i + j * n];
    }
    for (INTEGER i = 0; i < n; i++) {
        b[i] = 0.0;
        for (INTEGER k = 0; k < n; k++) b[i] = b[i] + a[i + k * n] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n, n, a, n); printf("\n");
    printf("AB = "); printmat(ldab, n, ab, ldab); printf("\n");
    Rgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, a, n, b, ldb, borg, ldb)); printf("\n");
    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] ab;
    delete[] a;
    return info != 0 ? 1 : 0;
}
