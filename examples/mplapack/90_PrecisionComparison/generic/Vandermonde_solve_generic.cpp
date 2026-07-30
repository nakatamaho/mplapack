REAL max_solution_error(INTEGER n, REAL *x, REAL *xexact) {
    REAL err = 0.0;
    for (INTEGER i = 0; i < n; i++) {
        REAL d = abs(x[i] - xexact[i]);
        if (err < d)
            err = d;
    }
    return err;
}
int main() {
    INTEGER n = 15, lda = n, ldb = n, info;
    REAL *a = new REAL[n * n];
    REAL *b = new REAL[n];
    REAL *xexact = new REAL[n];
    INTEGER *ipiv = new INTEGER[n];
    for (INTEGER j = 0; j < n; j++)
        xexact[j] = REAL(j % 3 - 1);
    for (INTEGER i = 0; i < n; i++) {
        REAL node = REAL(i + 1);
        REAL p = 1.0;
        for (INTEGER j = 0; j < n; j++) {
            a[i + j * lda] = p;
            p = p * node;
        }
    }
    for (INTEGER i = 0; i < n; i++) {
        b[i] = 0.0;
        for (INTEGER j = 0; j < n; j++)
            b[i] = b[i] + a[i + j * lda] * xexact[j];
    }
    Rgesv(n, (INTEGER)1, a, lda, ipiv, b, ldb, info);
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, b, xexact)); printf("\n");
    delete[] ipiv;
    delete[] xexact;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
