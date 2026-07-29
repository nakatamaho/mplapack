REAL max_eigen_residual(INTEGER n, REAL *a, REAL *b, REAL lambda, REAL *z) {
    REAL err = 0.0;
    for (INTEGER i = 0; i < n; i++) {
        REAL s = 0.0, t = 0.0;
        for (INTEGER j = 0; j < n; j++) {
            s = s + a[i + j * n] * z[j];
            t = t + b[i + j * n] * z[j];
        }
        REAL d = abs(s - lambda * t);
        if (err < d)
            err = d;
    }
    return err;
}

int main() {
    INTEGER n = 2, lda = n, ldb = n, info, lwork = -1;
    REAL *a = new REAL[n * n];
    REAL *b = new REAL[n * n];
    REAL *aorg = new REAL[n * n];
    REAL *borg = new REAL[n * n];
    REAL *w = new REAL[n];
    for (INTEGER i = 0; i < n * n; i++) {
        a[i] = 0;
        b[i] = 0;
    }
    a[0] = 2;
    a[3] = 6;
    b[0] = 1;
    b[3] = 2;
    for (INTEGER i = 0; i < n * n; i++) {
        aorg[i] = a[i];
        borg[i] = b[i];
    }
    REAL wk;
    Rsygv((INTEGER)1, "V", "U", n, a, lda, b, ldb, w, &wk, lwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rsygv((INTEGER)1, "V", "U", n, a, lda, b, ldb, w, work, lwork, info);
    printf("eigenvalues = "); printvec(w, n); printf("\n");
    printf("eigenvectors = "); printmat(n, n, a, lda); printf("\n");
    for (INTEGER j = 0; j < n; j++) {
        printf("residual[%ld] = ", (long)j); printnum(max_eigen_residual(n, aorg, borg, w[j], &a[j * lda])); printf("\n");
    }
    delete[] work;
    delete[] w;
    delete[] borg;
    delete[] aorg;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
