int main() {
    INTEGER n = 3, lda = n, ldb = n, ldv = 1, info, lwork = -1;
    REAL *a = new REAL[n * n];
    REAL *b = new REAL[n * n];
    REAL *alphar = new REAL[n];
    REAL *alphai = new REAL[n];
    REAL *beta = new REAL[n];
    REAL *vl = new REAL[1];
    REAL *vr = new REAL[1];
    for (INTEGER i = 0; i < n * n; i++) {
        a[i] = 0.0;
        b[i] = 0.0;
    }
    a[0] = 1.0;
    a[4] = 2.0;
    a[8] = 3.0;
    b[0] = 1.0;
    b[4] = 1.0;
    b[8] = 0.0;
    REAL wk;
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, work, lwork, info);
    for (INTEGER i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alphar[i]); printf(" + "); printnum(alphai[i]); printf("i, beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch("E")) {
            printf(", lambda = Inf\n");
        } else {
            printf(", lambda = "); printnum(alphar[i] / beta[i]); printf(" + "); printnum(alphai[i] / beta[i]); printf("i\n");
        }
    }
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] beta;
    delete[] alphai;
    delete[] alphar;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
