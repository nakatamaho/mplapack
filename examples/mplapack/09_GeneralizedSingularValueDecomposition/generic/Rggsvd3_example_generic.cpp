int main() {
    INTEGER m = 2, n = 3, p = 2, k, l, lda = m, ldb = p, ldu = m, ldv = p, ldq = n, info, lwork = -1;
    REAL *a = new REAL[lda * n];
    REAL *b = new REAL[ldb * n];
    REAL *alpha = new REAL[n];
    REAL *beta = new REAL[n];
    REAL *u = new REAL[ldu * m];
    REAL *v = new REAL[ldv * p];
    REAL *q = new REAL[ldq * n];
    INTEGER *iwork = new INTEGER[n];
    for (INTEGER i = 0; i < lda * n; i++)
        a[i] = 0.0;
    for (INTEGER i = 0; i < ldb * n; i++)
        b[i] = 0.0;
    a[0] = 1.0;
    a[1 + lda] = 2.0;
    a[0 + 2 * lda] = 1.0;
    b[0] = 1.0;
    b[1 + ldb] = 3.0;
    b[0 + 2 * ldb] = 1.0;
    REAL wk;
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, &wk, lwork, iwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
    printf("k = %ld, l = %ld\n", (long)k, (long)l);
    printf("alpha = "); printvec(alpha, n); printf("\n");
    printf("beta = "); printvec(beta, n); printf("\n");
    for (INTEGER i = k; i < k + l; i++) {
        printf("gsv[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch("E"))
            printf("Inf\n");
        else {
            printnum(alpha[i] / beta[i]);
            printf("\n");
        }
    }
    delete[] work;
    delete[] iwork;
    delete[] q;
    delete[] v;
    delete[] u;
    delete[] beta;
    delete[] alpha;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
