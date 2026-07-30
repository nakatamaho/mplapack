int main() {
    INTEGER m = 8, n = 4, lda = m, info, lwork = -1;
    REAL *a = new REAL[m * n];
    REAL *b = new REAL[m * n];
    for (INTEGER j = 0; j < n; j++)
        for (INTEGER i = 0; i < m; i++) {
            REAL scale = pow(REAL(10.0), REAL((double)(i - 4)));
            a[i + j * lda] = scale / (REAL((double)(i + j + 1)));
            b[i + j * lda] = a[i + j * lda];
        }
    REAL *s1 = new REAL[n];
    REAL *s2 = new REAL[n];
    REAL *u = new REAL[1];
    REAL *vt = new REAL[1];
    REAL wk;
    Rgesvd("N", "N", m, n, a, lda, s1, u, (INTEGER)1, vt, (INTEGER)1, &wk, lwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rgesvd("N", "N", m, n, a, lda, s1, u, (INTEGER)1, vt, (INTEGER)1, work, lwork, info);
    delete[] work;
    INTEGER *iwork = new INTEGER[m + 3 * n + 10];
    lwork = 5 * n * n + 9 * n + m;
    if (lwork < 2 * m + n)
        lwork = 2 * m + n;
    if (lwork < 4 * n + n * n)
        lwork = 4 * n + n * n;
    if (lwork < 7)
        lwork = 7;
    work = new REAL[lwork];
    Rgejsv("G", "N", "N", "N", "N", "N", m, n, b, lda, s2, u, (INTEGER)1, vt, (INTEGER)1, work, lwork, iwork, info);
    printf("Rgesvd singular values = "); printvec(s1, n); printf("\n");
    printf("Rgejsv singular values = "); printvec(s2, n); printf("\n");
    for (INTEGER i = 0; i < n; i++) {
        REAL rel = abs(s1[i] - s2[i]) / (abs(s2[i]) + Rlamch("S"));
        printf("relative difference[%ld] = ", (long)i); printnum(rel); printf("\n");
    }
    delete[] work;
    delete[] iwork;
    delete[] vt;
    delete[] u;
    delete[] s2;
    delete[] s1;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
