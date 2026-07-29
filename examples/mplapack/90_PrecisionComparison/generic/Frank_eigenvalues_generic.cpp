void sort_real(INTEGER n, REAL *x) {
    for (INTEGER i = 0; i < n; i++)
        for (INTEGER j = i + 1; j < n; j++)
            if (x[j] < x[i]) {
                REAL t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
}
int main() {
    INTEGER n = 20, lda = n, ldv = 1, info, lwork = -1;
    REAL *a = new REAL[n * n];
    for (INTEGER j = 0; j < n; j++)
        for (INTEGER i = 0; i < n; i++)
            a[i + j * lda] = (i <= j) ? REAL(n - j) : REAL(n - i);
    REAL *wr = new REAL[n];
    REAL *wi = new REAL[n];
    REAL *vl = new REAL[1];
    REAL *vr = new REAL[1];
    REAL wk;
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, work, lwork, info);
    sort_real(n, wr);
    printf("Frank eigenvalues = "); printvec(wr, n); printf("\n");
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] wi;
    delete[] wr;
    delete[] a;
    return info != 0 ? 1 : 0;
}
