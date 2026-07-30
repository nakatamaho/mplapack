int main() {
    INTEGER n = 12, m = n, lda = m, info, lwork = -1;
    REAL theta = REAL(0.1), c = cos(theta), s = sin(theta);
    REAL *a = new REAL[m * n];
    REAL *asvd = new REAL[m * n];
    REAL *aqr = new REAL[m * n];
    for (INTEGER j = 0; j < n; j++)
        for (INTEGER i = 0; i < m; i++) {
            REAL scale = pow(s, REAL((double)i));
            REAL val = (i == j) ? REAL(1.0) : ((i < j) ? -c : REAL(0.0));
            a[i + j * lda] = scale * val;
            asvd[i + j * lda] = a[i + j * lda];
            aqr[i + j * lda] = a[i + j * lda];
        }
    REAL *sigma = new REAL[n];
    REAL *u = new REAL[1];
    REAL *vt = new REAL[1];
    REAL wk;
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (INTEGER)1, vt, (INTEGER)1, &wk, lwork, info);
    lwork = castInTEGER(wk);
    REAL *work = new REAL[lwork];
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (INTEGER)1, vt, (INTEGER)1, work, lwork, info);
    delete[] work;
    INTEGER *jpvt = new INTEGER[n];
    REAL *tau = new REAL[n];
    for (INTEGER i = 0; i < n; i++)
        jpvt[i] = 0;
    lwork = -1;
    Rgeqp3(m, n, aqr, lda, jpvt, tau, &wk, lwork, info);
    lwork = castInTEGER(wk);
    work = new REAL[lwork];
    Rgeqp3(m, n, aqr, lda, jpvt, tau, work, lwork, info);
    printf("smallest singular value = "); printnum(sigma[n - 1]); printf("\n");
    printf("|R(n,n)| from Rgeqp3 = "); printnum(abs(aqr[n - 1 + (n - 1) * lda])); printf("\n");
    delete[] work;
    delete[] tau;
    delete[] jpvt;
    delete[] vt;
    delete[] u;
    delete[] sigma;
    delete[] aqr;
    delete[] asvd;
    delete[] a;
    return info != 0 ? 1 : 0;
}
