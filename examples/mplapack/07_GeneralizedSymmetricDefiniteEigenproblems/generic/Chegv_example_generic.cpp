REAL max_eigen_residual(INTEGER n, COMPLEX *a, COMPLEX *b, REAL lambda, COMPLEX *z) {
    REAL err = 0.0;
    for (INTEGER i = 0; i < n; i++) {
        COMPLEX s = COMPLEX(0.0, 0.0), t = COMPLEX(0.0, 0.0);
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
    COMPLEX *a = new COMPLEX[n * n];
    COMPLEX *b = new COMPLEX[n * n];
    COMPLEX *aorg = new COMPLEX[n * n];
    COMPLEX *borg = new COMPLEX[n * n];
    REAL *w = new REAL[n];
    REAL *rwork = new REAL[3 * n];
    for (INTEGER i = 0; i < n * n; i++) {
        a[i] = COMPLEX(0.0, 0.0);
        b[i] = COMPLEX(0.0, 0.0);
    }
    a[0] = COMPLEX(2.0, 0.0);
    a[3] = COMPLEX(6.0, 0.0);
    b[0] = COMPLEX(1.0, 0.0);
    b[3] = COMPLEX(2.0, 0.0);
    for (INTEGER i = 0; i < n * n; i++) {
        aorg[i] = a[i];
        borg[i] = b[i];
    }
    COMPLEX wk;
    Chegv((INTEGER)1, "V", "U", n, a, lda, b, ldb, w, &wk, lwork, rwork, info);
    lwork = castInTEGER(wk.real());
    COMPLEX *work = new COMPLEX[lwork];
    Chegv((INTEGER)1, "V", "U", n, a, lda, b, ldb, w, work, lwork, rwork, info);
    printf("eigenvalues = "); printvec(w, n); printf("\n");
    printf("eigenvectors = "); printmat(n, n, a, lda); printf("\n");
    for (INTEGER j = 0; j < n; j++) {
        printf("residual[%ld] = ", (long)j); printnum(max_eigen_residual(n, aorg, borg, w[j], &a[j * lda])); printf("\n");
    }
    delete[] work;
    delete[] rwork;
    delete[] w;
    delete[] borg;
    delete[] aorg;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
