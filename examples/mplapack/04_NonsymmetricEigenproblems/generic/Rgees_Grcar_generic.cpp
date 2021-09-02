bool rselect(REAL ar, REAL ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    INTEGER n = 100;

    REAL *a = new REAL[n * n];
    REAL *vs = new REAL[n * n];
    INTEGER sdim = 0;
    INTEGER lwork = 3 * n;
    REAL *wr = new REAL[n];
    REAL *wi = new REAL[n];
    REAL *work = new REAL[lwork];
    bool bwork[n];
    INTEGER info;
    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = 0.0;
            if ( i <= j && j <= i + 3 ) a[(i - 1) + (j - 1) * n] = 1.0;
            if ( i - 1 == j ) a[(i - 1) + (j - 1) * n] = -1.0;
	}
    }
    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "N", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("vs*vs'\n");
    printf("eig(a)\n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }
    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vs;
    delete[] a;
}
