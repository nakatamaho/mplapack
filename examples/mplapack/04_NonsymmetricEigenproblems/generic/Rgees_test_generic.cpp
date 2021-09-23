bool rselect(REAL ar, REAL ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    INTEGER n = 4;

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
    a[0 + 0 * n] = -2.0;     a[0 + 1 * n] = 2.0;   a[0 + 2 * n] = 2.0;    a[0 + 3 * n] = 2.0;
    a[1 + 0 * n] = -3.0;     a[1 + 1 * n] = 3.0;   a[1 + 2 * n] = 2.0;    a[1 + 3 * n] = 2.0;
    a[2 + 0 * n] = -2.0;     a[2 + 1 * n] = 0.0;   a[2 + 2 * n] = 4.0;    a[2 + 3 * n] = 2.0;
    a[3 + 0 * n] = -1.0;     a[3 + 1 * n] = 0.0;   a[3 + 2 * n] = 0.0;    a[3 + 3 * n] = 5.0;

    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "S", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
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
