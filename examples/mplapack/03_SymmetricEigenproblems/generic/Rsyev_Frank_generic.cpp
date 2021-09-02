void Frank(INTEGER n) {
    INTEGER lwork, info;
    REAL *a = new REAL[n * n];
    REAL *w = new REAL[n];
    REAL *lambda = new REAL[n];
    REAL *reldiff = new REAL[n];
    REAL PI;
    PI = pi(PI);

    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = n - std::max(i, j) + 1;
        }
    }
    printf("a =");
    printmat(n, n, a, n);
    printf("\n");

    // work space query
    lwork = -1;
    REAL *work = new REAL[1];

    Rsyev("V", "U", n, a, n, w, work, lwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new REAL[std::max((INTEGER)1, lwork)];

    // diagonalize matrix
    Rsyev("N", "U", n, a, n, w, work, lwork, info);

    // print out
    printf("#eigenvalues \n");
    printf("w ="); printvec(w, n); printf("\n");

    // print out
    printf("# analytic eigenvalues\n");
    for (int i = 1; i <= n; i++) {
        lambda[(n - i)] = 0.5 * 1.0 / (1.0 - cos((2.0 * i - 1.0) * PI / castReAL(2 * n + 1)));
    }
    printf("lambda ="); printvec(lambda, n); printf("\n");

    for (int i = 1; i <= n; i++) {
        reldiff[i - 1] = abs((lambda[i - 1] - w[i - 1]) / lambda[i - 1]);
    }
    printf("reldiff ="); printvec(reldiff, n); printf("\n");

    REAL maxreldiff = 0.0;
    maxreldiff = reldiff[0]; 
    for (int i = 2; i <= n; i++) {
        maxreldiff = std::max(reldiff[i - 1], maxreldiff);
    }
    printf("maxreldiff_%d =", (int)n); printnum(maxreldiff); printf("\n");

    delete[] reldiff;
    delete[] lambda;
    delete[] work;
    delete[] w;
    delete[] a;
}

int main() {
    printf("split_long_rows(0)\n");
    for (int n = 1; n < 100; n++) {
        printf("# Eigenvalues of Frank matrix of order n=%d\n", n);
        Frank((INTEGER)n);
    }
    for (int n = 100; n < 10000; n = n + 100) {
        printf("# Eigenvalues of Frank matrix of order n=%d\n", n);
        Frank((INTEGER)n);
    }
}