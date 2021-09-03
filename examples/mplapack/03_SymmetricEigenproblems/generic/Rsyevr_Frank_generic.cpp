void Frank(INTEGER n) {
    INTEGER lwork, liwork, info, m;
    REAL *a = new REAL[n * n];
    REAL *z = new REAL[n * n]; //not used
    INTEGER *isuppz = new INTEGER[2 * n]; //not used
    REAL *w = new REAL[n];
    REAL *lambda = new REAL[n];
    REAL *reldiff = new REAL[n];
    REAL vldummy;
    REAL vudummy;
    INTEGER ildummy;
    INTEGER iudummy;
    REAL abstol = Rlamch("U");
    REAL PI;
    PI = pi(PI);

    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = n - std::max(i, j) + 1;
        }
    }
    printf("a ="); printmat(n, n, a, n); printf("\n");

    // work space query
    lwork = -1;
    REAL *work = new REAL[1];
    liwork = -1;
    INTEGER *iwork = new INTEGER[1];

    Rsyevr("N", "A", "U", n, a, n, vldummy, vudummy, ildummy, iudummy, abstol, m, w, z, n, isuppz, work, lwork, iwork, liwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new REAL[std::max((INTEGER)1, lwork)];
    liwork = iwork[0];
    delete[] iwork;
    iwork = new INTEGER[std::max((INTEGER)1, liwork)];

    // diagonalize matrix
    Rsyevr("N", "A", "U", n, a, n, vldummy, vudummy, ildummy, iudummy, abstol, m, w, z, n, isuppz, work, lwork, iwork, liwork, info);

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

    delete[] iwork;
    delete[] work;
    delete[] reldiff;
    delete[] lambda;
    delete[] w;
    delete[] isuppz;
    delete[] z;
    delete[] a;
}

int main(int argc, char *argv[]) {
    int STARTN = 100;
    int ENDN = 1000;
    int STEPN = 100;
    if (argc != 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STARTN", argv[i]) == 0) {
                STARTN = atoi(argv[++i]);
            } else if (strcmp("-ENDN", argv[i]) == 0) {
                ENDN = atoi(argv[++i]);
            }
        }
    }
    for (int n = STARTN; n <= ENDN; n = n + STEPN) {
        printf("# Eigenvalues of Frank matrix of order n=%d\n", n);
        Frank((INTEGER)n);
    }
}
