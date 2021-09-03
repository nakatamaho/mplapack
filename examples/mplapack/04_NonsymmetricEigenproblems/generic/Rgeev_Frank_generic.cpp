void Frank(INTEGER n) {
    INTEGER lwork, liwork, info, m;
    REAL *a = new REAL[n * n];
    REAL *vl = new REAL[n * n]; //not used
    REAL *vr = new REAL[n * n]; //not used
    REAL *wr = new REAL[n];
    REAL *wi = new REAL[n];

    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = 0;
	}
    }
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = n - std::max(i, j) + 1;
        }
    }
    for (int i = 1; i <= n - 1 ; i++) {
        a[i + (i - 1) * n] = n - i;
    }
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");

    // work space query
    lwork = -1;
    REAL *work = new REAL[1];
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new REAL[std::max((INTEGER)1, lwork)];

    // diagonalize matrix
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);

    // print out
    printf("#eigenvalues \n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }

    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vr;
    delete[] vl;
    delete[] a;
}

bool rselect(REAL ar, REAL ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main(int argc, char *argv[]) {
    int STARTN = 5;
    int ENDN = 10;
    int STEPN = 1;
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
