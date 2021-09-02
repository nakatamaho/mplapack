void inv_hilbert_matrix(int n) {
    INTEGER lwork, info;
    REAL *ainv = new REAL[n * n];
    REAL *aorg = new REAL[n * n];
    REAL *c = new REAL[n * n];
    INTEGER *ipiv = new INTEGER[n];
    REAL one = 1.0, zero = 0.0, mtmp;

    // setting A matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mtmp = (i + 1) + (j + 1) - 1;
            ainv[i + j * n] = one / mtmp;
            aorg[i + j * n] = one / mtmp;
        }
    }

    printf("a = "); printmat(n, n, ainv, n); printf("\n");
    // work space query
    lwork = -1;
    REAL *work = new REAL[1];
    Rgetri(n, ainv, n, ipiv, work, lwork, info);
    lwork = castInTEGER(work[0]);
    delete[] work;
    work = new REAL[std::max(1, (int)lwork)];

    // inverse matrix
    Rgetrf(n, n, ainv, n, ipiv, info);
    Rgetri(n, ainv, n, ipiv, work, lwork, info);
    printf("ainv = "); printmat(n, n, ainv, n); printf("\n");

    // Left residual |ainv * a -I|/|ainv||a|
    // is usually accurate than Right residual |a * ainv -I|/|ainv||a|  See chap.14 of
    // Accuracy and Stability of Numerical Algorithms by Nicholas J. Higham
    // https://doi.org/10.1137/1.9780898718027
    one = 1.0, zero = 0.0;
    Rgemm("N", "N", n, n, n, one, aorg, n, ainv, n, zero, c, n);
    printf("a * ainv ="); printmat(n, n, c, n); printf("\n");
    printf("InfnormR:(a * ainv - I)=");
    mtmp = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (mtmp < abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0)))
                mtmp = abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0));
        }
    }
    printnum(mtmp); printf("\n");

    Rgemm("N", "N", n, n, n, one, ainv, n, aorg, n, zero, c, n);
    printf("ainv * a ="); printmat(n, n, c, n); printf("\n");
    printf("InfnormL:(ainv * a - I)=");
    mtmp = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (mtmp < abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0)))
                mtmp = abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0));
        }
    }
    printnum(mtmp); printf("\n");

    delete[] ainv;
    delete[] aorg;
    delete[] c;
    delete[] ipiv;
    delete[] work;
}

int main()
{
    for (int n = 1; n < 50; n++) {
	printf("# inversion of Hilbert matrix of order n=%d\n", n);
	inv_hilbert_matrix(n);
    }
}
