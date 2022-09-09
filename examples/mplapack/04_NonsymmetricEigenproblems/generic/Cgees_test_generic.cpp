bool cselect(COMPLEX a) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    INTEGER n = 4;

    COMPLEX *a = new COMPLEX[n * n];
    INTEGER sdim = 0;
    INTEGER lwork = 2 * n;
    COMPLEX *w = new COMPLEX[n];
    COMPLEX *vs = new COMPLEX[n * n];
    COMPLEX *work = new COMPLEX[lwork];
    REAL *rwork = new REAL[n];
    bool bwork[n];
    INTEGER info;

    // setting A matrix
    a[0 + 0 * n] = COMPLEX(5.0,  9.0); a[0 + 1 * n] = COMPLEX(5.0, 5.0);   a[0 + 2 * n] = COMPLEX(-6.0, -6.0); a[0 + 3 * n] = COMPLEX(-7.0,-7.0);
    a[1 + 0 * n] = COMPLEX(3.0,  3.0); a[1 + 1 * n] = COMPLEX(6.0,10.0);   a[1 + 2 * n] = COMPLEX(-5.0, -5.0); a[1 + 3 * n] = COMPLEX(-6.0,-6.0);
    a[2 + 0 * n] = COMPLEX(2.0,  2.0); a[2 + 1 * n] = COMPLEX(3.0, 3.0);   a[2 + 2 * n] = COMPLEX(-1.0, 3.0);  a[2 + 3 * n] = COMPLEX(-5.0,-5.0);
    a[3 + 0 * n] = COMPLEX(1.0,  1.0); a[3 + 1 * n] = COMPLEX(2.0, 2.0);   a[3 + 2 * n] = COMPLEX(-3.0,-3.0);  a[3 + 3 * n] = COMPLEX(0.0, 4.0);

    printf("# Ex. 6.5 p. 116, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgees("V", "S", cselect, n, a, n, sdim, w, vs, n, work, lwork, rwork, bwork, info);
    printf("w ="); printvec(w, n); printf("\n");
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
    printf("eig(a)\n");

    delete[] rwork;
    delete[] work;
    delete[] vs;
    delete[] w;
    delete[] a;
}
