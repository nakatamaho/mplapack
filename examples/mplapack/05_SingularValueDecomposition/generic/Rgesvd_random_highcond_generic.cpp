#include <random>

int compare_real(const void *a, const void *b)
{
    return *(REAL*)a > *(REAL*)b;
}

int main(int argc, char *argv[]) {
    INTEGER n = 3;
    INTEGER dispersion = 3;
    if (argc != 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp("-DIMN", argv[i]) == 0) {
                n = atoi(argv[++i]);
            }
            if (strcmp("-DISPERSION", argv[i]) == 0) {
                dispersion = atoi(argv[++i]);
            }
        }
    }
    printf("#dimension %d, dispersion = %d \n", (int)n, (int)dispersion);

    REAL *a = new REAL[n * n];
    REAL *aorg = new REAL[n * n];
    REAL *ainv = new REAL[n * n];
    REAL *at_a = new REAL[n * n];
    REAL *I_ = new REAL[n * n]; //I is reserved for imaginary number 
    REAL *s = new REAL[n * n];
    REAL *sorg = new REAL[n];
    REAL *u = new REAL[n * n];
    REAL *vt = new REAL[n * n];
    REAL *w = new REAL[n * n];

    INTEGER lwork = std::max({(INTEGER)1, 3 * n + n, 5 * n});
    INTEGER liwork;
    INTEGER *ipiv = new INTEGER[n];
    INTEGER info;

    REAL *work = new REAL[lwork];

    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::normal_distribution<> dist(0, (double)dispersion);
    std::uniform_int_distribution<> dist1(-n*n, n*n);

    // Generation of high condition integer matrix.
    // Strategy. Generate a matrix whose elements are integer with determinant = 1.
    // Then the elements of inverse of the matrix is also integer.
    // Prepare a geometric integer series, and obtain a matrix with diagonal elements by the series.
    // Calculate A-1 S A.
    // 1. Set Hessenberg matirx with a[1,n] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j > i) {
                int r = dist(engine);
                a[i + j * n] = r;
                aorg[i + j * n] = r;
            } else if (j == i) {
                int r = dist(engine);
                while (r == 0) {
                    r = dist(engine);
                }
                a[i + j * n] = r;
                aorg[i + j * n] = r;
            } else if (j == i - 1) {
                a[i + j * n] = 1.0;
                aorg[i + j * n] = 1.0;
            } else {
                a[i + j * n] = 0.0;
                aorg[i + j * n] = 0.0;
            }
        }
    }
    a[0 + (n - 1) * n] = 0.0;
    aorg[0 + (n - 1) * n] = 0.0;
    printf("split_long_rows(0)\n");
    printf("aorg ="); printmat(n, n, aorg, n); printf("\n");

    // 2. get determinant via LU factorization
    Rgetrf(n, n, a, n, ipiv, info);
    // printf("aLU ="); printmat(n, n, a, n); printf("\n");
    REAL det = 1;
    for (int i = 0; i < n; i++) {
        det = det * a[i + i * n];
        if (ipiv[i] != i + 1)
            det = det * -1.0;
    }
    printf("det="); printnum(det); printf("\n");

    // 3. Set Hessenberg matirx with a[1,n] = 0 to make det of matrix a = 1.
    aorg[0 + (n - 1) * n] = -det + 1;
    printf("anew ="); printmat(n, n, aorg, n); printf("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            ainv[i + j * n] = aorg[i + j * n];
            a[i + j * n] = aorg[i + j * n];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s[i + j * n] = 0.0;
        }
    }
    // 4. genrarate eigenvalues
    for (int i = 0; i < n; i++) {
         int r = dist1(engine);
         s[i + i * n] =  r;
         while ( r == 0 ) {
             r = dist1(engine);
             s[i + i * n] = r;
	 }
    }
    s[0] = 1.0;
    for (int i = 0; i < n; i++)
        sorg[i] = s[i + i * n];
    printf("s = ["); for (int i = 0; i < n; i++) { printnum(s[i + i * n]);  printf(" "); } printf("]\n");

    // 5. inverse matrix. All the elements are integers.
    Rgetrf(n, n, ainv, n, ipiv, info);
    Rgetri(n, ainv, n, ipiv, work, lwork, info);

    // 5.5. verify Ainv * A = I
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            REAL rtmp = 0.0;
            for (int k = 0; k < n; k++) {
                rtmp = rtmp + ainv[i + k * n] * a[k + j * n];
            }
            I_[i + j * n] = rtmp;
        }
    }
    printf("I ="); printmat(n, n, I_, n); printf("\n");

    // 6. Make a  A <- Ainv * S * A
    printf("ainv ="); printmat(n, n, ainv, n); printf("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            REAL rtmp = 0.0;
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    rtmp = rtmp + ainv[i + k * n] * s[k + l * n] * aorg[l + j * n];
                }
            }
            a[i + j * n] = rtmp;
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aorg[i + j * n] = a[i + j * n];
        }
    }
    // 7. svd(A)
    Rgesvd("A", "A", n, n, a, n, s, u, n, vt, n, work, lwork, info);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i + j * n] = aorg[i + j * n];
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            REAL rtmp = 0.0;
            for (int k = 0; k < n; k++) {
                rtmp = rtmp + a[k + i * n] * a[k + j * n];
            }
            at_a[i + j * n] = rtmp;
        }
    }
    // 7. eig(A^t A).
    printf("at_a ="); printmat(n, n, at_a, n); printf("\n");
    printf("eig(at_a)\n");
    REAL rtmp = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ( abs(at_a[i + j * n]) > abs(rtmp) )
                rtmp = at_a[i + j * n];
	}
    }
    printf("abs_of_max_at_a="); printnum(rtmp); printf("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ( abs(at_a[i + j * n]) < abs(rtmp) )
                rtmp = at_a[i + j * n];
	}
    }
    printf("abs_min_at_a="); printnum(rtmp); printf("\n");
    delete[] work;

    // work space query
    lwork = -1;
    work = new REAL[1];
    liwork = -1;
    INTEGER *iwork = new INTEGER[1];

    Rsyevd("N", "U", n, at_a, n, w, work, lwork, iwork, liwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new REAL[std::max((INTEGER)1, lwork)];
    liwork = iwork[0];
    delete[] iwork;
    iwork = new INTEGER[std::max((INTEGER)1, liwork)];

    // diagonalize matrix
    Rsyevd("N", "U", n, at_a, n, w, work, lwork, iwork, liwork, info);

    qsort(s, n, sizeof(REAL), compare_real);
    printf("s=[");
    for (int i = 0; i < n; i++) { printnum(s[i]); printf(" "); } printf(" ] \n");
    printf("s_squared=["); for (int i = 0; i < n; i++) { printnum(s[i] * s[i]); printf(" "); } printf(" ] \n");
    printf("      w = ["); for (int i = 0; i < n; i++) { printnum(w[i]); printf(" "); } printf(" ] \n");

    // 8. There is a relation \lambda_i of eig(A^t A) and \sigma_i svd(A)
    // \lambda_i = \sigma_i^2
    // 9. Relative error

    REAL relerror;
    for (int i = 0; i < n; i = i + 1) {
        relerror = abs ( (w[i] - s[i] * s[i]) / (s[i] * s[i]) ) ;
        printf("Relative_error_%d = ", (int)i); printnum(relerror); printf("\n");
    }

    delete[] work;
    delete[] ipiv;
    delete[] w;
    delete[] u;
    delete[] sorg;
    delete[] s;
    delete[] I_;
    delete[] at_a;
    delete[] ainv;
    delete[] aorg;
    delete[] a;
}
