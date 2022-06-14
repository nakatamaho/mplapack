int main() {
    INTEGER n = 4;
    INTEGER m = 4;

    COMPLEX *a = new COMPLEX[m * n];
    REAL *s = new REAL[std::min(m, n)];
    COMPLEX *u = new COMPLEX[m * m];
    COMPLEX *vt = new COMPLEX[n * n];
    INTEGER lwork = std::max((INTEGER)1, 2 * std::min(m, n) + std::max(m, n));
    COMPLEX *work = new COMPLEX[lwork];
    REAL *rwork = new REAL[5 * std::min(m, n)];
    INTEGER info;

    // setting A matrix
    a[0 + 0 * n] = COMPLEX(0.9, -1.0); a[0 + 1 * n] = COMPLEX(20.0, -2.25);  a[0 + 2 * n] = COMPLEX(1.75, -0.5);  a[0 + 3 * n] = COMPLEX(0.0, 0.5);
    a[1 + 0 * n] = COMPLEX(8.0,-2.25); a[1 + 1 * n] = COMPLEX(-0.25, 0.0);   a[1 + 2 * n] = COMPLEX(1.25, -0.25); a[1 + 3 * n] = COMPLEX(-3.75, 0.0);
    a[2 + 0 * n] = COMPLEX(-1.75,0.0); a[2 + 1 * n] = COMPLEX(-80.0,  1.25); a[2 + 2 * n] = COMPLEX(1.5, 0.0);    a[2 + 3 * n] = COMPLEX(30.0, 2.25);
    a[3 + 0 * n] = COMPLEX(3.0, 0.25); a[3 + 1 * n] = COMPLEX(1.75, 0.0);    a[3 + 2 * n] = COMPLEX(0.0, 2.25);   a[3 + 3 * n] = COMPLEX(-0.25, -80.0);
    
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Cgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, rwork, info);
    printf("s="); printvec(s, std::min(m, n)); printf("\n");
    if (m < n)
        printf("padding=zeros(%d, %d-%d)\n", (int)m, (int)n, (int)m);
    if (n < m)
        printf("padding=zeros(%d-%d,%d)\n", (int)m, (int)n, (int)n);
    printf("u ="); printmat(m, m, u, m); printf("\n");
    printf("vt ="); printmat(n, n, vt, n); printf("\n");
    printf("svd(a)\n");
    if (m < n)
        printf("sigma=[diag(s) padding] \n");
    if (n < m)
        printf("sigma=[diag(s); padding] \n");
    if (n == m)
        printf("sigma=[diag(s)] \n");
    printf("sigma \n");
    printf("u * sigma  * vt\n");
    delete[] rwork;
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
