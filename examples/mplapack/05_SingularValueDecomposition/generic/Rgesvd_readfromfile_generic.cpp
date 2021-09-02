using namespace std;

int main() {
    INTEGER m, n;

    string str;
    getline(cin, str);
    cout << str << endl;
    getline(cin, str);
    stringstream ss(str);
    ss >> m;
    ss >> n;
    printf("# m n %d %d \n", (int)m, (int)m);

    REAL *a = new REAL[m * n];
    REAL *s = new REAL[std::min(m, n)];
    REAL *u = new REAL[m * m];
    REAL *vt = new REAL[n * n];
    INTEGER lwork = std::max({(INTEGER)1, 3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n)});
    REAL *work = new REAL[lwork];
    INTEGER info;
    double dtmp;
    for (int i = 0; i < m; i++) {
        getline(cin, str);
        stringstream ss(str);
        for (int j = 0; j < n; j++) {
            ss >> dtmp;
            a[i + j * m] = dtmp;
        }
    }
    printf("# octave check\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Rgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, info);
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
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
