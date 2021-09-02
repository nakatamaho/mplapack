//public domain
#include <iostream>
#include <string>
#include <sstream>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(qd_real *a, int len) {
    qd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, qd_real * a, int lda)
{
    qd_real mtmp;
    printf("[ ");
    for (int i = 0; i < n; i++) {
        printf("[ ");
        for (int j = 0; j < m; j++) {
            mtmp = a[i + j * lda];
            printnum(mtmp);     
            if (j < m - 1)
                printf(", ");
        }
        if (i < n - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}
bool rselect(qd_real ar, qd_real ai) {
    // sorting rule for eigenvalues.
    return false;
}

using namespace std;

int main() {
    mplapackint n;

    string str;
    getline(cin, str);
    cout << str << endl;
    getline(cin, str);
    stringstream ss(str);
    ss >> n;
    printf("# n %d\n", (int)n);

    qd_real *a = new qd_real[n * n];
    qd_real *vl = new qd_real[n * n];
    qd_real *vr = new qd_real[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 4 * n;
    qd_real *wr = new qd_real[n];
    qd_real *wi = new qd_real[n];
    qd_real *work = new qd_real[lwork];
    mplapackint info;
    double dtmp;
    for (int i = 0; i < n; i++) {
        getline(cin, str);
        stringstream ss(str);
        for (int j = 0; j < n; j++) {
            ss >> dtmp;
            a[i + j * n] = dtmp;
        }
    }
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgeev("V", "V", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    printf("# left vectors\n");
    printf("vl ="); printmat(n, n, vl, n); printf("\n");
    printf("# left vectors\n");
    printf("vr ="); printmat(n, n, vr, n); printf("\n");
    printf("[vl, d, w] = eig(a)\n");
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
