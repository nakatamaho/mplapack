//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

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

int main() {
    mplapackint n = 4;

    qd_real *a = new qd_real[n * n];
    qd_real *vs = new qd_real[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 3 * n;
    qd_real *wr = new qd_real[n];
    qd_real *wi = new qd_real[n];
    qd_real *work = new qd_real[lwork];
    bool bwork[n];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = 0.35;     a[0 + 1 * n] = -0.1160;   a[0 + 2 * n] = -0.3886;    a[0 + 3 * n] = -0.2942;
    a[1 + 0 * n] = -0.5140;  a[1 + 1 * n] = 0.1225;    a[1 + 2 * n] = 0.1004;     a[1 + 3 * n] = 0.1126;
    a[2 + 0 * n] = 0.0;      a[2 + 1 * n] = 0.6443;    a[2 + 2 * n] = -0.1357;    a[2 + 3 * n] = -0.0977;
    a[3 + 0 * n] = -3.0;     a[3 + 1 * n] = 0.0;       a[3 + 2 * n] = 0.4262;     a[3 + 3 * n] = 0.1632;

    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "N", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("vs*vs'\n");
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
