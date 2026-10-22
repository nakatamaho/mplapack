//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_td.h>
#include <mplapack_td.h>

#define TD_PRECISION_SHORT 16

inline void printnum(td_real rtmp) {
    std::cout.precision(TD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(td_real *a, int len) {
    td_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, td_real * a, int lda)
{
    td_real mtmp;
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
bool rselect(td_real ar, td_real ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    td_real *a = new td_real[n * n];
    td_real *vs = new td_real[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 3 * n;
    td_real *wr = new td_real[n];
    td_real *wi = new td_real[n];
    td_real *work = new td_real[lwork];
    bool bwork[n];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = -2.0;     a[0 + 1 * n] = 2.0;   a[0 + 2 * n] = 2.0;    a[0 + 3 * n] = 2.0;
    a[1 + 0 * n] = -3.0;     a[1 + 1 * n] = 3.0;   a[1 + 2 * n] = 2.0;    a[1 + 3 * n] = 2.0;
    a[2 + 0 * n] = -2.0;     a[2 + 1 * n] = 0.0;   a[2 + 2 * n] = 4.0;    a[2 + 3 * n] = 2.0;
    a[3 + 0 * n] = -1.0;     a[3 + 1 * n] = 0.0;   a[3 + 2 * n] = 0.0;    a[3 + 3 * n] = 5.0;

    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "S", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
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
