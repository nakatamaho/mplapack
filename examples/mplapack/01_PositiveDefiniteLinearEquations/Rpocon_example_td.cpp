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
td_real maxabs(td_real a, td_real b) {
    td_real d = abs(a - b);
    return d;
}

td_real max_solution_error(mplapackint n, mplapackint nrhs, td_real *x, mplapackint ldx, td_real *xexact, mplapackint ldxexact) {
    td_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            td_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

td_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, td_real *a, mplapackint lda, td_real *x, mplapackint ldx, td_real *b, mplapackint ldb) {
    td_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            td_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            td_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

td_real one_norm(mplapackint n, td_real *a, mplapackint lda) {
    td_real anorm = 0.0;
    for (mplapackint j = 0; j < n; j++) {
        td_real s = 0.0;
        for (mplapackint i = 0; i < n; i++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
int main() {
    mplapackint n = 2, lda = n, info;
    td_real *a = new td_real[n * n];
    td_real *aorg = new td_real[n * n];
    td_real *work = new td_real[3 * n];
    mplapackint *iwork = new mplapackint[n];
    td_real rcond = 0.0;
    a[0] = 4.0;
    a[1] = 2.0;
    a[2] = 2.0;
    a[3] = 5.0;
    for (mplapackint i = 0; i < n * n; i++)
        aorg[i] = a[i];
    Rpotrf("L", n, a, lda, info);
    if (info == 0)
        Rpocon("L", n, a, lda, one_norm(n, aorg, lda), rcond, work, iwork, info);
    printf("A = "); printmat(n, n, aorg, lda); printf("\n");
    printf("rcond_1 = "); printnum(rcond); printf("\n");
    delete[] iwork;
    delete[] work;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
