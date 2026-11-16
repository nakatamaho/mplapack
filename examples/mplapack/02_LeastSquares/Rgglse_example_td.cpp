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

int main() {
    mplapackint m = 3, n = 2, p = 1, lda = m, ldb = p, info, lwork = -1;
    td_real *a = new td_real[lda * n];
    td_real *bmat = new td_real[ldb * n];
    td_real *c = new td_real[m];
    td_real *d = new td_real[p];
    td_real *x = new td_real[n];
    td_real *xexact = new td_real[n];
    a[0] = 1.0;
    a[1] = 0.0;
    a[2] = 1.0;
    a[0 + lda] = 0.0;
    a[1 + lda] = 1.0;
    a[2 + lda] = 1.0;
    bmat[0] = 1.0;
    bmat[0 + ldb] = 1.0;
    xexact[0] = 1.0;
    xexact[1] = 2.0;
    for (mplapackint i = 0; i < m; i++)
        c[i] = a[i] * xexact[0] + a[i + lda] * xexact[1];
    d[0] = 3.0;
    td_real wk;
    Rgglse(m, n, p, a, lda, bmat, ldb, c, d, x, &wk, lwork, info);
    lwork = castINTEGER_td(wk);
    td_real *work = new td_real[lwork];
    Rgglse(m, n, p, a, lda, bmat, ldb, c, d, x, work, lwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("constraint B*x-d = "); printnum(x[0] + x[1] - d[0]); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, (mplapackint)1, x, n, xexact, n)); printf("\n");
    delete[] work;
    delete[] xexact;
    delete[] x;
    delete[] d;
    delete[] c;
    delete[] bmat;
    delete[] a;
    return info != 0 ? 1 : 0;
}
