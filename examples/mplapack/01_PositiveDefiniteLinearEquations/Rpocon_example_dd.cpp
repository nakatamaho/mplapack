//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_dd.h>
#include <mplapack_dd.h>

#define DD_PRECISION_SHORT 16

inline void printnum(dd_real rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(dd_real *a, int len) {
    dd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, dd_real * a, int lda)
{
    dd_real mtmp;
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
dd_real maxabs(dd_real a, dd_real b) {
    dd_real d = abs(a - b);
    return d;
}

dd_real max_solution_error(mplapackint n, mplapackint nrhs, dd_real *x, mplapackint ldx, dd_real *xexact, mplapackint ldxexact) {
    dd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            dd_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

dd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, dd_real *a, mplapackint lda, dd_real *x, mplapackint ldx, dd_real *b, mplapackint ldb) {
    dd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            dd_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            dd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

dd_real one_norm(mplapackint n, dd_real *a, mplapackint lda) {
    dd_real anorm = 0.0;
    for (mplapackint j = 0; j < n; j++) {
        dd_real s = 0.0;
        for (mplapackint i = 0; i < n; i++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
int main() {
    mplapackint n = 2, lda = n, info;
    dd_real *a = new dd_real[n * n];
    dd_real *aorg = new dd_real[n * n];
    dd_real *work = new dd_real[3 * n];
    mplapackint *iwork = new mplapackint[n];
    dd_real rcond = 0.0;
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
