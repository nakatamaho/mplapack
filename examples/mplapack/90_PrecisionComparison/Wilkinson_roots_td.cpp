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
void sort_real(mplapackint n, td_real *x) {
    for (mplapackint i = 0; i < n; i++)
        for (mplapackint j = i + 1; j < n; j++)
            if (x[j] < x[i]) {
                td_real t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
}
int main() {
    mplapackint n = 20, lda = n, ldv = 1, info, lwork = -1;
    td_real *coef = new td_real[n + 1];
    for (mplapackint i = 0; i <= n; i++)
        coef[i] = 0.0;
    coef[0] = 1.0;
    for (mplapackint k = 1; k <= n; k++) {
        for (mplapackint j = k; j >= 1; j--)
            coef[j] = coef[j] - td_real((double)k) * coef[j - 1];
    }
    td_real *a = new td_real[n * n];
    for (mplapackint i = 0; i < n * n; i++)
        a[i] = 0.0;
    for (mplapackint i = 1; i < n; i++)
        a[i + (i - 1) * lda] = 1.0;
    for (mplapackint j = 0; j < n; j++)
        a[j + (n - 1) * lda] = -coef[n - j] / coef[0];
    td_real *wr = new td_real[n];
    td_real *wi = new td_real[n];
    td_real *vl = new td_real[1];
    td_real *vr = new td_real[1];
    td_real wk;
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_td(wk);
    td_real *work = new td_real[lwork];
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, work, lwork, info);
    sort_real(n, wr);
    td_real maxerr = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        td_real err = abs(wr[i] - td_real((double)(i + 1)));
        if (maxerr < err)
            maxerr = err;
        printf("root[%ld] = ", (long)i); printnum(wr[i]); printf(", error = "); printnum(err); printf("\n");
    }
    printf("max root error = "); printnum(maxerr); printf("\n");
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] wi;
    delete[] wr;
    delete[] a;
    delete[] coef;
    return info != 0 ? 1 : 0;
}
