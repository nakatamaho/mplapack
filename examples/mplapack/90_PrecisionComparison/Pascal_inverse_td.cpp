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
td_real binom(mplapackint n, mplapackint k) {
    td_real r = 1.0;
    for (mplapackint i = 1; i <= k; i++)
        r = r * td_real((double)(n - k + i)) / td_real((double)i);
    return r;
}
td_real nearest_integer_error(td_real x) {
    td_real rounded;
    if (x >= td_real(0.0)) {
        rounded = x + td_real(0.5);
    } else {
        rounded = x - td_real(0.5);
    }
    mplapackint nearest = castINTEGER_td(rounded);
    return abs(x - td_real((double)nearest));
}
int main() {
    mplapackint n = 8, lda = n, info, lwork = -1;
    td_real *a = new td_real[n * n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < n; i++)
            a[i + j * lda] = binom(i + j, i);
    Rgetrf(n, n, a, lda, ipiv, info);
    td_real wk;
    if (info == 0)
        Rgetri(n, a, lda, ipiv, &wk, lwork, info);
    lwork = castINTEGER_td(wk);
    td_real *work = new td_real[lwork];
    if (info == 0)
        Rgetri(n, a, lda, ipiv, work, lwork, info);
    td_real err = 0.0;
    for (mplapackint i = 0; i < n * n; i++) {
        td_real d = nearest_integer_error(a[i]);
        if (err < d)
            err = d;
    }
    printf("P inverse = "); printmat(n, n, a, lda); printf("\n");
    printf("max distance to integer = "); printnum(err); printf("\n");
    delete[] work;
    delete[] ipiv;
    delete[] a;
    return info != 0 ? 1 : 0;
}
