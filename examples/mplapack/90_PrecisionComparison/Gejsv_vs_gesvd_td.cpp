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
int main() {
    mplapackint m = 8, n = 4, lda = m, info, lwork = -1;
    td_real *a = new td_real[m * n];
    td_real *b = new td_real[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            td_real scale = pow(td_real(10.0), td_real((double)(i - 4)));
            a[i + j * lda] = scale / (td_real((double)(i + j + 1)));
            b[i + j * lda] = a[i + j * lda];
        }
    td_real *s1 = new td_real[n];
    td_real *s2 = new td_real[n];
    td_real *u = new td_real[1];
    td_real *vt = new td_real[1];
    td_real wk;
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_td(wk);
    td_real *work = new td_real[lwork];
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, info);
    delete[] work;
    mplapackint *iwork = new mplapackint[m + 3 * n + 10];
    lwork = 5 * n * n + 9 * n + m;
    if (lwork < 2 * m + n)
        lwork = 2 * m + n;
    if (lwork < 4 * n + n * n)
        lwork = 4 * n + n * n;
    if (lwork < 7)
        lwork = 7;
    work = new td_real[lwork];
    Rgejsv("G", "N", "N", "N", "N", "N", m, n, b, lda, s2, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, iwork, info);
    printf("Rgesvd singular values = "); printvec(s1, n); printf("\n");
    printf("Rgejsv singular values = "); printvec(s2, n); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        td_real rel = abs(s1[i] - s2[i]) / (abs(s2[i]) + Rlamch_td("S"));
        printf("relative difference[%ld] = ", (long)i); printnum(rel); printf("\n");
    }
    delete[] work;
    delete[] iwork;
    delete[] vt;
    delete[] u;
    delete[] s2;
    delete[] s1;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
