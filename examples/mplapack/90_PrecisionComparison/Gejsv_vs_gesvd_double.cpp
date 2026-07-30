//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }

// Matlab/Octave format
void printvec(double *a, int len) {
    double tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, double *a, int lda)
{
    double mtmp;

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
    double *a = new double[m * n];
    double *b = new double[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            double scale = pow(double(10.0), double((double)(i - 4)));
            a[i + j * lda] = scale / (double((double)(i + j + 1)));
            b[i + j * lda] = a[i + j * lda];
        }
    double *s1 = new double[n];
    double *s2 = new double[n];
    double *u = new double[1];
    double *vt = new double[1];
    double wk;
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
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
    work = new double[lwork];
    Rgejsv("G", "N", "N", "N", "N", "N", m, n, b, lda, s2, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, iwork, info);
    printf("Rgesvd singular values = "); printvec(s1, n); printf("\n");
    printf("Rgejsv singular values = "); printvec(s2, n); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        double rel = abs(s1[i] - s2[i]) / (abs(s2[i]) + Rlamch_double("S"));
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
