//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
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
bool rselect(double ar, double ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    double *a = new double[n * n];
    double *vl = new double[n * n];
    double *vr = new double[n * n];
    mplapackint lwork = 4 * n;
    double *wr = new double[n];
    double *wi = new double[n];
    double *work = new double[lwork];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = 0.35;     a[0 + 1 * n] = -0.1160;   a[0 + 2 * n] = -0.3886;    a[0 + 3 * n] = -0.2942;
    a[1 + 0 * n] = -0.5140;  a[1 + 1 * n] = 0.1225;    a[1 + 2 * n] = 0.1004;     a[1 + 3 * n] = 0.1126;
    a[2 + 0 * n] = 0.0;      a[2 + 1 * n] = 0.6443;    a[2 + 2 * n] = -0.1357;    a[2 + 3 * n] = -0.0977;
    a[3 + 0 * n] = -3.0;     a[3 + 1 * n] = 0.0;       a[3 + 2 * n] = 0.4262;     a[3 + 3 * n] = 0.1632;

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
