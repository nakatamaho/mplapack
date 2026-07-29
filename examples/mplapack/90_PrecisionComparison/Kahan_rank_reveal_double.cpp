//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
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
    mplapackint n = 12, m = n, lda = m, info, lwork = -1;
    double theta = double(0.1), c = cos(theta), s = sin(theta);
    double *a = new double[m * n];
    double *asvd = new double[m * n];
    double *aqr = new double[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            double scale = pow(s, double(i));
            double val = (i == j) ? double(1.0) : ((i < j) ? -c : double(0.0));
            a[i + j * lda] = scale * val;
            asvd[i + j * lda] = a[i + j * lda];
            aqr[i + j * lda] = a[i + j * lda];
        }
    double *sigma = new double[n];
    double *u = new double[1];
    double *vt = new double[1];
    double wk;
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, info);
    delete[] work;
    mplapackint *jpvt = new mplapackint[n];
    double *tau = new double[n];
    for (mplapackint i = 0; i < n; i++)
        jpvt[i] = 0;
    lwork = -1;
    Rgeqp3(m, n, aqr, lda, jpvt, tau, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    work = new double[lwork];
    Rgeqp3(m, n, aqr, lda, jpvt, tau, work, lwork, info);
    printf("smallest singular value = "); printnum(sigma[n - 1]); printf("\n");
    printf("|R(n,n)| from Rgeqp3 = "); printnum(abs(aqr[n - 1 + (n - 1) * lda])); printf("\n");
    delete[] work;
    delete[] tau;
    delete[] jpvt;
    delete[] vt;
    delete[] u;
    delete[] sigma;
    delete[] aqr;
    delete[] asvd;
    delete[] a;
    return info != 0 ? 1 : 0;
}
