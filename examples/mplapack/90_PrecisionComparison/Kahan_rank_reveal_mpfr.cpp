//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpfr_class rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpfr_class rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpfr_class *a, int len) {
    mpfr_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpfr_class *a, int lda) {
    mpfr_class mtmp;
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
    mpfr_class theta = mpfr_class(0.1), c = cos(theta), s = sin(theta);
    mpfr_class *a = new mpfr_class[m * n];
    mpfr_class *asvd = new mpfr_class[m * n];
    mpfr_class *aqr = new mpfr_class[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            mpfr_class scale = pow(s, mpfr_class(i));
            mpfr_class val = (i == j) ? mpfr_class(1.0) : ((i < j) ? -c : mpfr_class(0.0));
            a[i + j * lda] = scale * val;
            asvd[i + j * lda] = a[i + j * lda];
            aqr[i + j * lda] = a[i + j * lda];
        }
    mpfr_class *sigma = new mpfr_class[n];
    mpfr_class *u = new mpfr_class[1];
    mpfr_class *vt = new mpfr_class[1];
    mpfr_class wk;
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_mpfr(wk);
    mpfr_class *work = new mpfr_class[lwork];
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, info);
    delete[] work;
    mplapackint *jpvt = new mplapackint[n];
    mpfr_class *tau = new mpfr_class[n];
    for (mplapackint i = 0; i < n; i++)
        jpvt[i] = 0;
    lwork = -1;
    Rgeqp3(m, n, aqr, lda, jpvt, tau, &wk, lwork, info);
    lwork = castINTEGER_mpfr(wk);
    work = new mpfr_class[lwork];
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
