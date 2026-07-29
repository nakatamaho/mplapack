//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(mplapack_binary128_t rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(mplapack_binary128_t *a, int len) {
    mplapack_binary128_t tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mplapack_binary128_t *a, int lda)
{
    mplapack_binary128_t mtmp;

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
    mplapack_binary128_t theta = mplapack_binary128_t(0.1), c = cos(theta), s = sin(theta);
    mplapack_binary128_t *a = new mplapack_binary128_t[m * n];
    mplapack_binary128_t *asvd = new mplapack_binary128_t[m * n];
    mplapack_binary128_t *aqr = new mplapack_binary128_t[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            mplapack_binary128_t scale = pow(s, mplapack_binary128_t(i));
            mplapack_binary128_t val = (i == j) ? mplapack_binary128_t(1.0) : ((i < j) ? -c : mplapack_binary128_t(0.0));
            a[i + j * lda] = scale * val;
            asvd[i + j * lda] = a[i + j * lda];
            aqr[i + j * lda] = a[i + j * lda];
        }
    mplapack_binary128_t *sigma = new mplapack_binary128_t[n];
    mplapack_binary128_t *u = new mplapack_binary128_t[1];
    mplapack_binary128_t *vt = new mplapack_binary128_t[1];
    mplapack_binary128_t wk;
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_binary128(wk);
    mplapack_binary128_t *work = new mplapack_binary128_t[lwork];
    Rgesvd("N", "N", m, n, asvd, lda, sigma, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, info);
    delete[] work;
    mplapackint *jpvt = new mplapackint[n];
    mplapack_binary128_t *tau = new mplapack_binary128_t[n];
    for (mplapackint i = 0; i < n; i++)
        jpvt[i] = 0;
    lwork = -1;
    Rgeqp3(m, n, aqr, lda, jpvt, tau, &wk, lwork, info);
    lwork = castINTEGER_binary128(wk);
    work = new mplapack_binary128_t[lwork];
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
