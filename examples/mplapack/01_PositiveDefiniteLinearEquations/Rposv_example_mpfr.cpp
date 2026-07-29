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

mpfr_class maxabs(mpfr_class a, mpfr_class b) {
    mpfr_class d = abs(a - b);
    return d;
}

mpfr_class max_solution_error(mplapackint n, mplapackint nrhs, mpfr_class *x, mplapackint ldx, mpfr_class *xexact, mplapackint ldxexact) {
    mpfr_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpfr_class d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpfr_class max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpfr_class *a, mplapackint lda, mpfr_class *x, mplapackint ldx, mpfr_class *b, mplapackint ldb) {
    mpfr_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpfr_class s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpfr_class d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 2, nrhs = 1, lda = n, ldb = n, info;
    mpfr_class *a = new mpfr_class[n * n];
    mpfr_class *aorg = new mpfr_class[n * n];
    mpfr_class *b = new mpfr_class[n];
    mpfr_class *borg = new mpfr_class[n];
    mpfr_class *xexact = new mpfr_class[n];
    a[0]=4; a[1]=2; a[2]=2; a[3]=5;
    xexact[0]=1; xexact[1]=2;
    for (mplapackint i=0;i<n*n;i++) aorg[i]=a[i];
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint k = 0; k < n; k++)
            b[i] += aorg[i + k * lda] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    Rposv("L", n, nrhs, a, lda, b, ldb, info);
    printf("x = "); printvec(b,n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n");
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
