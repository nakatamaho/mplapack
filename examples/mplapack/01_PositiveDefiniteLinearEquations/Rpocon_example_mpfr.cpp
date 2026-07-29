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

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpreal *a, int len) {
    mpreal tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpreal *a, int lda) {
    mpreal mtmp;
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

mpreal maxabs(mpreal a, mpreal b) {
    mpreal d = abs(a - b);
    return d;
}

mpreal max_solution_error(mplapackint n, mplapackint nrhs, mpreal *x, mplapackint ldx, mpreal *xexact, mplapackint ldxexact) {
    mpreal err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpreal d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpreal max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpreal *a, mplapackint lda, mpreal *x, mplapackint ldx, mpreal *b, mplapackint ldb) {
    mpreal err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpreal s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpreal d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpreal one_norm(mplapackint n, mpreal *a, mplapackint lda) { mpreal anorm=0.0; for(mplapackint j=0;j<n;j++){ mpreal s=0.0; for(mplapackint i=0;i<n;i++) s=s+abs(a[i+j*lda]); if(anorm<s) anorm=s; } return anorm; }
int main() {
    mplapackint n = 2, lda = n, info;
    mpreal *a = new mpreal[n*n]; mpreal *aorg = new mpreal[n*n]; mpreal *work = new mpreal[3*n]; mplapackint *iwork = new mplapackint[n]; mpreal rcond=0.0;
    a[0]=4; a[1]=2; a[2]=2; a[3]=5; for(mplapackint i=0;i<n*n;i++) aorg[i]=a[i];
    Rpotrf("L", n, a, lda, info);
    if (info == 0) Rpocon("L", n, a, lda, one_norm(n,aorg,lda), rcond, work, iwork, info);
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    printf("rcond_1 = "); printnum(rcond); printf("\n");
    delete[] iwork; delete[] work; delete[] aorg; delete[] a; return info != 0 ? 1 : 0;
}
