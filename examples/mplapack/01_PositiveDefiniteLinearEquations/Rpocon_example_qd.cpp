//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(qd_real *a, int len) {
    qd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, qd_real * a, int lda)
{
    qd_real mtmp;
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
qd_real maxabs(qd_real a, qd_real b) {
    qd_real d = abs(a - b);
    return d;
}

qd_real max_solution_error(mplapackint n, mplapackint nrhs, qd_real *x, mplapackint ldx, qd_real *xexact, mplapackint ldxexact) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            qd_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

qd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, qd_real *a, mplapackint lda, qd_real *x, mplapackint ldx, qd_real *b, mplapackint ldb) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            qd_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            qd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

qd_real one_norm(mplapackint n, qd_real *a, mplapackint lda) { qd_real anorm=0.0; for(mplapackint j=0;j<n;j++){ qd_real s=0.0; for(mplapackint i=0;i<n;i++) s=s+abs(a[i+j*lda]); if(anorm<s) anorm=s; } return anorm; }
int main() {
    mplapackint n = 2, lda = n, info;
    qd_real *a = new qd_real[n*n]; qd_real *aorg = new qd_real[n*n]; qd_real *work = new qd_real[3*n]; mplapackint *iwork = new mplapackint[n]; qd_real rcond=0.0;
    a[0]=4; a[1]=2; a[2]=2; a[3]=5; for(mplapackint i=0;i<n*n;i++) aorg[i]=a[i];
    Rpotrf("L", n, a, lda, info);
    if (info == 0) Rpocon("L", n, a, lda, one_norm(n,aorg,lda), rcond, work, iwork, info);
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    printf("rcond_1 = "); printnum(rcond); printf("\n");
    delete[] iwork; delete[] work; delete[] aorg; delete[] a; return info != 0 ? 1 : 0;
}
