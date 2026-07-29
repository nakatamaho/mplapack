//public domain
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>

#define GMP_FORMAT "%+68.64Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"

inline void printnum(mpf_class rtmp) { gmp_printf(GMP_FORMAT, rtmp.get_mpf_t()); }
inline void printnum_short(mpf_class rtmp) { gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t()); }

//Matlab/Octave format
void printvec(mpf_class *a, int len) {
    mpf_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpf_class * a, int lda)
{
    mpf_class mtmp;

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
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, info, lwork = -1;
    mpf_class *a = new mpf_class[n * n];
    mpf_class *b = new mpf_class[n * n];
    mpf_class *alphar = new mpf_class[n];
    mpf_class *alphai = new mpf_class[n];
    mpf_class *beta = new mpf_class[n];
    mpf_class *vl = new mpf_class[1];
    mpf_class *vr = new mpf_class[1];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0;
        b[i] = 0;
    }
    a[0] = 1;
    a[4] = 2;
    a[8] = 3;
    b[0] = 1;
    b[4] = 1;
    b[8] = 0;
    mpf_class wk;
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_gmp(wk);
    mpf_class *work = new mpf_class[lwork];
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, work, lwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alphar[i]); printf(" + "); printnum(alphai[i]); printf("i, beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_gmp("E")) {
            printf(", lambda = Inf\n");
        } else {
            printf(", lambda = "); printnum(alphar[i] / beta[i]); printf(" + "); printnum(alphai[i] / beta[i]); printf("i\n");
        }
    }
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] beta;
    delete[] alphai;
    delete[] alphar;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
