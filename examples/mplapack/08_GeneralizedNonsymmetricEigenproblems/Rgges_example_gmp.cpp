//public domain
#include <sstream>
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
bool select_none(mpf_class ar, mpf_class ai, mpf_class beta) {
    return false;
}
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, sdim, info, lwork = -1;
    mpf_class *a = new mpf_class[n * n];
    mpf_class *b = new mpf_class[n * n];
    mpf_class *alphar = new mpf_class[n];
    mpf_class *alphai = new mpf_class[n];
    mpf_class *beta = new mpf_class[n];
    mpf_class *vsl = new mpf_class[1];
    mpf_class *vsr = new mpf_class[1];
    bool *bwork = new bool[n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0.0;
        b[i] = 0.0;
    }
    a[0] = 1.0;
    a[4] = 2.0;
    a[8] = 3.0;
    b[0] = 1.0;
    b[4] = 1.0;
    b[8] = 0.0;
    mpf_class wk;
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, &wk, lwork, bwork, info);
    lwork = castINTEGER_gmp(wk);
    mpf_class *work = new mpf_class[lwork];
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, work, lwork, bwork, info);
    printf("S = "); printmat(n, n, a, lda); printf("\n");
    printf("T = "); printmat(n, n, b, ldb); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        printf("lambda[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_gmp("E"))
            printf("Inf\n");
        else {
            printnum(alphar[i] / beta[i]);
            printf(" + "); printnum(alphai[i] / beta[i]); printf("i\n");
        }
    }
    delete[] work;
    delete[] bwork;
    delete[] vsr;
    delete[] vsl;
    delete[] beta;
    delete[] alphai;
    delete[] alphar;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
