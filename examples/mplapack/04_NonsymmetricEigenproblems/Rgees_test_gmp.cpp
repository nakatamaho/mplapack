//public domain
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>
#include <iostream>
#include <cstring>
#include <algorithm>

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
bool rselect(mpf_class ar, mpf_class ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    mpf_class *a = new mpf_class[n * n];
    mpf_class *vs = new mpf_class[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 3 * n;
    mpf_class *wr = new mpf_class[n];
    mpf_class *wi = new mpf_class[n];
    mpf_class *work = new mpf_class[lwork];
    bool bwork[n];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = 0.35;     a[0 + 1 * n] = -0.1160;   a[0 + 2 * n] = -0.3886;    a[0 + 3 * n] = -0.2942;
    a[1 + 0 * n] = -0.5140;  a[1 + 1 * n] = 0.1225;    a[1 + 2 * n] = 0.1004;     a[1 + 3 * n] = 0.1126;
    a[2 + 0 * n] = 0.0;      a[2 + 1 * n] = 0.6443;    a[2 + 2 * n] = -0.1357;    a[2 + 3 * n] = -0.0977;
    a[3 + 0 * n] = -3.0;     a[3 + 1 * n] = 0.0;       a[3 + 2 * n] = 0.4262;     a[3 + 3 * n] = 0.1632;

    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "N", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("vs*vs'\n");
    printf("eig(a)\n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }
    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vs;
    delete[] a;
}
