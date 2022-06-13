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
inline void printnum(mpc_class ctmp) { gmp_printf(GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t()); }

//Matlab/Octave format
template <class X> void printvec(X *a, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

template <class X> void printmat(int n, int m, X *a, int lda)
{
    X mtmp;

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
    mplapackint n = 4;
    mplapackint m = 4;

    mpc_class *a = new mpc_class[m * n];
    mpf_class *s = new mpf_class[std::min(m, n)];
    mpc_class *u = new mpc_class[m * m];
    mpc_class *vt = new mpc_class[n * n];
    mplapackint lwork = std::max((mplapackint)1, 2 * std::min(m, n) + std::max(m, n));
    mpc_class *work = new mpc_class[lwork];
    mpf_class *rwork = new mpf_class[5 * std::min(m, n)];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = mpc_class(0.25, 0.0); a[0 + 1 * n] = mpc_class(0.0, -2.25); a[0 + 2 * n] = mpc_class(-1.75,0.0);   a[0 + 3 * n] = mpc_class(0.0, 0.25);
    a[1 + 0 * n] = mpc_class(0.0,-2.25); a[1 + 1 * n] = mpc_class(-0.25, 0.0); a[1 + 2 * n] = mpc_class(0.0, -0.25);  a[1 + 3 * n] = mpc_class(-1.75, 0.0);
    a[2 + 0 * n] = mpc_class(-1.75,0.0); a[2 + 1 * n] = mpc_class(0.0, -0.25); a[2 + 2 * n] = mpc_class(0.25, 0.0);   a[2 + 3 * n] = mpc_class(0.0, 2.25);
    a[3 + 0 * n] = mpc_class(0.0, 0.25); a[3 + 1 * n] = mpc_class(-1.75, 0.0); a[3 + 2 * n] = mpc_class(0.0, 2.25);   a[3 + 3 * n] = mpc_class(-0.25, 0.0);
    
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Cgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, rwork, info);
    printf("s="); printvec(s, std::min(m, n)); printf("\n");
    if (m < n)
        printf("padding=zeros(%d, %d-%d)\n", (int)m, (int)n, (int)m);
    if (n < m)
        printf("padding=zeros(%d-%d,%d)\n", (int)m, (int)n, (int)n);
    printf("u ="); printmat(m, m, u, m); printf("\n");
    printf("vt ="); printmat(n, n, vt, n); printf("\n");
    printf("svd(a)\n");
    if (m < n)
        printf("sigma=[diag(s) padding] \n");
    if (n < m)
        printf("sigma=[diag(s); padding] \n");
    if (n == m)
        printf("sigma=[diag(s)] \n");
    printf("sigma \n");
    printf("u * sigma  * vt\n");
    delete[] rwork;
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
