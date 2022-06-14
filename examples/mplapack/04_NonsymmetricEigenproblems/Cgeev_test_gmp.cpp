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
bool rselect(mpf_class ar, mpf_class ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    mpc_class *a = new mpc_class[n * n];
    mpc_class *w = new mpc_class[n];
    mpc_class *vl = new mpc_class[n * n];
    mpc_class *vr = new mpc_class[n * n];
    mplapackint lwork = 4 * n;
    mpc_class *work = new mpc_class[lwork];    
    mpf_class *rwork = new mpf_class[lwork];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = mpc_class(5.0, 9.0); a[0 + 1 * n] = mpc_class(5.0, 5.0);   a[0 + 2 * n] = mpc_class(-6.0, -6.0); a[0 + 3 * n] = mpc_class(-7.0, -7.0);
    a[1 + 0 * n] = mpc_class(3.0, 3.0); a[1 + 1 * n] = mpc_class(6.0, 10.0);  a[1 + 2 * n] = mpc_class(-5.0, -5.0); a[1 + 3 * n] = mpc_class(-6.0, -6.0);
    a[2 + 0 * n] = mpc_class(2.0, 2.0); a[2 + 1 * n] = mpc_class(3.0, 3.0);   a[2 + 2 * n] = mpc_class(-1.0,  3.0); a[2 + 3 * n] = mpc_class(-5.0, -5.0);
    a[3 + 0 * n] = mpc_class(1.0, 1.0); a[3 + 1 * n] = mpc_class(2.0, 2.0);   a[3 + 2 * n] = mpc_class(-3.0, -3.0); a[3 + 3 * n] = mpc_class(0.0, 4.0); 

    printf("# Ex. 6.5 p. 116, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgeev("V", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
    printf("lambda ="); printvec(w,n); printf("\n");    
    printf("vr ="); printmat(n,n,vr,n); printf("\n");    

    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
}
