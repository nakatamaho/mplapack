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
//taking from Collection of Matrices for Testing Computational Algorithms 1969 Robert T. Gregory, David L. Karney pp.30
int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    mpc_class *a = new mpc_class[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix


    a[0 + 0 * n] = mpc_class(1.0, 0.0);   a[0 + 1 * n] = mpc_class(1.0, 2.0);    a[0 + 2 * n] = mpc_class(2.0, 10.0);
    a[1 + 0 * n] = mpc_class(1.0, 1.0);   a[1 + 1 * n] = mpc_class(0.0, 3.0);    a[1 + 2 * n] = mpc_class(-5.0, 14.0);
    a[2 + 0 * n] = mpc_class(1.0, 1.0);   a[2 + 1 * n] = mpc_class(0.0, 5.0);    a[2 + 2 * n] = mpc_class(-8.0, 20.0);

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    mpc_class *work = new mpc_class[1];

    Cgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_gmp (work[0].real());
    delete[]work;
    work = new mpc_class[std::max(1, (int) lwork)];

//inverse matrix
    Cgetrf(n, n, a, n, ipiv, info);
    Cgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
