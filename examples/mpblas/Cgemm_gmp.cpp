//public domain
#include <mpblas_gmp.h>
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
int main()
{
    mplapackint n = 3;

    mpc_class *a = new mpc_class[n * n];
    mpc_class *b = new mpc_class[n * n];
    mpc_class *c = new mpc_class[n * n];
    mpc_class alpha, beta;

//setting A matrix
    a[0 + 0 * n] = mpc_class(1.0,-1.0);    a[0 + 1 * n] = mpc_class(8.0, 2.2);    a[0 + 2 * n] = mpc_class(0.0, -10.0);
    a[1 + 0 * n] = mpc_class(2.0, 0.0);    a[1 + 1 * n] = mpc_class(10.0,0.0);    a[1 + 2 * n] = mpc_class(8.1, 2.2);
    a[2 + 0 * n] = mpc_class(-9.0,3.0);    a[2 + 1 * n] = mpc_class(-5.0,3.0);    a[2 + 2 * n] = mpc_class(-1.0, 0.0);

    b[0 + 0 * n] = mpc_class(9.0, 0.0);    b[0 + 1 * n] = mpc_class(8.0, -0.01);  b[0 + 2 * n] = mpc_class(3.0, 1.001);
    b[1 + 0 * n] = mpc_class(3.0, -8.0);   b[1 + 1 * n] = mpc_class(-11.0, 0.1);  b[1 + 2 * n] = mpc_class(8.0, 0.00001);
    b[2 + 0 * n] = mpc_class(-8.0, 1.0);   b[2 + 1 * n] = mpc_class(6.0, 0.0);    b[2 + 2 * n] = mpc_class(1.1, 1.0);

    c[0 + 0 * n] = mpc_class(3.0, 1.0);   c[0 + 1 * n] = mpc_class(-3.0, 9.99);   c[0 + 2 * n] = mpc_class(-9.0, -11.0);
    c[1 + 0 * n] = mpc_class(8.0, -1.0);  c[1 + 1 * n] = mpc_class(4.0, 4.44);    c[1 + 2 * n] = mpc_class(8.0, 9.0);
    c[2 + 0 * n] = mpc_class(6.0, 0.0);   c[2 + 1 * n] = mpc_class(-1.0, 0.0);    c[2 + 2 * n] = mpc_class(-2.0, 1.0);

    printf("# Cgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = mpc_class(3.0,-1.2);
    beta = mpc_class(-2.0, -2.0);
    Cgemm("n", "n", n, n, n, alpha, a, n, b, n, beta, c, n);

    printf("alpha = "); printnum(alpha); printf("\n");
    printf("beta = "); printnum(beta); printf("\n");
    printf("ans ="); printmat(n, n, c, n); printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * a * b + beta * c \n");
    delete[]c;
    delete[]b;
    delete[]a;
}
