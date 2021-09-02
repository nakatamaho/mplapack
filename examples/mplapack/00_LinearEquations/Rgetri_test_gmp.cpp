//public domain
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>
#include <iostream>
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
int main()
{
    mplapackint n = 4;
    mplapackint lwork, info;

    mpf_class *a = new mpf_class[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
//https://academic.oup.com/qjmam/article/1/1/253/1883468
//The Quarterly Journal of Mechanics and Applied Mathematics, Volume 1, Issue 1, 1948, Pages 253280, https://doi.org/10.1093/qjmam/1.1.253
//https://babel.hathitrust.org/cgi/pt?id=mdp.39015023899019&view=1up&seq=30
//Marcus, M. (1960). Basic theorems in matrix theory. Washington: U.S. Govt. Print. Off.

    a[0 + 0 * n] = 5;    a[0 + 1 * n] = 7;    a[0 + 2 * n] = 6;      a[0 + 3 * n] = 5;
    a[1 + 0 * n] = 7;    a[1 + 1 * n] = 10;   a[1 + 2 * n] = 8;      a[1 + 3 * n] = 7;
    a[2 + 0 * n] = 6;    a[2 + 1 * n] = 8;    a[2 + 2 * n] = 10;     a[2 + 3 * n] = 9;
    a[3 + 0 * n] = 5;    a[3 + 1 * n] = 7;    a[3 + 2 * n] = 9;      a[3 + 3 * n] = 10;

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    mpf_class *work = new mpf_class[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_gmp (work[0]);
    delete[]work;
    work = new mpf_class[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
