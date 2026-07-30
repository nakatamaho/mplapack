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

inline void printnum(mpfr_class rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpfr_class rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpfr_class *a, int len) {
    mpfr_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpfr_class *a, int lda) {
    mpfr_class mtmp;
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

    mpfr_class *a = new mpfr_class[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
//https://academic.oup.com/qjmam/article/1/1/253/1883468
//The Quarterly Journal of Mechanics and Applied Mathematics, Volume 1, Issue 1, 1948, Pages 253280, https://doi.org/10.1093/qjmam/1.1.253
//https://babel.hathitrust.org/cgi/pt?id=mdp.39015023899019&view=1up&seq=30
//Marcus, M. (1960). Basic theorems in matrix theory. Washington: U.S. Govt. Print. Off.

    a[0 + 0 * n] = 5.0;    a[0 + 1 * n] = 7.0;    a[0 + 2 * n] = 6.0;      a[0 + 3 * n] = 5.0;
    a[1 + 0 * n] = 7.0;    a[1 + 1 * n] = 10.0;   a[1 + 2 * n] = 8.0;      a[1 + 3 * n] = 7.0;
    a[2 + 0 * n] = 6.0;    a[2 + 1 * n] = 8.0;    a[2 + 2 * n] = 10.0;     a[2 + 3 * n] = 9.0;
    a[3 + 0 * n] = 5.0;    a[3 + 1 * n] = 7.0;    a[3 + 2 * n] = 9.0;      a[3 + 3 * n] = 10.0;

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    mpfr_class *work = new mpfr_class[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_mpfr (work[0]);
    delete[]work;
    work = new mpfr_class[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    printf("a * ainv - eye(%d)\n", (int)n);
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
