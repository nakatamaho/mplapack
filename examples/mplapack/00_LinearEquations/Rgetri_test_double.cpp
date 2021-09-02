//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }

// Matlab/Octave format
void printvec(double *a, int len) {
    double tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, double *a, int lda)
{
    double mtmp;

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

    double *a = new double[n * n];
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
    double *work = new double[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_double (work[0]);
    delete[]work;
    work = new double[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("inv_a ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
