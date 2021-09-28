//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(qd_real *a, int len) {
    qd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, qd_real * a, int lda)
{
    qd_real mtmp;
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
    mplapackint info;
    qd_real *a = new qd_real[n * n];
    qd_real *b = new qd_real[n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix
///Collection of Matrices for Testing Computational Algorithms 1969/1/21 Robert T. Gregory,  David L. Karney
    a[0 + 0 * n] = 1;   a[0 + 1 * n] = -2;  a[0 + 2 * n] = 3;    a[0 + 3 * n] = 1;
    a[1 + 0 * n] = -2;  a[1 + 1 * n] = 1;   a[1 + 2 * n] = -2;   a[1 + 3 * n] = -1;
    a[2 + 0 * n] = 3;   a[2 + 1 * n] = -2;  a[2 + 2 * n] = 1;    a[2 + 3 * n] = 5;
    a[3 + 0 * n] = 1;   a[3 + 1 * n] = -1;  a[3 + 2 * n] = 5;    a[3 + 3 * n] = 3;

    b[0] = 3;
    b[1] = -4;
    b[2] = 7;
    b[3] = 8;

//answer is [1, 1, 1, 1]
    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printvec(b, n); printf("\n");

//Solve linear equation
    Rgesv(n, (mplapackint)1, a, n, ipiv, b, n, info);

    printf("x ="); printvec(b, n); printf("\n");
    printf("a*x'\n");
    delete[]ipiv;
    delete[]b;
    delete[]a;
}
