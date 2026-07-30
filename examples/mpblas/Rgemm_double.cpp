//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
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
    mplapackint n = 3;

    double *a = new double[n * n];
    double *b = new double[n * n];
    double *c = new double[n * n];
    double alpha, beta;

//setting A matrix
    a[0 + 0 * n] = 1.0;    a[0 + 1 * n] = 8.0;    a[0 + 2 * n] = 3.0;
    a[1 + 0 * n] = 2.0;    a[1 + 1 * n] = 10.0;   a[1 + 2 * n] = 8.0;
    a[2 + 0 * n] = 9.0;    a[2 + 1 * n] = -5.0;   a[2 + 2 * n] = -1.0;

    b[0 + 0 * n] = 9.0;    b[0 + 1 * n] = 8.0;    b[0 + 2 * n] = 3.0;
    b[1 + 0 * n] = 3.0;    b[1 + 1 * n] = -11.0;  b[1 + 2 * n] = 8.0;
    b[2 + 0 * n] = -8.0;   b[2 + 1 * n] = 6.0;    b[2 + 2 * n] = 1.0;

    c[0 + 0 * n] = 3.0;    c[0 + 1 * n] = 3.0;    c[0 + 2 * n] = -9.0;
    c[1 + 0 * n] = 8.0;    c[1 + 1 * n] = 4.0;    c[1 + 2 * n] = 8.0;
    c[2 + 0 * n] = 6.0;    c[2 + 1 * n] = 1.0;    c[2 + 2 * n] = -2.0;

    printf("# Rgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = 3.0;
    beta = -2.0;
    Rgemm("n", "n", n, n, n, alpha, a, n, b, n, beta, c, n);

    printf("alpha = "); printnum(alpha); printf("\n");
    printf("beta = "); printnum(beta); printf("\n");
    printf("ans ="); printmat(n, n, c, n); printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * a * b + beta * c \n");
    delete[]c;
    delete[]b;
    delete[]a;
}
