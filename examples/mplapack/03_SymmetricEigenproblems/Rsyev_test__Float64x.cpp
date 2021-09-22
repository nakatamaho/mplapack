//public domain
#include <mpblas__Float64x.h>
#include <mplapack__Float64x.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define FLOAT64X_FORMAT "%+25.21Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

void printnum(_Float64x rtmp)
{
    printf(FLOAT64X_FORMAT, rtmp);
    return;
}

//Matlab/Octave format
void printvec(_Float64x *a, int len) {
    _Float64x tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, _Float64x *a, int lda)
{
    _Float64x mtmp;

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

    _Float64x *A = new _Float64x[n * n];
    _Float64x *w = new _Float64x[n];

//setting A matrix
    A[0 + 0 * n] = 5;    A[0 + 1 * n] = 4;    A[0 + 2 * n] = 1;    A[0 + 3 * n] = 1;
    A[1 + 0 * n] = 4;    A[1 + 1 * n] = 5;    A[1 + 2 * n] = 1;    A[1 + 3 * n] = 1;
    A[2 + 0 * n] = 1;    A[2 + 1 * n] = 1;    A[2 + 2 * n] = 4;    A[2 + 3 * n] = 2;
    A[3 + 0 * n] = 1;    A[3 + 1 * n] = 1;    A[3 + 2 * n] = 2;    A[3 + 3 * n] = 4;

    printf("A ="); printmat(n, n, A, n); printf("\n");
//work space query
    lwork = -1;
    _Float64x *work = new _Float64x[1];

    Rsyev("V", "U", n, A, n, w, work, lwork, info);
    lwork = (int) cast2double (work[0]);
    delete[]work;
    work = new _Float64x[std::max((mplapackint) 1, lwork)];
//inverse matrix
    Rsyev("V", "U", n, A, n, w, work, lwork, info);
//print out some results.
    printf("#eigenvalues \n");
    printf("w ="); printmat(n, 1, w, 1); printf("\n");

    printf("#eigenvecs \n");
    printf("U ="); printmat(n, n, A, n); printf("\n");
    printf("#you can check eigenvalues using octave/Matlab by:\n");
    printf("eig(A)\n");
    printf("#you can check eigenvectors using octave/Matlab by:\n");
    printf("U'*A*U\n");

    delete[]work;
    delete[]w;
    delete[]A;
}
