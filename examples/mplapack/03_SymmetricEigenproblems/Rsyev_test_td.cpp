//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_td.h>
#include <mplapack_td.h>

#define TD_PRECISION_SHORT 16

inline void printnum(td_real rtmp) {
    std::cout.precision(TD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(td_real *a, int len) {
    td_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, td_real * a, int lda)
{
    td_real mtmp;
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

    td_real *A = new td_real[n * n];
    td_real *w = new td_real[n];

//setting A matrix
    A[0 + 0 * n] = 5.0;    A[0 + 1 * n] = 4.0;    A[0 + 2 * n] = 1.0;    A[0 + 3 * n] = 1.0;
    A[1 + 0 * n] = 4.0;    A[1 + 1 * n] = 5.0;    A[1 + 2 * n] = 1.0;    A[1 + 3 * n] = 1.0;
    A[2 + 0 * n] = 1.0;    A[2 + 1 * n] = 1.0;    A[2 + 2 * n] = 4.0;    A[2 + 3 * n] = 2.0;
    A[3 + 0 * n] = 1.0;    A[3 + 1 * n] = 1.0;    A[3 + 2 * n] = 2.0;    A[3 + 3 * n] = 4.0;

    printf("A ="); printmat(n, n, A, n); printf("\n");
//work space query
    lwork = -1;
    td_real *work = new td_real[1];

    Rsyev("V", "U", n, A, n, w, work, lwork, info);
    lwork = (int) cast2double (work[0]);
    delete[]work;
    work = new td_real[std::max((mplapackint) 1, lwork)];
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
