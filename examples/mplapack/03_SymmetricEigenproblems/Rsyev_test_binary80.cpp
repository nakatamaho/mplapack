//public domain
#include <mpblas_binary80.h>
#include <mplapack_binary80.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#define BUFLEN 1024
void printnum(mplapack_binary80_t rtmp)
{
    int width = 25;
    char buf[BUFLEN];
#if MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_STRFROMF64X
    strfromf64x(buf, sizeof(buf), "%.21e", rtmp);
#elif MPLAPACK_BINARY80_IO == MPLAPACK_BINARY80_IO_SNPRINTF_LDBL
    snprintf(buf, sizeof(buf), "%*.21Le", width, rtmp);
#else
    #error "unsupported binary80 type"
#endif
    if (rtmp >= 0.0)
        printf("+%s", buf);
    else
        printf("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(mplapack_binary80_t *a, int len) {
    mplapack_binary80_t tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mplapack_binary80_t *a, int lda)
{
    mplapack_binary80_t mtmp;

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

    mplapack_binary80_t *A = new mplapack_binary80_t[n * n];
    mplapack_binary80_t *w = new mplapack_binary80_t[n];

//setting A matrix
    A[0 + 0 * n] = 5.0;    A[0 + 1 * n] = 4.0;    A[0 + 2 * n] = 1.0;    A[0 + 3 * n] = 1.0;
    A[1 + 0 * n] = 4.0;    A[1 + 1 * n] = 5.0;    A[1 + 2 * n] = 1.0;    A[1 + 3 * n] = 1.0;
    A[2 + 0 * n] = 1.0;    A[2 + 1 * n] = 1.0;    A[2 + 2 * n] = 4.0;    A[2 + 3 * n] = 2.0;
    A[3 + 0 * n] = 1.0;    A[3 + 1 * n] = 1.0;    A[3 + 2 * n] = 2.0;    A[3 + 3 * n] = 4.0;

    printf("A ="); printmat(n, n, A, n); printf("\n");
//work space query
    lwork = -1;
    mplapack_binary80_t *work = new mplapack_binary80_t[1];

    Rsyev("V", "U", n, A, n, w, work, lwork, info);
    lwork = (int) cast2double (work[0]);
    delete[]work;
    work = new mplapack_binary80_t[std::max((mplapackint) 1, lwork)];
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
