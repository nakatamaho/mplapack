// Get inverse of Matrix A via Rgetri, using binary128
// This file is freely usable.
// written by Nakata Maho, 2012/5/30.

#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <stdio.h>

#define BUFLEN 1024

void printnum(binary128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp);
    if ((size_t) n < sizeof buf)
    printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printmat(int N, int M, binary128 *A, int LDA)
{
    binary128 mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%5.2e", (double)mtmp);
//quadmath libc support is needed to enable this, otherwise, use printnum above or
//cast to double.
//          printf("%5.2Qe", mtmp);
	    if (j < M - 1)
		printf(", ");
	}
	if (i < N - 1)
	    printf("]; ");
	else
	    printf("] ");
    }
    printf("]");
}

int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    binary128 *A = new binary128[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting A matrix
    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 4;    A[0 + 2 * n] = 6;
    A[1 + 0 * n] = 2;    A[1 + 1 * n] = 9;    A[1 + 2 * n] = 14;
    A[2 + 0 * n] = 3;    A[2 + 1 * n] = 14;   A[2 + 2 * n] = 23;

    printf("A = ");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    binary128 *work = new binary128[1];

    Rgetri(n, A, n, ipiv, work, lwork, &info);
    lwork = int (work[0]);
    delete[]work;
    work = new binary128[std::max(1, (int) lwork)];
//inverse matrix
    Rgetrf(n, n, A, n, ipiv, &info);
    Rgetri(n, A, n, ipiv, work, lwork, &info);

    printf("invA = ");
    printmat(n, n, A, n);
    printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]A;
}
