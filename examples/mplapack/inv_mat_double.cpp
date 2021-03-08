// Get inverse of Matrix A via Rgetri, using double
// This file is freely usable.
// written by Nakata Maho, 2009/9/23.

#include <mpblas_double.h>
#include <mplapack_double.h>
#include <stdio.h>

//Matlab/Octave format
void printmat(int N, int M, double *A, int LDA)
{
    double mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%5.2e", mtmp);
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

    double *A = new double[n * n];
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
    double *work = new double[1];

    Rgetri(n, A, n, ipiv, work, lwork, &info);
    lwork = int (work[0]);
    delete[]work;
    work = new double[std::max(1, (int) lwork)];
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
