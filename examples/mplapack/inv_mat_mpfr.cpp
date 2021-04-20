// Get inverse of Matrix A via Rgetri, using MPFR.
// This file is freely usable.
// written by Nakata Maho, 2009/9/23.

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

//Matlab/Octave format
void printmat(int N, int M, mpreal * A, int LDA)
{
    mpreal mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    mpfr_printf("%5.2Re", mpfr_ptr(mtmp));
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
//initialization of MPFR
    int default_prec = 256;
    mpfr_set_default_prec(default_prec);

    mpreal *A = new mpreal[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting A matrix
    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 4;    A[0 + 2 * n] = 6;
    A[1 + 0 * n] = 2;    A[1 + 1 * n] = 9;    A[1 + 2 * n] = 14;
    A[2 + 0 * n] = 3;    A[2 + 1 * n] = 14;   A[2 + 2 * n] = 23;

    printf("A =");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    mpreal *work = new mpreal[1];

    Rgetri(n, A, n, ipiv, work, lwork, info);
    lwork = (int) double (work[0]);
    delete[]work;
    work = new mpreal[std::max(1, (int) lwork)];
//inverse matrix
    Rgetrf(n, n, A, n, ipiv, info);
    Rgetri(n, A, n, ipiv, work, lwork, info);

    printf("invA =");
    printmat(n, n, A, n);
    printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]A;
}
