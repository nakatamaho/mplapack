// Get inverse of Hilbert matrix via Rgetri, using quad-double
// Very difficult to calculate the inverse due to very large
// condition number.
// For details, see http://en.wikipedia.org/wiki/Hilbert_matrix
// This file is freely usable.
// written by Nakata Maho, 2009/9/23.

#include <mblas_qd.h>
#include <mlapack_qd.h>
#include <stdio.h>

//Matlab/Octave format
void printmat(int N, int M, qd_real * A, int LDA)
{
    qd_real mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%3.1f", mtmp.x[0]);
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

void inv_hilbert_matrix(int n)
{
    mpackint lwork, info;
    qd_real *A = new qd_real[n * n];	//A is overwritten
    qd_real *Aorg = new qd_real[n * n];
    qd_real *C = new qd_real[n * n];
    mpackint *ipiv = new mpackint[n];
    qd_real One = 1.0, Zero = 0.0, mtmp;

//setting A matrix
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtmp = (i + 1) + (j + 1) - 1;
	    A[i + j * n] = One / mtmp;
	    Aorg[i + j * n] = One / mtmp;
	}
    }

    printf("A = ");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    qd_real *work = new qd_real[1];

    Rgetri(n, A, n, ipiv, work, lwork, &info);
    lwork = int (work[0].x[0]);
    delete[]work;
    work = new qd_real[std::max(1, (int) lwork)];
//inverse matrix
    Rgetrf(n, n, A, n, ipiv, &info);
    Rgetri(n, A, n, ipiv, work, lwork, &info);

    printf("invA = ");
    printmat(n, n, A, n);
    printf("\n");

    One = 1.0, Zero = 0.0;
    Rgemm("N", "N", n, n, n, One, Aorg, n, A, n, Zero, C, n);
    printf("A * invA =");
    printmat(n, n, C, n);
    printf("\n");

    printf("Infnorm(A*invA)=");
    mtmp = 0.0;
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    if (mtmp < abs(C[i + j * n] - ((i == j) ? 1.0 : 0.0)))
		mtmp = abs(C[i + j * n] - ((i == j) ? 1.0 : 0.0));
	}
    }
    printf("%.16e\n", mtmp.x[0]);
    delete[]A;
    delete[]Aorg;
    delete[]C;
    delete[]ipiv;
    delete[]work;
}

int main()
{
    for (int n = 1; n < 40; n++) {
	printf("# inversion of Hilbert matrix of order n=%d\n", n);
	inv_hilbert_matrix(n);
    }
}
