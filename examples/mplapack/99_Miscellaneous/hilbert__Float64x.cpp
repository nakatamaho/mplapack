// Get inverse of Hilbert matrix via Rgetri, using _Float128
// Very difficult to calculate the inverse due to very large
// condition number.
// For details, see http://en.wikipedia.org/wiki/Hilbert_matrix
// This file is freely usable.
// written by Nakata Maho, 2012/5/30.

#include <mpblas__Float64x.h>
#include <mplapack__Float64x.h>
#include <stdio.h>

#define FLOAT64X_FORMAT "%+25.21Le"
#define FLOAT64X_SHORT_FORMAT "%+20.16Le"

void printnum(_Float64x rtmp)
{
    printf(FLOAT64X_FORMAT, rtmp);
    return;
}
//Matlab/Octave format
void printmat(int N, int M, _Float64x *A, int LDA)
{
    _Float64x mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
            printnum(mtmp);
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
    mplapackint lwork, info;
    _Float64x *A = new _Float64x[n * n];	//A is overwritten
    _Float64x *Aorg = new _Float64x[n * n];
    _Float64x *C = new _Float64x[n * n];
    mplapackint *ipiv = new mplapackint[n];
    _Float64x One = 1.0, Zero = 0.0, mtmp;

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
    _Float64x *work = new _Float64x[1];

    Rgetri(n, A, n, ipiv, work, lwork, info);
    lwork = int (work[0]);
    delete[]work;
    work = new _Float64x[std::max(1, (int) lwork)];
//inverse matrix
    Rgetrf(n, n, A, n, ipiv, info);
    Rgetri(n, A, n, ipiv, work, lwork, info);

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
    printf("%.16e\n", (double)mtmp);
//quadmath libc support is needed to enable this, otherwise, use printnum above or
//cast to double.
//  printf("%.16Qe\n", mtmp);
    delete[]A;
    delete[]Aorg;
    delete[]C;
    delete[]ipiv;
    delete[]work;
}

int main()
{
    for (int n = 1; n < 21; n++) {
	printf("# inversion of Hilbert matrix of order n=%d\n", n);
	inv_hilbert_matrix(n);
    }
}
