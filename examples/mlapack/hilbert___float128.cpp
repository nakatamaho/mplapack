// Get inverse of Hilbert matrix via Rgetri, using __float128
// Very difficult to calculate the inverse due to very large
// condition number.
// For details, see http://en.wikipedia.org/wiki/Hilbert_matrix
// This file is freely usable.
// written by Nakata Maho, 2012/5/30.

#include <mblas___float128.h>
#include <mlapack___float128.h>
#include <stdio.h>
#define BUFLEN 1024

void printnum(__float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp);
    if ((size_t) n < sizeof buf)
    printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printmat(int N, int M, __float128 *A, int LDA)
{
    __float128 mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
            printf("%3.1f", (double)mtmp);
//quadmath libc support is needed to enable this, otherwise, use printnum above or
//cast to double.
//          printf("%3.1Qe", mtmp);
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
    __float128 *A = new __float128[n * n];	//A is overwritten
    __float128 *Aorg = new __float128[n * n];
    __float128 *C = new __float128[n * n];
    mpackint *ipiv = new mpackint[n];
    __float128 One = 1.0, Zero = 0.0, mtmp;

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
    __float128 *work = new __float128[1];

    Rgetri(n, A, n, ipiv, work, lwork, &info);
    lwork = int (work[0]);
    delete[]work;
    work = new __float128[std::max(1, (int) lwork)];
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
	    if (mtmp < fabs(C[i + j * n] - ((i == j) ? 1.0 : 0.0)))
		mtmp = fabs(C[i + j * n] - ((i == j) ? 1.0 : 0.0));
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
