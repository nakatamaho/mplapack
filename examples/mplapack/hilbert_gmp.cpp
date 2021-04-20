// Get inverse of Hilbert matrix via Rgetri, using GMP
// Very difficult to calculate the inverse due to very large
// condition number.
// For details, see http://en.wikipedia.org/wiki/Hilbert_matrix
// This file is freely usable.
// written by Nakata Maho, 2009/9/23.


#include <mpblas_gmp.h>
#include <mplapack_gmp.h>

//Matlab/Octave format
void printmat(int N, int M, mpf_class * A, int LDA)
{
    mpf_class mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    gmp_printf("%5.2Fe", mtmp.get_mpf_t());
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
    mpf_class *A = new mpf_class[n * n];	//A is overwritten
    mpf_class *Aorg = new mpf_class[n * n];
    mpf_class *C = new mpf_class[n * n];
    mplapackint *ipiv = new mplapackint[n];
    mpf_class One = 1.0, Zero = 0.0, mtmp;

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
    mpf_class *work = new mpf_class[1];

    Rgetri(n, A, n, ipiv, work, lwork, info);
    lwork = int (work[0].get_d());
    delete[]work;
    work = new mpf_class[std::max(1, (int) lwork)];
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
    printf("%.16e\n", mtmp.get_d());
    delete[]A;
    delete[]Aorg;
    delete[]C;
    delete[]ipiv;
    delete[]work;
}

int main()
{
//initialization of GMP
    int default_prec = 1024;
    mpf_set_default_prec(default_prec);

    for (int n = 1; n < 50; n++) {
	printf("# inversion of Hilbert matrix of order n=%d\n", n);
	inv_hilbert_matrix(n);
    }
}
