// Get eigenvalues and eigenvecs
// of Matrix A via Rsyev, using GMP.
// This file is freely usable.
// written by Nakata Maho, 2009/9/24.

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
	    gmp_printf("%8.6Fe", mtmp.get_mpf_t());
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
//initialization of GMP
    int default_prec = 256;
    mpf_set_default_prec(default_prec);

    mpf_class *A = new mpf_class[n * n];
    mpf_class *w = new mpf_class[n];

//setting A matrix
    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 2;    A[0 + 2 * n] = 3;
    A[1 + 0 * n] = 2;    A[1 + 1 * n] = 5;    A[1 + 2 * n] = 4;
    A[2 + 0 * n] = 3;    A[2 + 1 * n] = 4;    A[2 + 2 * n] = 6;

    printf("A =");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    mpf_class *work = new mpf_class[1];

    Rsyev("V", "U", n, A, n, w, work, lwork, info);
    lwork = (int) work[0].get_d();
    delete[]work;
    work = new mpf_class[std::max((mplapackint) 1, lwork)];
//inverse matrix
    Rsyev("V", "U", n, A, n, w, work, lwork, info);
//print out some results.
    printf("#eigenvalues \n");
    printf("w =");
    printmat(n, 1, w, 1);
    printf("\n");
    printf("#eigenvecs \n");
    printf("U =");
    printmat(n, n, A, n);
    printf("\n");
    printf("#you can check eigenvalues using octave/Matlab by:\n");
    printf("eig(A)\n");
    printf("#you can check eigenvectors using octave/Matlab by:\n");
    printf("U'*A*U\n");

    delete[]work;
    delete[]w;
    delete[]A;
}
