// Rgemm demo using DD
// This file is freely usable.
// written by Nakata Maho, 2009/9/24.

#include <mblas_dd.h>
#include <stdio.h>

//Matlab/Octave format
void printmat(int N, int M, dd_real * A, int LDA)
{
    dd_real mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%5.2e", mtmp.x[0]);
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
    mpackint n = 3;

    dd_real *A = new dd_real[n * n];
    dd_real *B = new dd_real[n * n];
    dd_real *C = new dd_real[n * n];
    dd_real alpha, beta;

//setting A matrix
    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 8;    A[0 + 2 * n] = 3;
    A[1 + 0 * n] = 2.5;  A[1 + 1 * n] = 10;   A[1 + 2 * n] = 8;
    A[2 + 0 * n] = 9;    A[2 + 1 * n] = -5;   A[2 + 2 * n] = -1;

    B[0 + 0 * n] = 9;    B[0 + 1 * n] = 8;    B[0 + 2 * n] = 3;
    B[1 + 0 * n] = 3;    B[1 + 1 * n] = -11;  B[1 + 2 * n] = 4.8;
    B[2 + 0 * n] = -8;   B[2 + 1 * n] = 6;    B[2 + 2 * n] = 1;

    C[0 + 0 * n] = 3;    C[0 + 1 * n] = 3;    C[0 + 2 * n] = 1.2;
    C[1 + 0 * n] = 8;    C[1 + 1 * n] = 4;    C[1 + 2 * n] = 8;
    C[2 + 0 * n] = 6;    C[2 + 1 * n] = 1;    C[2 + 2 * n] = -2;

    printf("Rgemm demo...\n");

    printf("A =");
    printmat(n, n, A, n);
    printf("\n");

    printf("B =");
    printmat(n, n, B, n);
    printf("\n");

    printf("C =");
    printmat(n, n, C, n);
    printf("\n");
    alpha = 3.0;
    beta = -2.0;
    Rgemm("n", "n", n, n, n, alpha, A, n, B, n, beta, C, n);

    printf("alpha = %5.3e\n", alpha.x[0]);
    printf("beta  = %5.3e\n", beta.x[0]);
    printf("ans =");
    printmat(n, n, C, n);
    printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * A * B + beta * C =\n");
    delete[]C;
    delete[]B;
    delete[]A;
}
