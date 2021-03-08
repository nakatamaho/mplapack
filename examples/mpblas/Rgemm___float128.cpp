// Rgemm demo using __float128
// This file is freely usable.
// written by Nakata Maho, 2012/5/30.

#include <mpblas___float128.h>
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
void printmat(int N, int M, __float128 * A, int LDA)
{
    __float128 mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%5.2e", (double)mtmp);
//quadmath libc support is needed to enable this, otherwise, use printnum above or
//cast to double.
//	    printf("%5.2Qe", mtmp);
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

    __float128 *A = new __float128[n * n];
    __float128 *B = new __float128[n * n];
    __float128 *C = new __float128[n * n];
    __float128 alpha, beta;

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

//quadmath libc support is needed to enable this, otherwise, use printnum above or
//cast to double.
//    printf("alpha = %5.3Qe\n", alpha);
//    printf("beta  = %5.3Qe\n", beta);
    printf("alpha = %5.3e\n", (double)alpha);
    printf("beta  = %5.3e\n", (double)beta);
    printf("ans =");
    printmat(n, n, C, n);
    printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * A * B + beta * C =\n");
    delete[]C;
    delete[]B;
    delete[]A;
}
