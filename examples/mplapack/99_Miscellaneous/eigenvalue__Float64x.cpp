// Get eigenvalues and eigenvecs
// of Matrix A via Rsyev, using _Float128.
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

int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    _Float64x *A = new _Float64x[n * n];
    _Float64x *w = new _Float64x[n];

//setting A matrix
    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 2;    A[0 + 2 * n] = 3;
    A[1 + 0 * n] = 2;    A[1 + 1 * n] = 5;    A[1 + 2 * n] = 4;
    A[2 + 0 * n] = 3;    A[2 + 1 * n] = 4;    A[2 + 2 * n] = 6;

    printf("A =");
    printmat(n, n, A, n);
    printf("\n");
//work space query
    lwork = -1;
    _Float64x *work = new _Float64x[1];

    Rsyev("V", "U", n, A, n, w, work, lwork, info);
    lwork = (int) work[0];
    delete[]work;
    work = new _Float64x[std::max((mplapackint) 1, lwork)];
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
