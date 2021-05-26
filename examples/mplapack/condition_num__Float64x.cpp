// Get condition number of Hibert matrix via Rgecon
// written by Nakata Maho, 2021/5.

#include <mpblas__Float64x.h>
#include <mplapack__Float64x.h>
#include <stdio.h>
#define BUFLEN 1024

//Matlab/Octave format
void printmat(int N, int M, _Float64x *A, int LDA)
{
    _Float64x mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
            printf("%5.2e", (double)mtmp);
//quadmath libc support is needed to enable this, otherwise, use printnum above or
//cast to double.
//          printf("%5.2Qe", mtmp);
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

void getAinv(mplapackint n, mplapackint lda, _Float64x *A)
{
    mplapackint info;
    mplapackint lwork;
    /* pivot vector allocation */
    mplapackint *ipiv = new mplapackint[n];
    lwork = -1;
    _Float64x *work = new _Float64x[1];
    /* query work space */
    Rgetri(n, A, lda, ipiv, work, lwork, info);
    lwork = (int) (work[0]);
    delete[]work;
    work = new _Float64x[std::max((mplapackint) 1, lwork)];
    /* do inversion */
    Rgetrf(n, n, A, lda, ipiv, info);
    Rgetri(n, A, lda, ipiv, work, lwork, info);
    delete[]ipiv;

    if (info == 0)
	return;
    if (info > 0)
	printf("matrix is singular\n");
    exit(1);
    if (info < 0)
	printf("%d th argument had an illegal value\n", (int) info);
    exit(1);
}

_Float64x get_estimated_condition_num(const char *norm, mplapackint n, mplapackint lda, _Float64x *A)
{
    _Float64x anorm, cond, rcond;
    _Float64x *work = new _Float64x[std::max((mplapackint) 1, n * 4)];
    mplapackint *iwork = new mplapackint[std::max((mplapackint) 1, n)];
    mplapackint info;

    /* First, calculate norm */
    anorm = Rlange(norm, n, n, A, lda, work);
    /* Second, do LU factorization */
    Rgetrf(n, n, A, lda, iwork, info);
    /* Third, calculate estimated condition number */
    Rgecon(norm, n, A, lda, anorm, rcond, work, iwork, info);

    cond = 1.0 / rcond;

    delete[]work;
    delete[]iwork;
    return cond;
}

_Float64x get_exact_condition_num(const char *norm, mplapackint n, mplapackint lda, _Float64x *A)
{
    _Float64x *Ainv = new _Float64x[n * n];
    _Float64x *work = new _Float64x[std::max((mplapackint) 1, n)];
    _Float64x anorm, ainvnorm, cond;
//save A matrix
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    Ainv[i + j * n] = A[i + j * n];
	}
    }
//get inversion of A matrix
    getAinv(n, lda, Ainv);
//get norm of A and Ainv
    anorm = Rlange(norm, n, n, A, lda, work);
    ainvnorm = Rlange(norm, n, n, Ainv, lda, work);
//condition number
    cond = anorm * ainvnorm;

    delete[]Ainv;
    delete[]work;
    return cond;
}

void condition_number_demo(mplapackint n)
{
    mplapackint lwork, info;
    _Float64x *A = new _Float64x[n * n];
    _Float64x mtmp;
    _Float64x cond_est, cond;

//setting Hilbert matrix in A.
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtmp = (i + 1) + (j + 1) - 1;
	    A[i + j * n] = 1.0 / mtmp;
	}
    }
    cond = get_exact_condition_num("1", n, n, A);
    cond_est = get_estimated_condition_num("1", n, n, A);

    printf("Hilbert matrix of order %d\n", (int) n);
    printf("Estimated condition number: %10.6e\n", (double)cond_est);
    printf("    Exact Condition number: %10.6e\n", (double)cond);
    delete[]A;
}

int main()
{
    for (int n = 2; n < 25; n++) {
	condition_number_demo(n);
    }
}
