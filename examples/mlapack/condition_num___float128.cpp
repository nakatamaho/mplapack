// Get condition number of Hibert matrix via Rgecon
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

void getAinv(mpackint n, mpackint lda, __float128 *A)
{
    mpackint info;
    mpackint lwork;
    /* pivot vector allocation */
    mpackint *ipiv = new mpackint[n];
    lwork = -1;
    __float128 *work = new __float128[1];
    /* query work space */
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
    lwork = (int) (work[0]);
    delete[]work;
    work = new __float128[std::max((mpackint) 1, lwork)];
    /* do inversion */
    Rgetrf(n, n, A, lda, ipiv, &info);
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
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

__float128 get_estimated_condition_num(const char *norm, mpackint n, mpackint lda, __float128 *A)
{
    __float128 anorm, cond, rcond;
    __float128 *work = new __float128[std::max((mpackint) 1, n * 4)];
    mpackint *iwork = new mpackint[std::max((mpackint) 1, n)];
    mpackint info;

    /* First, calculate norm */
    anorm = Rlange(norm, n, n, A, lda, work);
    /* Second, do LU factorization */
    Rgetrf(n, n, A, lda, iwork, &info);
    /* Third, calculate estimated condition number */
    Rgecon(norm, n, A, lda, anorm, &rcond, work, iwork, &info);

    cond = 1.0 / rcond;

    delete[]work;
    delete[]iwork;
    return cond;
}

__float128 get_exact_condition_num(const char *norm, mpackint n, mpackint lda, __float128 *A)
{
    __float128 *Ainv = new __float128[n * n];
    __float128 *work = new __float128[std::max((mpackint) 1, n)];
    __float128 anorm, ainvnorm, cond;
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

void condition_number_demo(mpackint n)
{
    mpackint lwork, info;
    __float128 *A = new __float128[n * n];
    __float128 mtmp;
    __float128 cond_est, cond;

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
