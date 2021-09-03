int main()
{
    INTEGER n = 4;
    INTEGER lwork, info;

    REAL *a = new REAL[n * n];
    INTEGER *ipiv = new INTEGER[n];

//setting a matrix
//https://www.tuhh.de/ti3/paper/rump/NiRuOi11.pdf
//Nonlinear Theory and Its Applications, IEICE, vol. 2, no. 2, pp. 226-245
//DOI: 10.1588/nolta.2.226
    a[0 + 0 * n] = 17;   a[0 + 1 * n] = -864; a[0 + 2 * n] = 716;    a[0 + 3 * n] = -799;
    a[1 + 0 * n] = 1;    a[1 + 1 * n] = -50;  a[1 + 2 * n] = 0.0;    a[1 + 3 * n] = 0.0;
    a[2 + 0 * n] = 0.0;  a[2 + 1 * n] = 1;    a[2 + 2 * n] = -50;    a[2 + 3 * n] = 0.0;
    a[3 + 0 * n] = 0.0;  a[3 + 1 * n] = 0.0;  a[3 + 2 * n] = 1;      a[3 + 3 * n] = -50;

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    REAL *work = new REAL[1];

    Rgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castInTEGER (work[0]);
    delete[]work;
    work = new REAL[std::max(1, (int) lwork)];

//inverse matrix
    Rgetrf(n, n, a, n, ipiv, info);
    Rgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
