int main()
{
    INTEGER n = 3;

    REAL *a = new REAL[n * n];
    REAL *b = new REAL[n * n];
    REAL *c = new REAL[n * n];
    REAL alpha, beta;

//setting A matrix
    a[0 + 0 * n] = 1;    a[0 + 1 * n] = 8;    a[0 + 2 * n] = 3;
    a[1 + 0 * n] = 2;    a[1 + 1 * n] = 10;   a[1 + 2 * n] = 8;
    a[2 + 0 * n] = 9;    a[2 + 1 * n] = -5;   a[2 + 2 * n] = -1;

    b[0 + 0 * n] = 9;    b[0 + 1 * n] = 8;    b[0 + 2 * n] = 3;
    b[1 + 0 * n] = 3;    b[1 + 1 * n] = -11;  b[1 + 2 * n] = 8;
    b[2 + 0 * n] = -8;   b[2 + 1 * n] = 6;    b[2 + 2 * n] = 1;

    c[0 + 0 * n] = 3;    c[0 + 1 * n] = 3;    c[0 + 2 * n] = -9;
    c[1 + 0 * n] = 8;    c[1 + 1 * n] = 4;    c[1 + 2 * n] = 8;
    c[2 + 0 * n] = 6;    c[2 + 1 * n] = 1;    c[2 + 2 * n] = -2;

    printf("# Rgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = 3.0;
    beta = -2.0;
    Rgemm("n", "n", n, n, n, alpha, a, n, b, n, beta, c, n);

    printf("alpha = "); printnum(alpha); printf("\n");
    printf("beta = "); printnum(beta); printf("\n");
    printf("ans ="); printmat(n, n, c, n); printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * a * b + beta * c \n");
    delete[]c;
    delete[]b;
    delete[]a;
}
