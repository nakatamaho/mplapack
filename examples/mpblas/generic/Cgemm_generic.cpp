int main()
{
    INTEGER n = 3;

    COMPLEX *a = new COMPLEX[n * n];
    COMPLEX *b = new COMPLEX[n * n];
    COMPLEX *c = new COMPLEX[n * n];
    COMPLEX alpha, beta;

//setting A matrix
    a[0 + 0 * n] = COMPLEX(1.0,-1.0);    a[0 + 1 * n] = COMPLEX(8.0, 2.2);    a[0 + 2 * n] = COMPLEX(0.0, -10.0);
    a[1 + 0 * n] = COMPLEX(2.0, 0.0);    a[1 + 1 * n] = COMPLEX(10.0,0.0);    a[1 + 2 * n] = COMPLEX(8.1, 2.2);
    a[2 + 0 * n] = COMPLEX(-9.0,3.0);    a[2 + 1 * n] = COMPLEX(-5.0,3.0);    a[2 + 2 * n] = COMPLEX(-1.0, 0.0);

    b[0 + 0 * n] = COMPLEX(9.0, 0.0);    b[0 + 1 * n] = COMPLEX(8.0, -0.01);  b[0 + 2 * n] = COMPLEX(3.0, 1.001);
    b[1 + 0 * n] = COMPLEX(3.0, -8.0);   b[1 + 1 * n] = COMPLEX(-11.0, 0.1);  b[1 + 2 * n] = COMPLEX(8.0, 0.00001);
    b[2 + 0 * n] = COMPLEX(-8.0, 1.0);   b[2 + 1 * n] = COMPLEX(6.0, 0.0);    b[2 + 2 * n] = COMPLEX(1.1, 1.0);

    c[0 + 0 * n] = COMPLEX(3.0, 1.0);   c[0 + 1 * n] = COMPLEX(-3.0, 9.99);   c[0 + 2 * n] = COMPLEX(-9.0, -11.0);
    c[1 + 0 * n] = COMPLEX(8.0, -1.0);  c[1 + 1 * n] = COMPLEX(4.0, 4.44);    c[1 + 2 * n] = COMPLEX(8.0, 9.0);
    c[2 + 0 * n] = COMPLEX(6.0, 0.0);   c[2 + 1 * n] = COMPLEX(-1.0, 0.0);    c[2 + 2 * n] = COMPLEX(-2.0, 1.0);

    printf("# Cgemm demo...\n");

    printf("a ="); printmat(n, n, a, n); printf("\n");
    printf("b ="); printmat(n, n, b, n); printf("\n");
    printf("c ="); printmat(n, n, c, n); printf("\n");
    alpha = COMPLEX(3.0,-1.2);
    beta = COMPLEX(-2.0, -2.0);
    Cgemm("n", "n", n, n, n, alpha, a, n, b, n, beta, c, n);

    printf("alpha = "); printnum(alpha); printf("\n");
    printf("beta = "); printnum(beta); printf("\n");
    printf("ans ="); printmat(n, n, c, n); printf("\n");
    printf("#please check by Matlab or Octave following and ans above\n");
    printf("alpha * a * b + beta * c \n");
    delete[]c;
    delete[]b;
    delete[]a;
}
