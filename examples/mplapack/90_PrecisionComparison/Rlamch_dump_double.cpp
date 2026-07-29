//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }

// Matlab/Octave format
void printvec(double *a, int len) {
    double tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, double *a, int lda)
{
    double mtmp;

    printf("[ ");
    for (int i = 0; i < n; i++) {
        printf("[ ");
        for (int j = 0; j < m; j++) {
            mtmp = a[i + j * lda];
            printnum(mtmp);
            if (j < m - 1)
                printf(", ");
        }
        if (i < n - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}
int main() {
    printf("backend = double\n");
    const char *names = "ESBPNRMO";
    const char *labels[] = {"eps", "sfmin", "base", "precision", "t", "rmin", "rmax"};
    const char codes[] = {'E', 'S', 'B', 'P', 'N', 'U', 'O'};
    for (int i = 0; i < 7; i++) {
        printf("%s = ", labels[i]);
        char c[2];
        c[0] = codes[i];
        c[1] = '\0';
        printnum(Rlamch_double(c));
        printf("\n");
    }
    return 0;
}
