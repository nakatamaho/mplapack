#include <stdio.h>
#include <string.h>

//  please use as following: Do not forget adding //CHAR(0) at the end.
//         CALL PRINTSTR('#dgelsd.f l.486'//CHAR(0))
//         CALL PRINTVEC('s='//CHAR(0),s,n)
//         CALL PRINTMAT('e='//CHAR(0),m,m,e,lde)
//         for complex matrix
//         CALL PRINTMATC('e='//CHAR(0),m,m,e,lde)
//         CALL PRINTMATC('x='//CHAR(0),x,n)

void printdummy_() { printf("\n"); }

void printstr_(char *s) {
    printf("%s\n", s);
}

void printnumi_(char *s, int *A) {
    printf("%s%d\n", s, *A);
}

void printnum_(char *s, double *A) {
    printf("%s%+21.16e\n", s, *A);
}

void printnumc(double _Complex rtmp) {
    printf("%+21.16e%+21.16ei", __real__ rtmp, __imag__ rtmp);
}

void printmatc_(char *s, int *N, int *M, double _Complex *A, int *LDA) {
    double _Complex tmp;
    printf("%s[ ", s);
    for (int i = 0; i < *N; i++) {
        printf("[ ");
        for (int j = 0; j < *M; j++) {
            tmp = A[i + j * (*LDA)];
            printnumc(tmp);
            if (j < *M - 1)
                printf(", ");
        }
        if (i < *N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
    printf("\n");
}

void printmat_(char *s, int *N, int *M, double *A, int *LDA) {
    double tmp;
    printf("%s[ ", s);
    for (int i = 0; i < *N; i++) {
        printf("[ ");
        for (int j = 0; j < *M; j++) {
            tmp = A[i + j * (*LDA)];
            printf("%+21.16e", tmp);
            if (j < *M - 1)
                printf(", ");
        }
        if (i < *N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
    printf("\n");
}

void printmati_(char *s, int *N, int *M, int *A, int *LDA) {
    int tmp;
    printf("%s[ ", s);
    for (int i = 0; i < *N; i++) {
        printf("[ ");
        for (int j = 0; j < *M; j++) {
            tmp = A[i + j * (*LDA)];
            printf("%d", tmp);
            if (j < *M - 1)
                printf(", ");
        }
        if (i < *N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
    printf("\n");
}

void printvec_(char *s, double *A, int *lenvec) {
    double tmp;
    printf("%s[ ", s);
    for (int i = 0; i < *lenvec; i++) {
        tmp = A[i];
        printf("%+21.16e", tmp);
        if (i < *lenvec - 1)
            printf(", ");
    }
    printf("]");
    printf("\n");
}

void printvecc_(char *s, double _Complex*A, int *lenvec) {
    double _Complex tmp;
    printf("%s[ ", s);
    for (int i = 0; i < *lenvec; i++) {
        printnumc(tmp);
        if (i < *lenvec - 1)
            printf(", ");
    }
    printf("]");
    printf("\n");
}

void printveci_(char *s, int *A, int *lenvec) {
    int tmp;
    printf("%s[ ", s);
    for (int i = 0; i < *lenvec; i++) {
        tmp = A[i];
        printf("%d", tmp);
        if (i < *lenvec - 1)
            printf(", ");
    }
    printf("]");
    printf("\n");
}
