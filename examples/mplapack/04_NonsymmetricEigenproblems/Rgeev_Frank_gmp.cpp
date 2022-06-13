//public domain
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>

#define GMP_FORMAT "%+68.64Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"

inline void printnum(mpf_class rtmp) { gmp_printf(GMP_FORMAT, rtmp.get_mpf_t()); }
inline void printnum_short(mpf_class rtmp) { gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t()); }

//Matlab/Octave format
void printvec(mpf_class *a, int len) {
    mpf_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpf_class * a, int lda)
{
    mpf_class mtmp;

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
void Frank(mplapackint n) {
    mplapackint lwork, liwork, info, m;
    mpf_class *a = new mpf_class[n * n];
    mpf_class *vl = new mpf_class[n * n]; //not used
    mpf_class *vr = new mpf_class[n * n]; //not used
    mpf_class *wr = new mpf_class[n];
    mpf_class *wi = new mpf_class[n];

    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = 0.0;
	}
    }
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = n - std::max(i, j) + 1;
        }
    }
    for (int i = 1; i <= n - 1 ; i++) {
        a[i + (i - 1) * n] = n - i;
    }
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");

    // work space query
    lwork = -1;
    mpf_class *work = new mpf_class[1];
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new mpf_class[std::max((mplapackint)1, lwork)];

    // diagonalize matrix
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);

    // print out
    printf("#eigenvalues \n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }

    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vr;
    delete[] vl;
    delete[] a;
}

bool rselect(mpf_class ar, mpf_class ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main(int argc, char *argv[]) {
    int STARTN = 5;
    int ENDN = 25;
    int STEPN = 1;
    if (argc != 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STARTN", argv[i]) == 0) {
                STARTN = atoi(argv[++i]);
            } else if (strcmp("-ENDN", argv[i]) == 0) {
                ENDN = atoi(argv[++i]);
            }
        }
    }
    for (int n = STARTN; n <= ENDN; n = n + STEPN) {
        printf("# Eigenvalues of Frank matrix of order n=%d\n", n);
        Frank((mplapackint)n);
    }
}
