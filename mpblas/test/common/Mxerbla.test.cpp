#include <mpblas.h>
#include <blas.h>
#include <mplapack_debug.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

void Mxerbla_test() {
    Mxerbla("Fasum", 10);
    Mxerbla("Maho", 100);
}

int main(int argc, char *argv[]) {
    printf("*** Testing Mxerbla start ***\n");
    Mxerbla_test();
    printf("*** Testing Mxerbla successful ***\n");
    return (0);
}
