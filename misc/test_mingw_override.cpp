// Exercise the actual test error handler without a numerical workload.
#include <mpblas.h>
#include <fem.hpp>
#include <mplapack_debug.h>

int main() {
    nout = 6;
    infot = 1;
    srnamt = "Rgetrf";
    ok = true;
    lerr = false;
    Mxerbla("Rgetrf", 1);
    return (!lerr || !ok) ? 2 : 0;
}
