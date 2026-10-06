// Verify executable-owned test overrides without running numerical suites.
#include <mpblas.h>
#include <mplapack.h>
#include <fem.hpp>
#include <mplapack_lin.h>
#include <mplapack_debug.h>

int main() {
    nout = 6;
    infot = 1;
    srnamt = "Rgetrf";
    ok = true;
    lerr = false;
    REAL a[1];
    INTEGER pivot[1], info = 0;
    Rgetrf(-1, 1, a, 1, pivot, info);
    if (info != -1 || !lerr || !ok)
        return 1;
#ifdef TEST_DMD
    // Test common state must not replace the production DMD tuning defaults.
    Mxlaenv(9, 7);
    if (iMlaenv(9, "Rgesvd", " ", 1, 1, -1, -1) != 25)
#else
    Mxlaenv(1, 7);
    if (iMlaenv(1, "Rgetrf", " ", 1, 1, -1, -1) != 7)
#endif
        return 2;
    return 0;
}
