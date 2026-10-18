// Verify executable-owned test overrides without running numerical suites.
#include <mpblas.h>
#include <mplapack.h>
#include <fem.hpp>
#include <mplapack_lin.h>
#include <mplapack_debug.h>
#include <cstdlib>

int main() {
#if defined(___MPLAPACK_BUILD_WITH_GMP___)
    // Archive extraction must not silently leave GMP at its 64-bit default.
    // Honor the same environment override as the backend constructor.
    const char *precision = std::getenv("MPLAPACK_GMP_PRECISION");
    const unsigned long requested = precision ? std::strtoul(precision, nullptr, 10)
                                              : ___MPLAPACK_GMP_DEFAULT_PRECISION___;
    const unsigned long limb = GMP_NUMB_BITS;
    const unsigned long expected = ((requested + limb - 1) / limb) * limb;
    if (mpf_get_default_prec() != expected)
        return 3;
#endif
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
