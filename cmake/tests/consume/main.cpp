#include <mpblas_dd.h>
#include <mplapack_dd.h>

int main() {
    dd_real a[1] = {dd_real(2.0)};
    dd_real b[1] = {dd_real(3.0)};
    dd_real c[1] = {dd_real(0.0)};

    Rgemm("N", "N", 1, 1, 1, dd_real(1.0), a, 1, b, 1,
          dd_real(0.0), c, 1);

    mplapackint info = 0;
    Rpotrf("U", 1, c, 1, info);
    return info == 0 ? 0 : 1;
}
