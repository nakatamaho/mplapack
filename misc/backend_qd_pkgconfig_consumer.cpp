#include <mplapack_qd.h>

int main() {
    qd_real value(1.0);
    return Risnan(value) ? 1 : 0;
}
