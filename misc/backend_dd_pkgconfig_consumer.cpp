#include <mplapack_dd.h>

int main() {
    dd_real value(1.0);
    return Risnan(value) ? 1 : 0;
}
