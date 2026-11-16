#include <mplapack_td.h>

int main() {
    td_real value(1.0);
    return Risnan(value) ? 1 : 0;
}
