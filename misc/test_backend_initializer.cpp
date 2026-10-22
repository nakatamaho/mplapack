// Start with extended x87 precision before the backend constructor runs.
// Linking the initializer must select 53-bit precision for libQD3 arithmetic.
#include <cstdio>

#if defined(__i386__) || defined(__x86_64__)
static void __attribute__((constructor(101))) preset_x87_precision() {
    unsigned short cw;
    asm volatile("fnstcw %0" : "=m"(cw));
    cw = (cw & ~0x0300u) | 0x0300u;
    asm volatile("fldcw %0" : : "m"(cw));
}
#endif

int main() {
#if defined(__i386__) || defined(__x86_64__)
    unsigned short cw;
    asm volatile("fnstcw %0" : "=m"(cw));
    std::printf("x87 precision control: 0x%03x\n", cw & 0x0300u);
    return (cw & 0x0300u) == 0x0200u ? 0 : 1;
#else
    return 77;
#endif
}
