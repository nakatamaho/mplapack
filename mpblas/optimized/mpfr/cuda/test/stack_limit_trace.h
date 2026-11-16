// Test-executable-only tracing of successful stack-limit API transitions.
#ifndef MPLAPACK_MPFR_STACK_LIMIT_TRACE_H
#define MPLAPACK_MPFR_STACK_LIMIT_TRACE_H
#include <cstddef>
namespace mpfr_stack_trace {
struct counts { size_t peak; unsigned raises,lowers; };
void begin(size_t initial);
counts end();
}
#endif
