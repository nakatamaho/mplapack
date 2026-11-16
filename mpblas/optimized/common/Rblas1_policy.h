// Per-call opt-in and overflow-safe BLAS1 indexing (host only).
#ifndef MPLAPACK_BLAS1_POLICY_H
#define MPLAPACK_BLAS1_POLICY_H
#include <cstdlib>
#include <limits>
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {
inline bool fast_enabled(const char *name="MPLAPACK_MPBLAS_FAST") {
    const char *s=std::getenv(name);
    return s && s[0]=='1' && s[1]=='\0';
}
inline bool stride_ok(long n,long inc) {
    if(n<=0 || inc==0 || inc==std::numeric_limits<long>::min()) return false;
    const unsigned long step=inc<0?-inc:inc;
    return static_cast<unsigned long>(n-1)<=static_cast<unsigned long>(std::numeric_limits<long>::max())/step;
}
}
#endif
