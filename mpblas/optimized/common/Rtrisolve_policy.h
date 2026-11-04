#ifndef MPLAPACK_TRISOLVE_POLICY_H
#define MPLAPACK_TRISOLVE_POLICY_H
#include "Rtrisolve_jobs.h"
#include "Rblas1_policy.h"
#include <limits>
namespace mplapack_trisolve __attribute__((visibility("hidden"))) {
inline char flag(const char *s) { char c=s?*s:0;return c>='a' && c<='z'?c-32:c; }
// Validate all index and byte arithmetic before inspecting caller memory.
inline bool validate(shape &s,int storage,const char *u,const char *t,const char *d,
    long n,long k,long lda,long inc,size_t host_bytes,size_t device_bytes) {
    char U=flag(u),T=flag(t),D=flag(d);
    if((U!='U' && U!='L') || (T!='N' && T!='T' && T!='C') || (D!='N' && D!='U') ||
       !mplapack_blas1::stride_ok(n,inc)) return false;
    const long max=std::numeric_limits<long>::max();
    if(storage==1 && (k<0 || k==max || lda<k+1)) return false;
    if(storage==0 && lda<n) return false;
    if(n==max || (storage==2 && n>max/(n+1))) return false;
    if(storage!=2 && (lda<=0 || n>max/lda)) return false;
    s={n,k,lda,storage,U=='U',T=='N',T=='C',D=='U'};
    size_t count=s.count(),lim=std::numeric_limits<size_t>::max();
    size_t step=inc<0?size_t(-inc):size_t(inc);
    if(count>lim/host_bytes || (size_t(n-1)*step+1)>lim/host_bytes) return false;
    return count<lim/device_bytes && size_t(n)<=((lim/device_bytes)-count-1)/2;
}
}
#endif
