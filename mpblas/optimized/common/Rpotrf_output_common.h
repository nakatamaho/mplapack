// Publish only the region a left-looking Lower potrf has reached on failure.
#ifndef MPLAPACK_POTRF_OUTPUT_COMMON_H
#define MPLAPACK_POTRF_OUTPUT_COMMON_H
namespace mplapack_lapack_gpu {
static inline bool potrf_output(bool upper,long n,long nb,long info,long i,long j) {
    if(upper || info<=0 || nb<=1 || nb>=n) return true;
    const long first=((info-1)/nb)*nb;
    const long end=first+nb<n?first+nb:n;
    return j<first || (j<end && i<end);
}
}
#endif
