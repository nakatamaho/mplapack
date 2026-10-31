// Linux-only test interposition: model staged, unfinished pageable H2D copies.
// The default mode delays MPFR pivots; BF mode delays the {1,-1} constants.
#include <cuda_runtime_api.h>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <thread>
#include <vector>
namespace {
using Copy=cudaError_t(*)(void*,const void*,size_t,cudaMemcpyKind);
using Sync=cudaError_t(*)(cudaStream_t);
Copy real_copy() { static auto f=reinterpret_cast<Copy>(dlsym(RTLD_NEXT,"cudaMemcpy")); return f; }
Sync real_sync() { static auto f=reinterpret_cast<Sync>(dlsym(RTLD_NEXT,"cudaStreamSynchronize")); return f; }
bool pending=false;
int injections=0,waits=0,skips=0;
std::vector<void*> staging;
void CUDART_CB delay(void*) {
    // CUDA APIs are forbidden in a host callback; only delay the default stream.
    std::this_thread::sleep_for(std::chrono::seconds(1));
}
}
extern "C" cudaError_t cudaMemcpy(void *dst,const void *src,size_t bytes,cudaMemcpyKind kind) {
    constexpr long n=65;
    if(kind!=cudaMemcpyHostToDevice) return real_copy()(dst,src,bytes,kind);
    const char *mode=std::getenv("MPLAPACK_TEST_DELAY_BF_CONSTANTS");
    bool bf=mode && *mode=='1';
    std::vector<unsigned char> stale;
    if(bf) {
        if(bytes!=32 && bytes!=48 && bytes!=64) return real_copy()(dst,src,bytes,kind);
        const auto *v=static_cast<const double*>(src); size_t words=bytes/sizeof(double);
        for(size_t i=0;i<words;++i) if(v[i]!=(i==0?1.0:i==words/2?-1.0:0.0)) return real_copy()(dst,src,bytes,kind);
        stale.assign(bytes,0);
    } else {
        if(bytes!=n*sizeof(long)) return real_copy()(dst,src,bytes,kind);
        const auto *p=static_cast<const long*>(src);
        if(p[0]!=n) return real_copy()(dst,src,bytes,kind);
        for(long i=1;i<n;++i) if(p[i]!=i+1) return real_copy()(dst,src,bytes,kind);
        stale.resize(bytes);
        // Valid but stale identity pivots avoid an illegal-memory control.
        for(long i=0;i<n;++i) { long value=i+1; std::memcpy(stale.data()+i*sizeof(long),&value,sizeof(value)); }
    }
    auto rc=real_copy()(dst,stale.data(),bytes,kind);
    if(rc==cudaSuccess) rc=real_sync()(nullptr);
    void *copy=nullptr;
    if(rc==cudaSuccess) rc=cudaMallocHost(&copy,bytes);
    if(rc!=cudaSuccess) return rc;
    staging.push_back(copy); std::memcpy(copy,src,bytes);
    rc=cudaLaunchHostFunc(nullptr,delay,nullptr);
    if(rc==cudaSuccess) rc=cudaMemcpyAsync(dst,copy,bytes,kind,nullptr);
    if(rc==cudaSuccess) { pending=true; ++injections; }
    return rc;
}
extern "C" cudaError_t cudaStreamSynchronize(cudaStream_t stream) {
    if(!stream && pending) {
        pending=false;
        if(const char *s=std::getenv("MPLAPACK_TEST_SKIP_COPY_FENCE")) if(*s=='1') { ++skips; return cudaSuccess; }
        ++waits;
    }
    return real_sync()(stream);
}
extern "C" int mplapack_test_pivot_injections() { return injections; }
extern "C" int mplapack_test_pivot_waits() { return waits; }
extern "C" int mplapack_test_pivot_skips() { return skips; }
extern "C" void mplapack_test_pivot_release() {
    real_sync()(nullptr);
    for(void *p:staging) cudaFreeHost(p);
    staging.clear();
}
