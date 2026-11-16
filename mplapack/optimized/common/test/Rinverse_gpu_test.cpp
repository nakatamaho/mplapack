// Direct BF inverse entries: no dispatcher fallback is allowed for n>0.
#include "../Rinverse_cuda_entry.cuh"
#include <string>
namespace {
int inverse_calls=0,inverse_rejections=0;
#define INVERSE_DIRECT_WRAPPERS(H,P) \
static void direct_##P##trtri(const char *u,const char *d,INTEGER n,H *a,INTEGER ld,INTEGER &info) { \
    if(!n) {P##trtri(u,d,n,a,ld,info);return;} ++inverse_calls;if(!TEST_ENTRY(P##trtri)(u,d,n,a,ld,info))++inverse_rejections; } \
static void direct_##P##lauum(const char *u,INTEGER n,H *a,INTEGER ld,INTEGER &info) { \
    if(!n) {P##lauum(u,n,a,ld,info);return;} ++inverse_calls;if(!TEST_ENTRY(P##lauum)(u,n,a,ld,info))++inverse_rejections; } \
static void direct_##P##potri(const char *u,INTEGER n,H *a,INTEGER ld,INTEGER &info) { \
    if(!n) {P##potri(u,n,a,ld,info);return;} ++inverse_calls;if(!TEST_ENTRY(P##potri)(u,n,a,ld,info))++inverse_rejections; } \
static void direct_##P##getri(INTEGER n,H *a,INTEGER ld,INTEGER *p,H *w,INTEGER lw,INTEGER &info) { \
    if(!n || lw==-1) {P##getri(n,a,ld,p,w,lw,info);return;} ++inverse_calls;if(!TEST_ENTRY(P##getri)(n,a,ld,p,w,lw,info))++inverse_rejections; }
#if defined(MPLAPACK_BUILD_WITH_DD)
#define TEST_ENTRY(R) R##_dd_cuda
#define TEST_PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define TEST_ENTRY(R) R##_td_cuda
#define TEST_PREFIX "MPLAPACK_TD"
#else
#define TEST_ENTRY(R) R##_qd_cuda
#define TEST_PREFIX "MPLAPACK_QD"
#endif
INVERSE_DIRECT_WRAPPERS(REAL,R)
INVERSE_DIRECT_WRAPPERS(COMPLEX,C)
#undef INVERSE_DIRECT_WRAPPERS

template<class T> struct inverse_contract_api;
#define INVERSE_CONTRACT_API(H,P) \
template<> struct inverse_contract_api<H> { \
static bool tri(const char *u,const char *d,INTEGER n,H *a,INTEGER ld,INTEGER &i){return TEST_ENTRY(P##trtri)(u,d,n,a,ld,i);} \
static bool lau(const char *u,INTEGER n,H *a,INTEGER ld,INTEGER &i){return TEST_ENTRY(P##lauum)(u,n,a,ld,i);} \
static bool pot(const char *u,INTEGER n,H *a,INTEGER ld,INTEGER &i){return TEST_ENTRY(P##potri)(u,n,a,ld,i);} \
static bool get(INTEGER n,H *a,INTEGER ld,INTEGER *p,H *w,INTEGER lw,INTEGER &i){return TEST_ENTRY(P##getri)(n,a,ld,p,w,lw,i);} \
};
INVERSE_CONTRACT_API(REAL,R)
INVERSE_CONTRACT_API(COMPLEX,C)
#undef INVERSE_CONTRACT_API
static void inverse_test_setting(const char *suffix,const char *value) {std::string key=std::string(TEST_PREFIX)+suffix;setenv(key.c_str(),value,1);}
template<class T> static int inverse_contracts() {
    using API=inverse_contract_api<T>;
    std::vector<T> a(28,T(3)),w(128,T(17));std::vector<INTEGER> p={1,2,3,4};
    const auto before=a,beforew=w;const auto beforep=p;INTEGER info=-77;int bad=0,cases=0;
    auto reject=[&](bool computed) {
        ++cases;
        if(computed || info!=-77 || std::memcmp(a.data(),before.data(),a.size()*sizeof(T)) ||
           std::memcmp(w.data(),beforew.data(),w.size()*sizeof(T)) || p!=beforep) {
            ++bad;std::fprintf(stderr,"FAIL inverse rejection/no-publication case %d\n",cases);
        }
        a=before;w=beforew;p=beforep;info=-77;
    };
    for(INTEGER n:{INTEGER(-1),INTEGER(0)}) {
        reject(API::tri("U","N",n,a.data(),7,info));reject(API::lau("L",n,a.data(),7,info));
        reject(API::pot("U",n,a.data(),7,info));reject(API::get(n,a.data(),7,p.data(),w.data(),128,info));
    }
    reject(API::tri("X","N",4,a.data(),7,info));reject(API::tri("U","X",4,a.data(),7,info));
    reject(API::lau("X",4,a.data(),7,info));reject(API::pot("X",4,a.data(),7,info));
    reject(API::tri("U","N",4,a.data(),3,info));reject(API::lau("L",4,a.data(),3,info));
    reject(API::pot("U",4,a.data(),3,info));reject(API::get(4,a.data(),3,p.data(),w.data(),128,info));
    reject(API::get(4,a.data(),7,p.data(),w.data(),-1,info));reject(API::get(4,a.data(),7,p.data(),w.data(),3,info));
    reject(API::get(4,a.data(),7,nullptr,w.data(),128,info));reject(API::get(4,a.data(),7,p.data(),nullptr,128,info));
    p[2]=0;const auto invalidp=p;
    // Restore p after checking its actual invalid input, not the valid baseline.
    bool result=API::get(4,a.data(),7,p.data(),w.data(),128,info);
    ++cases;if(result || p!=invalidp || info!=-77 || std::memcmp(a.data(),before.data(),a.size()*sizeof(T)) ||
        std::memcmp(w.data(),beforew.data(),w.size()*sizeof(T)))++bad;p=beforep;
    reject(API::tri("U","N",4,nullptr,7,info));
    const INTEGER huge=std::numeric_limits<INTEGER>::max();reject(API::tri("U","N",huge,a.data(),huge,info));
    for(const char *suffix:{"_LAPACK_CUDA","_LAPACK_INVERSE_CUDA_MIN"}) {
        inverse_test_setting(suffix,!std::strcmp(suffix,"_LAPACK_CUDA")?"0":"1000");
        reject(API::tri("U","N",4,a.data(),7,info));reject(API::lau("L",4,a.data(),7,info));
        reject(API::pot("U",4,a.data(),7,info));reject(API::get(4,a.data(),7,p.data(),w.data(),128,info));
        inverse_test_setting(suffix,!std::strcmp(suffix,"_LAPACK_CUDA")?"1":"0");
    }
    std::printf("inverse rejection contracts: %d cases, %d failures\n",cases,bad);return bad;
}
}
#define MPLAPACK_INVERSE_TEST_ACTUAL(R) direct_##R
#define main inverse_numeric_main
#include "Rfactor_test.cpp"
#undef main
#undef MPLAPACK_INVERSE_TEST_ACTUAL
int main(int argc,char **argv) {
#if defined(__CUDACC__)
    int count=0;if(cudaGetDeviceCount(&count)!=cudaSuccess || !count)return 77;
#endif
    inverse_test_setting("_LAPACK_CUDA_MIN","0");inverse_test_setting("_LAPACK_INVERSE_CUDA_MIN","0");
    inverse_test_setting("_LAPACK_CUDA","1");
    int bad=inverse_contracts<REAL>()+inverse_contracts<COMPLEX>();
    const bool quick=argc>1 && (!std::strcmp(argv[1],"--quick") || !std::strcmp(argv[1],"--inverse-quick"));
    char name[]="inverse-test",q[]="--inverse-quick",full[]="--inverses";char *args[]={name,quick?q:full};
    int status=inverse_numeric_main(2,args);
#if defined(__CUDACC__)
    if(cudaDeviceSynchronize()!=cudaSuccess)++bad;
#endif
    std::printf("inverse direct calls: %d, unexpected rejections: %d\n",inverse_calls,inverse_rejections);
    return status || bad || inverse_rejections?1:0;
}
