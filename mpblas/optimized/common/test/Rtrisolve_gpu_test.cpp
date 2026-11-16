// All layouts/flags/increments, strict bit identity including signed zeros.
#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <type_traits>
#define Rtrsv ref_Rtrsv
#include "../../../reference/Rtrsv.cpp"
#undef Rtrsv
#define Rtbsv ref_Rtbsv
#include "../../../reference/Rtbsv.cpp"
#undef Rtbsv
#define Rtpsv ref_Rtpsv
#include "../../../reference/Rtpsv.cpp"
#undef Rtpsv
#define Ctrsv ref_Ctrsv
#include "../../../reference/Ctrsv.cpp"
#undef Ctrsv
#define Ctbsv ref_Ctbsv
#include "../../../reference/Ctbsv.cpp"
#undef Ctbsv
#define Ctpsv ref_Ctpsv
#include "../../../reference/Ctpsv.cpp"
#undef Ctpsv
#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
#if defined(MPLAPACK_BUILD_WITH_DD)
#include "../../dd/cuda/Rgemm_gpu_cuda.cu"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include "../../td/cuda/Rgemm_gpu_cuda.cu"
#else
#include "../../qd/cuda/Rgemm_gpu_cuda.cu"
#endif
#elif !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
#include <cuda_runtime.h>
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
#define PREFIX "MPLAPACK_MPFR"
#elif defined(MPLAPACK_BUILD_WITH_DD)
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define PREFIX "MPLAPACK_TD"
#else
#define PREFIX "MPLAPACK_QD"
#endif
namespace {
long checks=0,failures=0;
void check(bool ok,const char *name) { ++checks;if(!ok) { if(failures<30) std::printf("FAIL %s\n",name);++failures; } }
bool same(const REAL &a,const REAL &b) {
#ifdef MPLAPACK_BUILD_WITH_MPFR
    return mpfr_equal_p(a.get_mpfr_t(),b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t())==mpfr_signbit(b.get_mpfr_t());
#else
    return std::memcmp(&a,&b,sizeof a)==0;
#endif
}
bool same(const COMPLEX &a,const COMPLEX &b) { return same(a.real(),b.real()) && same(a.imag(),b.imag()); }
template<class H> bool same(const std::vector<H> &a,const std::vector<H> &b) { for(size_t i=0;i<a.size();++i) if(!same(a[i],b[i])) return false;return true; }
#define CALLS(P,H) \
bool entry(int kind,const char *u,const char *t,const char *d,long n,long k,H *a,long lda,H *x,long inc) { \
    if(kind==0) return MPLAPACK_CUDA_ENTRY(P##trsv)(u,t,d,n,a,lda,x,inc); \
    if(kind==1) return MPLAPACK_CUDA_ENTRY(P##tbsv)(u,t,d,n,k,a,lda,x,inc); \
    return MPLAPACK_CUDA_ENTRY(P##tpsv)(u,t,d,n,a,x,inc); } \
void reference(int kind,const char *u,const char *t,const char *d,long n,long k,H *a,long lda,H *x,long inc) { \
    if(kind==0) ref_##P##trsv(u,t,d,n,a,lda,x,inc); \
    else if(kind==1) ref_##P##tbsv(u,t,d,n,k,a,lda,x,inc); \
    else ref_##P##tpsv(u,t,d,n,a,x,inc); } \
void dispatch(int kind,const char *u,const char *t,const char *d,long n,long k,H *a,long lda,H *x,long inc) { \
    if(kind==0) P##trsv(u,t,d,n,a,lda,x,inc); \
    else if(kind==1) P##tbsv(u,t,d,n,k,a,lda,x,inc); \
    else P##tpsv(u,t,d,n,a,x,inc); }
CALLS(R,REAL)
CALLS(C,COMPLEX)
#undef CALLS
template<class H> H number(long i,long j) {
    REAL r=REAL((i*17+j*13)%19-9)/REAL(997.0);
    if constexpr(std::is_same<H,COMPLEX>::value) return H(r,REAL((i*11+j*7)%13-6)/REAL(991.0));
    else return r;
}
long at(int kind,bool upper,long n,long k,long lda,long i,long j) {
    if(kind==0) return i+j*lda;
    if(kind==1) return (upper?k+i-j:i-j)+j*lda;
    return upper?j*(j+1)/2+i:j*n-j*(j-1)/2+i-j;
}
template<class H> void run() {
    for(long n:{0L,1L,7L,64L,200L}) for(char u:{'U','L'}) for(char t:{'N','T','C'}) for(char d:{'U','N'})
    for(long inc:{-3L,-1L,1L,2L}) for(int kind=0;kind<3;++kind) for(long k: {0L,3L,203L}) {
        if(kind!=1 && k!=0) continue;
#ifdef MPLAPACK_BUILD_WITH_MPFR
        if(mpfrxx::default_precision_bits()>4096 && n>7) continue;
#endif
        long lda=(kind==1?k+1:std::max(1L,n))+3;
        long count=kind==2?n*(n+1)/2:lda*n;
        std::vector<H> a(count+2,H(77.0)),x(3+std::max(1L,n)*std::labs(inc),H(91.0));
        for(long j=0;j<n;++j) for(long i=0;i<n;++i) if((u=='U'?i<=j:i>=j) && (kind!=1 || std::labs(i-j)<=k))
            a[1+at(kind,u=='U',n,k,lda,i,j)]=i==j?H(1.3)+number<H>(i,j):number<H>(i,j);
        // Unit diagonal must not be used (deliberately zero).
        if(d=='U') for(long j=0;j<n;++j) a[1+at(kind,u=='U',n,k,lda,j,j)]=H(0.0);
        for(long i=0;i<n;++i) x[1+i*std::labs(inc)]=i%7?number<H>(i,2):H(-0.0);
        auto before=a,expected=x,original=x;
        reference(kind,&u,&t,&d,n,k,a.data()+1,lda,expected.data()+1,inc);
        bool ok=entry(kind,&u,&t,&d,n,k,a.data()+1,lda,x.data()+1,inc);
        check(ok==(n>0),"status");
        if(!same(x,expected) && failures<30) {
            std::printf("identity case complex=%d n=%ld %c%c%c inc=%ld layout=%d k=%ld\n",int(std::is_same<H,COMPLEX>::value),n,u,t,d,inc,kind,k);
            for(size_t z=0;z<x.size();++z) if(!same(x[z],expected[z])) { std::printf("first index=%zu\n",z);break; }
        }
        check(same(x,expected),"reference identity/guards");check(same(a,before),"matrix unchanged");
        if(n<=7 && inc==1) {
            setenv(PREFIX "_BLAS_CUDA_L12","0",1);setenv(PREFIX "_CUDA_L12","0",1);
            x=original;check(!entry(kind,&u,&t,&d,n,k,a.data()+1,lda,x.data()+1,inc) && same(x,original),"disabled untouched");
            dispatch(kind,&u,&t,&d,n,k,a.data()+1,lda,x.data()+1,inc);check(same(x,expected),"CPU fallback unchanged");
            setenv(PREFIX "_BLAS_CUDA_L12","1",1);setenv(PREFIX "_CUDA_L12","1",1);
        }
    }
    for(int kind=0;kind<3;++kind) {
        std::vector<H> a(64,H(0.0)),x(16,H(0.0)),before=x;
        // BLAS does not precheck singularity. Zero RHS/N skips a zero pivot.
        check(entry(kind,"U","N","N",1,0,a.data(),8,x.data(),1),"singular zero RHS skip");
        check(same(x,before),"zero skip identity");
        x[0]=H(1.0);before=x;
        check(!entry(kind,"U","N","N",1,0,a.data(),8,x.data(),1) && same(x,before),"division failure unpublished");
        for(const char *bad:{"?", ""}) {
            check(!entry(kind,bad,"N","N",7,0,a.data(),8,x.data(),1) && same(x,before),"bad uplo");
            check(!entry(kind,"U",bad,"N",7,0,a.data(),8,x.data(),1) && same(x,before),"bad trans");
            check(!entry(kind,"U","N",bad,7,0,a.data(),8,x.data(),1) && same(x,before),"bad diag");
        }
        check(!entry(kind,nullptr,"N","N",7,0,a.data(),8,x.data(),1),"null flag");
        check(!entry(kind,"U","N","N",-1,0,a.data(),8,x.data(),1) && same(x,before),"negative n");
        check(!entry(kind,"U","N","N",7,0,a.data(),8,x.data(),0) && same(x,before),"zero stride");
        check(!entry(kind,"U","N","N",7,0,a.data(),8,x.data(),std::numeric_limits<long>::min()),"stride overflow");
        check(!entry(kind,"U","N","N",std::numeric_limits<long>::max(),0,a.data(),8,x.data(),1),"size overflow");
        if(kind!=2) check(!entry(kind,"U","N","N",7,7,a.data(),0,x.data(),1),"lda invalid");
        if(kind==1) {
            check(!entry(kind,"U","N","N",7,-1,a.data(),8,x.data(),1),"band invalid");
            check(!entry(kind,"U","N","N",7,std::numeric_limits<long>::max(),a.data(),8,x.data(),1),"band overflow");
        }
        a[0]=H(std::numeric_limits<double>::quiet_NaN());
        check(!entry(kind,"U","N","N",1,0,a.data(),8,x.data(),1) && same(x,before),"nonfinite input rejected");
        check(entry(kind,"U","N","U",1,0,a.data(),8,x.data(),1) && same(x,before),"unit NaN diagonal ignored");
    }
}
}
int main() {
#if !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
    int devices=0;if(cudaGetDeviceCount(&devices) || !devices) return 77;
#endif
    setenv(PREFIX "_BLAS_CUDA_MIN","0",1);setenv(PREFIX "_BLAS_CUDA_L12","1",1);
    setenv(PREFIX "_CUDA_MIN_MNK","0",1);setenv(PREFIX "_CUDA_L12","1",1);
    run<REAL>();
#ifdef MPLAPACK_BUILD_WITH_MPFR
    if(mpfrxx::default_precision_bits()<=4096)
#endif
    run<COMPLEX>();
#ifdef MPLAPACK_BUILD_WITH_MPFR
    {
        REAL a=1.0,x=2.0,original=x;
        long p=mpfr_get_prec(x.get_mpfr_t());
        mpfr_prec_round(x.get_mpfr_t(),p+1,MPFR_RNDN);
        check(!MPLAPACK_CUDA_ENTRY(Rtrsv)("U","N","U",1,&a,1,&x,1) && x==original,"precision rejected");
        mpfr_prec_round(x.get_mpfr_t(),p,MPFR_RNDN);
        auto rounding=mpfr_get_default_rounding_mode();mpfr_set_default_rounding_mode(MPFR_RNDZ);
        check(!MPLAPACK_CUDA_ENTRY(Rtrsv)("U","N","U",1,&a,1,&x,1) && same(x,original),"rounding rejected");
        mpfr_set_default_rounding_mode(rounding);
    }
#endif
    std::printf("triangular GPU solves: %ld checks, %ld failures\n",checks,failures);return failures?1:0;
}
