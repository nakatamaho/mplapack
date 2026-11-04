// Focused inverse tests: direct entries, unchanged reference CPU oracles.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <limits>
#include <string>
#include <type_traits>
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include <cuda_runtime_api.h>
#endif
#define Rtrtri oracle_Rtrtri
#define Ctrtri oracle_Ctrtri
#define Rlauum oracle_Rlauum
#define Clauum oracle_Clauum
#define Rpotri oracle_Rpotri
#define Cpotri oracle_Cpotri
#define Rgetri oracle_Rgetri
#define Cgetri oracle_Cgetri
#include "../../../../../mplapack/reference/Rtrtri.cpp"
#include "../../../../../mplapack/reference/Ctrtri.cpp"
#include "../../../../../mplapack/reference/Rlauum.cpp"
#include "../../../../../mplapack/reference/Clauum.cpp"
#include "../../../../../mplapack/reference/Rpotri.cpp"
#include "../../../../../mplapack/reference/Cpotri.cpp"
#include "../../../../../mplapack/reference/Rgetri.cpp"
#include "../../../../../mplapack/reference/Cgetri.cpp"
#undef Rtrtri
#undef Ctrtri
#undef Rlauum
#undef Clauum
#undef Rpotri
#undef Cpotri
#undef Rgetri
#undef Cgetri
namespace {
long checks=0,errors=0,elements=0,calls=0;
bool dispatch=false;
void check(bool ok,const char *what,long n=0,int kind=0) {
    ++checks; if(!ok) { ++errors; if(errors<30) std::fprintf(stderr,"FAIL %s n=%ld kind=%d\n",what,n,kind); }
}
bool same_part(mpfr_srcptr a,mpfr_srcptr b) {
    ++elements;
    return mpfr_get_prec(a)==mpfr_get_prec(b) && mpfr_signbit(a)==mpfr_signbit(b) &&
        ((mpfr_nan_p(a) && mpfr_nan_p(b)) || mpfr_equal_p(a,b));
}
bool same(const REAL &a,const REAL &b) { return same_part(a.get_mpfr_t(),b.get_mpfr_t()); }
bool same(const COMPLEX &a,const COMPLEX &b) {
    return same_part(mpc_realref(a.mpc_data()),mpc_realref(b.mpc_data())) && same_part(mpc_imagref(a.mpc_data()),mpc_imagref(b.mpc_data()));
}
template<class T> bool same(const std::vector<T> &a,const std::vector<T> &b) {
    bool ok=a.size()==b.size();for(size_t i=0;i<a.size() && i<b.size();++i) if(!same(a[i],b[i])) ok=false;return ok;
}
REAL scalar(long i,long j) { return REAL((i*17+j*23)%31-15)/REAL(37); }
REAL input(REAL *,long i,long j) { return i==j?REAL(2):scalar(i,j); }
COMPLEX input(COMPLEX *,long i,long j) { return COMPLEX(i==j?REAL(2):scalar(i,j),scalar(j,i)/REAL(11)); }
void zero(REAL &x) { mpfr_set_zero(x.get_mpfr_t(),-1); }
void zero(COMPLEX &x) { mpfr_set_zero(mpc_realref(x.mpc_data()),-1);mpfr_set_zero(mpc_imagref(x.mpc_data()),-1); }
void huge(REAL &x) { mpfr_set_ui_2exp(x.get_mpfr_t(),1,500000000,MPFR_RNDN); }
void huge(COMPLEX &x) { mpfr_set_ui_2exp(mpc_realref(x.mpc_data()),1,500000000,MPFR_RNDN); }
void nonfinite(REAL &x) { mpfr_set_nan(x.get_mpfr_t()); }
void nonfinite(COMPLEX &x) { mpfr_set_nan(mpc_realref(x.mpc_data())); }
void different(REAL &x) { mpfr_prec_round(x.get_mpfr_t(),mpfr_get_prec(x.get_mpfr_t())+1,MPFR_RNDN); }
void different(COMPLEX &x) { mpfr_prec_round(mpc_imagref(x.mpc_data()),mpfr_get_prec(mpc_imagref(x.mpc_data()))+1,MPFR_RNDN); }
template<class T> struct api;
#define API(T,P) \
template<> struct api<T> { \
static bool run(int kind,const char *u,const char *d,INTEGER n,T *a,INTEGER ld,INTEGER *p,T *w,INTEGER lw,INTEGER &info) { \
    if(dispatch) { switch(kind) {case 0:P##trtri(u,d,n,a,ld,info);break;case 1:P##lauum(u,n,a,ld,info);break;case 2:P##potri(u,n,a,ld,info);break;default:P##getri(n,a,ld,p,w,lw,info);} return true; } \
    switch(kind) { case 0:return P##trtri_mpfr_cuda(u,d,n,a,ld,info);case 1:return P##lauum_mpfr_cuda(u,n,a,ld,info);case 2:return P##potri_mpfr_cuda(u,n,a,ld,info);default:return P##getri_mpfr_cuda(n,a,ld,p,w,lw,info); } } \
static void ref(int kind,const char *u,const char *d,INTEGER n,T *a,INTEGER ld,INTEGER *p,T *w,INTEGER lw,INTEGER &info) { \
    switch(kind) {case 0:oracle_##P##trtri(u,d,n,a,ld,info);break;case 1:oracle_##P##lauum(u,n,a,ld,info);break;case 2:oracle_##P##potri(u,n,a,ld,info);break;default:oracle_##P##getri(n,a,ld,p,w,lw,info);} } \
};
API(REAL,R)
API(COMPLEX,C)
#undef API
template<class T> void trial(long n,int kind,bool upper,bool unit,int pattern,long columns) {
    long ld=n+3,nb=iMlaenv(1,std::is_same<T,REAL>::value?"Rgetri":"Cgetri"," ",n,-1,-1,-1);
    long lw=n*(columns?columns:nb)+7;
    std::vector<T> a(ld*n,T(19)),w(lw,T(17));std::vector<INTEGER> piv(n,n),p0=piv;
    for(long j=0;j<n;++j) for(long i=0;i<n;++i) {
        a[i+j*ld]=input((T*)nullptr,i,j)/T(i==j?1:n);
        if(pattern==1 && i!=j) zero(a[i+j*ld]);
    }
    if(pattern==2) a[n/2+n/2*ld]=T(0);
    auto r=a,wr=w;INTEGER info=-77,expected=-99;
    api<T>::ref(kind,upper?"U":"L",unit?"U":"N",n,r.data(),ld,piv.data(),wr.data(),lw,expected);
    ++calls;
    bool ok=api<T>::run(kind,upper?"U":"L",unit?"U":"N",n,a.data(),ld,piv.data(),w.data(),lw,info);
    check(ok,"accepted",n,kind);check(info==expected,"info",n,kind);
    check(same(a,r),"matrix/signs/padding",n,kind);check(same(w,wr),"whole workspace",n,kind);check(piv==p0,"input pivots",n,kind);
}
template<class T> void contracts() {
    const long n=3,ld=6;std::vector<T> a(ld*n,T(2)),w(200,T(17));std::vector<INTEGER> p={1,2,3};
    for(int kind=0;kind<4;++kind) for(int test=0;test<16;++test) {
        a.assign(ld*n,T(2));w.assign(200,T(17));p={1,2,3};INTEGER info=-77;
        long nn=n,ll=ld,lw=200;const char *u="U",*d="N";T *ap=a.data(),*wp=w.data();INTEGER *pp=p.data();
        if(test==0) nn=0;
        if(test==1) nn=-1;
        if(test==2) ll=2;
        if(test==3) ap=nullptr;
        if(test==4) { if(kind==3) lw=-1;else u="X"; }
        if(test==5) { if(kind==3) lw=2;else nonfinite(a[0]); }
        if(test==6) { if(kind==3) pp=nullptr;else different(a[0]); }
        if(test==7) { if(kind==3) wp=nullptr;else huge(a[0]); }
        if(test==8) { if(kind==3) p[1]=0;else setenv("MPLAPACK_MPFR_CUDA","0",1); }
        if(test==9) setenv("MPLAPACK_MPFR_LAPACK_CUDA","0",1);
        if(test==10) setenv("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN","1000",1);
        if(test==11) {unsetenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN");unsetenv("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN");}
        if(test==12) { if(kind==3) different(w[1]);else huge(a[1+ld]); }
        if(test==13) { if(kind==3) nonfinite(a[1]);else nonfinite(a[0]); }
        if(test==14) { nn=ll=std::numeric_limits<INTEGER>::max(); }
        if(test==15) mpfr_set_default_rounding_mode(MPFR_RNDD);
        auto before=a,bw=w;auto bp=p;
        bool ok=api<T>::run(kind,u,d,nn,ap,ll,pp,wp,lw,info);
        check(!ok && info==-77 && same(a,before) && same(w,bw) && p==bp,"reject/atomic publication",test,kind);
        setenv("MPLAPACK_MPFR_CUDA","1",1);setenv("MPLAPACK_MPFR_LAPACK_CUDA","1",1);
        setenv("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN","0",1);mpfr_set_default_rounding_mode(MPFR_RNDN);
    }
}
}
int main(int argc,char **argv) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int devices=0;if(cudaGetDeviceCount(&devices) || !devices) return 77;
    size_t old_stack=0;cudaDeviceGetLimit(&old_stack,cudaLimitStackSize);
#endif
    setenv("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN","0",1);
    // Each process gets its precision through MPFRXX_DEFAULT_PRECISION_BITS,
    // including the CPU BLAS worker threads; no global precision mutation.
    long p=mpfrxx::default_precision_bits();
    const std::string mode=argc>1?argv[1]:"";
    dispatch=mode=="--fallback";
    bool real_only=mode=="--real-small" || mode=="--real-blocked";
    if(dispatch) setenv("MPLAPACK_MPFR_CUDA","0",1);
    else { contracts<REAL>();if(!real_only) contracts<COMPLEX>(); }
    bool small=mode=="--small" || mode=="--real-small";
    const auto sizes=mode=="--real-blocked"?std::vector<long>{65}:small?std::vector<long>{1,3}:std::vector<long>{1,2,3,63,64,65,131};
    for(long n:sizes) {
        for(int kind=0;kind<3;++kind) for(bool upper:{false,true}) for(int pattern=0;pattern<3;++pattern) {
            trial<REAL>(n,kind,upper,false,pattern,0);
            if(!real_only) trial<COMPLEX>(n,kind,upper,false,pattern,0);
            if(kind==0) {
                trial<REAL>(n,kind,upper,true,pattern,0);
                if(!real_only) trial<COMPLEX>(n,kind,upper,true,pattern,0);
            }
        }
        for(int pattern=0;pattern<3;++pattern) for(long columns:{0L,1L,2L}) {
            trial<REAL>(n,3,true,false,pattern,columns);
            if(!real_only) trial<COMPLEX>(n,3,true,false,pattern,columns);
        }
    }
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t now=0;check(!cudaDeviceGetLimit(&now,cudaLimitStackSize) && now==old_stack,"stack restored");
    check(cudaDeviceSynchronize()==cudaSuccess,"device sync");
#endif
    std::printf("MPFR inverse p=%ld: %ld %s calls, %ld checks, %ld component comparisons, %ld failures\n",p,calls,dispatch?"fallback":"direct",checks,elements,errors);
    return errors?1:0;
}
