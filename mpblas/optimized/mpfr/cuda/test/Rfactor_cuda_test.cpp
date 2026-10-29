// Direct-entry contract and whole padded matrix oracle for MPFR real factors.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include <cuda_runtime_api.h>
#endif
#define Rgetrf2 factor_oracle_getrf2
#include "../../../../../mplapack/reference/Rgetrf2.cpp"
#define Rgetrf factor_oracle_getrf
#include "../../../../../mplapack/reference/Rgetrf.cpp"
#undef Rgetrf
#undef Rgetrf2
#define Rpotrf2 factor_oracle_potrf2
#include "../../../../../mplapack/reference/Rpotrf2.cpp"
#define Rpotrf factor_oracle_potrf
#include "../../../../../mplapack/reference/Rpotrf.cpp"
#undef Rpotrf
#undef Rpotrf2
namespace {
long checks=0,errors=0,elements=0,zero_differences=0;
long routine_cases[2]{},routine_errors[2]{},routine_elements[2]{},routine_zeros[2]{};
gmp_randstate_t rng;
bool same(const REAL &a,const REAL &b) {
    auto x=a.get_mpfr_t(),y=b.get_mpfr_t();
    ++elements;
    if(mpfr_get_prec(x)!=mpfr_get_prec(y)) return false;
    if(mpfr_zero_p(x) && mpfr_zero_p(y)) {
        if(mpfr_signbit(x)!=mpfr_signbit(y)) { ++zero_differences; return false; }
        return true;
    }
    if(mpfr_nan_p(x) && mpfr_nan_p(y)) return true;
    if(mpfr_inf_p(x) || mpfr_inf_p(y)) return mpfr_inf_p(x) && mpfr_inf_p(y) && mpfr_signbit(x)==mpfr_signbit(y);
    if(!mpfr_regular_p(x) || !mpfr_regular_p(y) || mpfr_signbit(x)!=mpfr_signbit(y) || mpfr_get_exp(x)!=mpfr_get_exp(y)) return false;
    return !std::memcmp(mpfr_custom_get_significand(x),mpfr_custom_get_significand(y),((mpfr_get_prec(x)+63)/64)*sizeof(mp_limb_t));
}
void check(bool pass,const char *what,long m=0,long n=0,long kind=0) {
    ++checks; if(!pass) { ++errors; std::fprintf(stderr,"FAIL %s m=%ld n=%ld kind=%ld\n",what,m,n,kind); }
}
bool backward(bool chol,bool upper,long m,long n,long ld,const std::vector<REAL> &a,
              std::vector<REAL> original,const std::vector<INTEGER> &piv) {
    REAL norm=0,residual=0;
    for(long j=0;j<n;++j) { REAL sum=0; for(long i=0;i<m;++i) sum+=abs(original[i+j*ld]); if(sum>norm) norm=sum; }
    if(!chol) for(long k=0;k<std::min(m,n);++k) for(long j=0;j<n;++j)
        std::swap(original[k+j*ld],original[piv[k]-1+j*ld]);
    for(long j=0;j<n;++j) {
        REAL sum=0;
        for(long i=0;i<m;++i) {
            REAL product=0;
            if(chol) for(long k=0;k<=std::min(i,j);++k)
                product+=upper?a[k+i*ld]*a[k+j*ld]:a[i+k*ld]*a[j+k*ld];
            else for(long k=0;k<std::min(m,n);++k) {
                REAL l=i==k?REAL(1):(i>k?a[i+k*ld]:REAL(0));
                REAL u=k<=j?a[k+j*ld]:REAL(0); product+=l*u;
            }
            sum+=abs(product-original[i+j*ld]);
        }
        if(sum>residual) residual=sum;
    }
    return residual/(REAL(std::max(m,n))*norm*Rlamch("E"))<=REAL(30);
}
void trial(bool chol,bool upper,long m,long n,int kind) {
    long previous_errors=errors,previous_elements=elements,previous_zeros=zero_differences;
    ++routine_cases[chol?1:0];
    long lda=m+3; std::vector<REAL> a(lda*n),b; std::vector<INTEGER> p(std::min(m,n)+3,-71),q=p;
    for(long j=0;j<n;++j) for(long i=0;i<lda;++i) {
        REAL &v=a[i+j*lda];
        if(i>=m) { v=731; continue; }
        mpfr_urandomb(v.get_mpfr_t(),rng);
        v=(v-REAL(0.5))/REAL(8);
        if(i==j) v+=REAL(n+2);
        if(i!=j && (i+j)%7==0) mpfr_set_zero(v.get_mpfr_t(),(i+j)%2?1:-1);
    }
    if(chol) {
        for(long j=0;j<n;++j) for(long i=j+1;i<m;++i) a[j+i*lda]=a[i+j*lda];
        if(kind) a[(kind==1?0:n/2)*(lda+1)]=-1;
    } else {
        if(kind==1) for(long i=0;i<m;++i) a[i+(n/2)*lda]=0;
        if(kind==2 && n>1) for(long i=0;i<m;++i) a[i+lda]=a[i];
        if(kind==3) for(long j=0;j<n;++j) for(long i=0;i<m;++i) a[i+j*lda]=i+j+1;
        if(kind==4) for(long i=0;i<m;++i) a[i]=i%2?-1:1;
        if(kind==5) for(long i=0;i<m;++i) a[i]=REAL(1)/(i+1);
        if(kind==5) for(long i=0;i<m;++i) mpfr_mul_2si(a[i].get_mpfr_t(),a[i].get_mpfr_t(),-1073741823L,MPFR_RNDN);
    }
    b=a; auto original=a; INTEGER info=-99,ref=-98; bool ran;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t stack_before=0,stack_after=0; cudaDeviceGetLimit(&stack_before,cudaLimitStackSize);
#endif
    if(chol) { factor_oracle_potrf(upper?"U":"L",n,b.data(),lda,ref); ran=Rpotrf_mpfr_cuda(upper?"U":"L",n,a.data(),lda,info); }
    else { factor_oracle_getrf(m,n,b.data(),lda,q.data(),ref); ran=Rgetrf_mpfr_cuda(m,n,a.data(),lda,p.data(),info); }
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    cudaDeviceGetLimit(&stack_after,cudaLimitStackSize);
    check(stack_before==stack_after,"stack restored",m,n,kind);
#endif
    check(ran && info==ref && p==q,chol?"potrf status":"getrf status",m,n,kind);
    long bad=0; for(size_t i=0;i<a.size();++i) if(!same(a[i],b[i])) {
        if(!bad && errors<10) { std::fprintf(stderr,"first diff idx=%zu ",i); mpfr_out_str(stderr,16,0,a[i].get_mpfr_t(),MPFR_RNDN); std::fputs(" vs ",stderr); mpfr_out_str(stderr,16,0,b[i].get_mpfr_t(),MPFR_RNDN); std::fputc('\n',stderr); } ++bad;
    }
    check(!bad,chol?"potrf A":"getrf A",m,n,kind);
    if(!kind && n<=65 && ran && !info) check(backward(chol,upper,m,n,lda,a,original,p),"backward error",m,n,kind);
    routine_errors[chol?1:0]+=errors-previous_errors;
    routine_elements[chol?1:0]+=elements-previous_elements;
    routine_zeros[chol?1:0]+=zero_differences-previous_zeros;
}
void rejected() {
    std::vector<REAL> a(12,REAL(3)),b=a; INTEGER p[2]={-7,-8},info=-9;
    auto untouched=[&] { bool result=info==-9 && p[0]==-7 && p[1]==-8; for(size_t i=0;i<a.size();++i) result=same(a[i],b[i])&&result; return result; };
    check(!Rgetrf_mpfr_cuda(0,2,a.data(),3,p,info)&&untouched(),"quick return");
    check(!Rgetrf_mpfr_cuda(3,0,a.data(),3,p,info)&&untouched(),"zero columns");
    check(!Rpotrf_mpfr_cuda("L",0,a.data(),3,info)&&untouched(),"potrf quick return");
    check(!Rgetrf_mpfr_cuda(2,2,a.data(),3,nullptr,info)&&untouched(),"null pivots");
    check(!Rpotrf_mpfr_cuda(nullptr,2,a.data(),3,info)&&untouched(),"null flag");
    check(!Rgetrf_mpfr_cuda(3,2,a.data(),2,p,info)&&untouched(),"bad lda");
    check(!Rpotrf_mpfr_cuda("?",2,a.data(),3,info)&&untouched(),"bad uplo");
    check(!Rgetrf_mpfr_cuda(std::numeric_limits<INTEGER>::max(),2,a.data(),std::numeric_limits<INTEGER>::max(),p,info)&&untouched(),"size overflow");
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    check(!Rgetrf_mpfr_cuda(1048576,1048576,a.data(),1048576,p,info)&&untouched(),"memory budget");
#endif
    setenv("MPLAPACK_MPFR_LAPACK_CUDA","0",1);
    check(!Rgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"disabled");
    setenv("MPLAPACK_MPFR_LAPACK_CUDA","1",1); setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","3",1);
    check(!Rgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"threshold");
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    mpfr_set_default_rounding_mode(MPFR_RNDZ);
    check(!Rgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"rounding");
    mpfr_set_default_rounding_mode(MPFR_RNDN);
    auto emin=mpfr_get_emin(); mpfr_set_emin(emin+1);
    check(!Rpotrf_mpfr_cuda("L",2,a.data(),3,info)&&untouched(),"exponent range");
    mpfr_set_emin(emin);
    mpfr_set_prec(a[0].get_mpfr_t(),mpfr_get_prec(a[0].get_mpfr_t())+1);
    mpfr_set_ui(a[0].get_mpfr_t(),3,MPFR_RNDN);
    mpfr_set_prec(b[0].get_mpfr_t(),mpfr_get_prec(a[0].get_mpfr_t()));
    mpfr_set(b[0].get_mpfr_t(),a[0].get_mpfr_t(),MPFR_RNDN);
    check(!Rgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"mixed precision");
    mpfr_set_prec(a[0].get_mpfr_t(),mpfrxx::default_precision_bits());
    mpfr_set_prec(b[0].get_mpfr_t(),mpfrxx::default_precision_bits()); a[0]=3;
    mpfr_set_inf(a[0].get_mpfr_t(),1); b=a;
    check(!Rgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"nonfinite");
    a.assign(12,REAL(1)); b=a;
    mpfr_set_ui_2exp(a[1].get_mpfr_t(),1,mpfr_get_emax()-1,MPFR_RNDN); b=a;
    check(!Rpotrf_mpfr_cuda("L",2,a.data(),3,info)&&untouched(),"intermediate overflow");
}
void natural_defaults() {
    mpfrxx::set_default_precision_bits(512);
    unsetenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN");
    setenv("MPLAPACK_MPFR_CUDA","1",1); setenv("MPLAPACK_MPFR_LAPACK_CUDA","1",1);
    setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME","0",1);
    for(bool chol:{false,true}) {
        long n=chol?2048:512;
        std::vector<REAL> a(n*n,REAL(0)); std::vector<INTEGER> p(n,-11);
        for(long i=0;i<n;++i) a[i+i*n]=2;
        INTEGER info=-9;
        bool ran=chol?Rpotrf_mpfr_cuda("L",n,a.data(),n,info):Rgetrf_mpfr_cuda(n,n,a.data(),n,p.data(),info);
        check(ran && !info,"natural default eligible",n,n);
        REAL diagonal=2,zero=0;
        if(chol) mpfr_sqrt(diagonal.get_mpfr_t(),diagonal.get_mpfr_t(),MPFR_RNDN);
        bool correct=true;
        for(long j=0;j<n;++j) for(long i=0;i<n;++i) correct=same(a[i+j*n],i==j?diagonal:zero)&&correct;
        if(!chol) for(long i=0;i<n;++i) correct=correct && p[i]==i+1;
        check(correct,"natural default result",n,n);
        auto previous=p; info=-9;
        ran=chol?Rpotrf_mpfr_cuda("L",n-1,a.data(),n,info):Rgetrf_mpfr_cuda(n-1,n-1,a.data(),n,p.data(),info);
        correct=!ran && info==-9 && p==previous;
        for(long j=0;j<n;++j) for(long i=0;i<n;++i) correct=same(a[i+j*n],i==j?diagonal:zero)&&correct;
        check(correct,"below natural default untouched",n-1,n-1);
    }
    std::printf("natural defaults checks=%ld elements=%ld errors=%ld signed_zero_diffs=%ld\n",checks,elements,errors,zero_differences);
}
}
int main(int argc,char **argv) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
    if(argc>1 && !std::strcmp(argv[1],"--defaults")) { natural_defaults(); return errors?1:0; }
    bool quick=argc>2 && !std::strcmp(argv[2],"--quick");
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    gmp_randinit_default(rng); gmp_randseed_ui(rng,7);
    bool selected=false;
    for(long prec:{256,512,768,1024,2048,200,64,333,4096}) {
        if(argc>1 && std::atol(argv[1])!=prec) continue;
        if(argc==1 && (prec==64 || prec==333 || prec==4096)) continue;
        selected=true;
        mpfrxx::set_default_precision_bits(prec);
        long first=checks,ec=elements,zc=zero_differences,fail=errors;
        long cases_before[2]={routine_cases[0],routine_cases[1]};
        long errors_before[2]={routine_errors[0],routine_errors[1]};
        long elements_before[2]={routine_elements[0],routine_elements[1]};
        long zeros_before[2]={routine_zeros[0],routine_zeros[1]};
        for(long n:{1,2,3,31,32,33,63,64,65,127,128,129,131}) {
            if(quick && n>65) continue;
            for(int kind=0;kind<6;++kind) { trial(false,false,n,n,kind); trial(false,false,n+3,n,kind); trial(false,false,n,n+3,kind); }
            for(bool u:{false,true}) for(int k=0;k<3;++k) trial(true,u,n,n,k);
        }
        // Odd sizes through 599, with padded leading dimensions and both
        // rectangular orientations. Full edge-case families run above.
        for(long n:{257,599}) {
            if(quick) continue;
            trial(false,false,n,n,0); trial(false,false,n+3,n,0); trial(false,false,n,n+3,0);
            trial(true,false,n,n,0); trial(true,true,n,n,0);
        }
        if(quick) { // Three blocks exercise actual panel/far-update overlap.
            trial(false,false,134,131,0); trial(false,false,131,134,0);
            trial(true,false,131,131,0); trial(true,true,131,131,0);
        }
        setenv("MPLAPACK_MPFR_LAPACK_FAST","1",1);
        trial(false,false,35,33,0); trial(true,false,33,33,0); trial(true,true,33,33,0);
        unsetenv("MPLAPACK_MPFR_LAPACK_FAST");
        rejected();
        for(int r=0;r<2;++r) std::printf("p=%ld routine=%s cases=%ld elements=%ld errors=%ld signed_zero_diffs=%ld\n",
            prec,r?"Rpotrf":"Rgetrf",routine_cases[r]-cases_before[r],routine_elements[r]-elements_before[r],
            routine_errors[r]-errors_before[r],routine_zeros[r]-zeros_before[r]);
        std::printf("p=%ld checks=%ld elements=%ld errors=%ld signed_zero_diffs=%ld\n",prec,checks-first,elements-ec,errors-fail,zero_differences-zc); std::fflush(stdout);
    }
    gmp_randclear(rng); return !selected?2:errors?1:0;
}
