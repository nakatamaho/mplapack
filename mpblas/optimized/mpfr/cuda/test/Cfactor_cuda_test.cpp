// Direct-entry complex factor oracle; padded matrices and signs are strict.
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
#define Cgetrf2 complex_oracle_getrf2
#include "../../../../../mplapack/reference/Cgetrf2.cpp"
#define Cgetrf complex_oracle_getrf
#include "../../../../../mplapack/reference/Cgetrf.cpp"
#undef Cgetrf
#undef Cgetrf2
#define Cpotrf2 complex_oracle_potrf2
#include "../../../../../mplapack/reference/Cpotrf2.cpp"
#define Cpotrf complex_oracle_potrf
#include "../../../../../mplapack/reference/Cpotrf.cpp"
#undef Cpotrf
#undef Cpotrf2
namespace {
long checks=0,errors=0,elements=0,zeros=0,cases[2]{};
gmp_randstate_t rng;
void check(bool pass,const char *what,long m=0,long n=0,long k=0) {
    ++checks; if(!pass) { ++errors; if(errors<30) std::fprintf(stderr,"FAIL %s m=%ld n=%ld kind=%ld\n",what,m,n,k); }
}
bool same_part(mpfr_srcptr a,mpfr_srcptr b) {
    ++elements;
    if(mpfr_get_prec(a)!=mpfr_get_prec(b)) return false;
    if(mpfr_zero_p(a) && mpfr_zero_p(b)) {
        if(mpfr_signbit(a)!=mpfr_signbit(b)) { ++zeros; return false; }
        return true;
    }
    if(mpfr_nan_p(a) && mpfr_nan_p(b)) return true;
    if(!mpfr_equal_p(a,b) || mpfr_signbit(a)!=mpfr_signbit(b)) return false;
    return true;
}
bool same(const COMPLEX &a,const COMPLEX &b) {
    auto x=a.mpc_data(),y=b.mpc_data();
    bool re=same_part(mpc_realref(x),mpc_realref(y));
    bool im=same_part(mpc_imagref(x),mpc_imagref(y)); return re && im;
}
bool backward(bool chol,bool upper,long m,long n,long ld,const std::vector<COMPLEX> &a,
              std::vector<COMPLEX> original,const std::vector<INTEGER> &piv) {
    REAL norm=0,residual=0;
    for(long j=0;j<n;++j) { REAL sum=0; for(long i=0;i<m;++i) sum+=abs(original[i+j*ld]); if(sum>norm) norm=sum; }
    if(!chol) for(long k=0;k<std::min(m,n);++k) for(long j=0;j<n;++j)
        std::swap(original[k+j*ld],original[piv[k]-1+j*ld]);
    for(long j=0;j<n;++j) {
        REAL sum=0;
        for(long i=0;i<m;++i) {
            COMPLEX product=0;
            if(chol) for(long k=0;k<=std::min(i,j);++k) {
                if(upper) product+=conj(a[k+i*ld])*a[k+j*ld];
                else product+=a[i+k*ld]*conj(a[j+k*ld]);
            }
            else for(long k=0;k<std::min(m,n);++k) {
                COMPLEX l=i==k?COMPLEX(1):(i>k?a[i+k*ld]:COMPLEX(0));
                COMPLEX u=k<=j?a[k+j*ld]:COMPLEX(0); product+=l*u;
            }
            sum+=abs(product-original[i+j*ld]);
        }
        if(sum>residual) residual=sum;
    }
    return residual/(REAL(std::max(m,n))*norm*Rlamch("E"))<=REAL(30);
}
void trial(bool chol,bool upper,long m,long n,int kind) {
    ++cases[chol]; long lda=m+3;
    std::vector<COMPLEX> a(lda*n),b;
    std::vector<INTEGER> p(std::min(m,n)+3,-71),q=p;
    for(long j=0;j<n;++j) for(long i=0;i<lda;++i) {
        auto z=a[i+j*lda].mpc_data();
        if(i>=m) { a[i+j*lda]=COMPLEX(731,-17); continue; }
        mpfr_urandomb(mpc_realref(z),rng); mpfr_urandomb(mpc_imagref(z),rng);
        a[i+j*lda]=(a[i+j*lda]-COMPLEX(.5,.5))/COMPLEX(8);
        if(i==j) a[i+j*lda]+=COMPLEX(n+2);
        if(i!=j && (i+j)%7==0) {
            mpfr_set_zero(mpc_realref(z),(i+j)%2?1:-1);
            mpfr_set_zero(mpc_imagref(z),(i+j)%3?1:-1);
        }
    }
    if(chol) {
        for(long j=0;j<n;++j) {
            mpfr_set_zero(mpc_imagref(a[j+j*lda].mpc_data()),1);
            for(long i=j+1;i<m;++i) a[j+i*lda]=conj(a[i+j*lda]);
        }
        if(kind==1 || kind==2 || kind==4) a[(kind==1?0:n/2)*(lda+1)]=-1;
        // A Hermitian factorization ignores the input imaginary diagonal;
        // even a failed pivot must retain exactly what the oracle retains.
        if(kind>=3) for(long j=0;j<n;++j)
            mpfr_set_si(mpc_imagref(a[j+j*lda].mpc_data()),-3,MPFR_RNDN);
    } else {
        if(kind==1) for(long i=0;i<m;++i) a[i+(n/2)*lda]=0;
        if(kind==2 && n>1) for(long i=0;i<m;++i) a[i+lda]=a[i];
        if(kind==3) for(long j=0;j<n;++j) for(long i=0;i<m;++i) a[i+j*lda]=COMPLEX(i+j+1,1);
        if(kind==4) for(long i=0;i<m;++i) a[i]=COMPLEX(i%2?-1:1,i%3?-1:1);
        if(kind==5) for(long i=0;i<m;++i) {
            a[i]=COMPLEX(REAL(1)/(i+1),REAL(.5)/(i+1));
            auto z=a[i].mpc_data();
            mpfr_mul_2si(mpc_realref(z),mpc_realref(z),-1073741823L,MPFR_RNDN);
            mpfr_mul_2si(mpc_imagref(z),mpc_imagref(z),-1073741823L,MPFR_RNDN);
        }
    }
    b=a; auto original=a; INTEGER info=-99,ref=-98; bool ran;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t before=0,after=0; cudaDeviceGetLimit(&before,cudaLimitStackSize);
#endif
    if(chol) { complex_oracle_potrf(upper?"U":"L",n,b.data(),lda,ref); ran=Cpotrf_mpfr_cuda(upper?"U":"L",n,a.data(),lda,info); }
    else { complex_oracle_getrf(m,n,b.data(),lda,q.data(),ref); ran=Cgetrf_mpfr_cuda(m,n,a.data(),lda,p.data(),info); }
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    cudaDeviceGetLimit(&after,cudaLimitStackSize); check(before==after,"stack restored",m,n,kind);
#endif
    check(ran && info==ref && p==q,chol?"potrf status":"getrf status",m,n,kind);
    long bad=0;
    for(size_t i=0;i<a.size();++i) if(!same(a[i],b[i])) {
        if(!bad && errors<10) std::fprintf(stderr,"first differing complex element=%zu\n",i);
        ++bad;
    }
    check(!bad,chol?"potrf A":"getrf A",m,n,kind);
    if(!kind && n<=65 && ran && !info) check(backward(chol,upper,m,n,lda,a,original,p),"backward error",m,n,kind);
}
void rejected() {
    std::vector<COMPLEX> a(12,COMPLEX(3)),b=a; INTEGER p[2]={-7,-8},info=-9;
    auto untouched=[&] { bool ok=info==-9 && p[0]==-7 && p[1]==-8; for(size_t i=0;i<a.size();++i) ok=same(a[i],b[i])&&ok; return ok; };
    check(!Cgetrf_mpfr_cuda(0,2,a.data(),3,p,info)&&untouched(),"quick return");
    check(!Cgetrf_mpfr_cuda(3,0,a.data(),3,p,info)&&untouched(),"zero columns");
    check(!Cpotrf_mpfr_cuda("L",0,a.data(),3,info)&&untouched(),"potrf quick return");
    check(!Cgetrf_mpfr_cuda(2,2,a.data(),3,nullptr,info)&&untouched(),"null pivots");
    check(!Cpotrf_mpfr_cuda(nullptr,2,a.data(),3,info)&&untouched(),"null flag");
    check(!Cgetrf_mpfr_cuda(3,2,a.data(),2,p,info)&&untouched(),"bad lda");
    check(!Cpotrf_mpfr_cuda("?",2,a.data(),3,info)&&untouched(),"bad uplo");
    check(!Cgetrf_mpfr_cuda(std::numeric_limits<INTEGER>::max(),2,a.data(),std::numeric_limits<INTEGER>::max(),p,info)&&untouched(),"size overflow");
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    check(!Cgetrf_mpfr_cuda(1048576,1048576,a.data(),1048576,p,info)&&untouched(),"memory budget");
#endif
    setenv("MPLAPACK_MPFR_LAPACK_CUDA","0",1);
    check(!Cgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"disabled");
    setenv("MPLAPACK_MPFR_LAPACK_CUDA","1",1); setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","3",1);
    check(!Cgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"threshold");
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    mpfr_set_default_rounding_mode(MPFR_RNDZ);
    check(!Cgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"rounding");
    mpfr_set_default_rounding_mode(MPFR_RNDN);
    auto emin=mpfr_get_emin(); mpfr_set_emin(emin+1);
    check(!Cpotrf_mpfr_cuda("L",2,a.data(),3,info)&&untouched(),"exponent range"); mpfr_set_emin(emin);
    mpfr_set_prec(mpc_imagref(a[0].mpc_data()),mpfr_get_default_prec()+1);
    mpfr_set_prec(mpc_imagref(b[0].mpc_data()),mpfr_get_default_prec()+1); a[0]=3; b=a;
    check(!Cgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"mixed precision");
    mpfr_set_prec(mpc_imagref(a[0].mpc_data()),mpfr_get_default_prec());
    mpfr_set_prec(mpc_imagref(b[0].mpc_data()),mpfr_get_default_prec());
    a.assign(12,COMPLEX(3)); mpfr_set_inf(mpc_realref(a[0].mpc_data()),1); b=a;
    check(!Cgetrf_mpfr_cuda(2,2,a.data(),3,p,info)&&untouched(),"nonfinite");
}
}
int main(int argc,char **argv) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
    bool quick=argc>2 && !std::strcmp(argv[2],"--quick");
    bool frames=argc>2 && !std::strcmp(argv[2],"--stack-frames");
    bool probe=frames || (argc>2 && !std::strcmp(argv[2],"--stack-probe"));
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    gmp_randinit_default(rng); gmp_randseed_ui(rng,7); bool selected=false;
    for(long prec:{256,512,768,1024,2048,200}) {
        if(argc>1 && std::atol(argv[1])!=prec) continue; selected=true;
        mpfrxx::set_default_precision_bits(prec); mpfrxx::set_default_mpc_precision_bits(prec);
        long first=checks,ec=elements,zc=zeros,fail=errors,lc=cases[0],cc=cases[1];
        if(probe) {
            // Stack sizing only: every panel branch, blocked updates and both
            // streams. The full/quick acceptance sweeps remain unchanged.
            long n=frames?3:131;
            trial(false,false,4,1,0); trial(false,false,1,4,0);
            for(int k:{0,4,5}) trial(false,false,n+3,n,k);
            for(bool u:{false,true}) for(int k=0;k<5;++k) trial(true,u,n,n,k);
        } else {
        for(long n:{1,2,3,31,32,33,63,64,65,127,128,129,131}) {
            if(quick && n>65) continue;
            for(int kind=0;kind<6;++kind) { trial(false,false,n,n,kind); trial(false,false,n+3,n,kind); trial(false,false,n,n+3,kind); }
            for(bool u:{false,true}) for(int k=0;k<5;++k) trial(true,u,n,n,k);
        }
        for(long n:{257,599}) {
            if(quick) continue;
            trial(false,false,n,n,0); trial(false,false,n+3,n,0); trial(false,false,n,n+3,0);
            trial(true,false,n,n,0); trial(true,true,n,n,0);
        }
        if(quick) {
            trial(false,false,134,131,0); trial(false,false,131,134,0);
            trial(true,false,131,131,0); trial(true,true,131,131,0);
        }
        }
        setenv("MPLAPACK_MPFR_LAPACK_FAST","1",1);
        long fast_n=frames?3:33;
        trial(false,false,fast_n+2,fast_n,0); trial(true,false,fast_n,fast_n,0); trial(true,true,fast_n,fast_n,0);
        unsetenv("MPLAPACK_MPFR_LAPACK_FAST"); rejected();
        std::printf("p=%ld LU=%ld Cholesky=%ld checks=%ld components=%ld errors=%ld signed_zero_diffs=%ld\n",prec,cases[0]-lc,cases[1]-cc,checks-first,elements-ec,errors-fail,zeros-zc); std::fflush(stdout);
    }
    gmp_randclear(rng); return !selected?2:errors?1:0;
}
