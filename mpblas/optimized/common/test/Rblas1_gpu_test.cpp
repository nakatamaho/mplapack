// Focused standalone GPU entry tests, also usable without a CUDA device.
#define main blas1_cpu_checks
#include "Rblas1_fast_test.cpp"
#undef main
#include "../Rblas_cuda_entry_common.h"
#if defined(MPLAPACK_BUILD_WITH_MPFR)
#define PREFIX "MPLAPACK_MPFR"
#else
#if defined(MPLAPACK_BUILD_WITH_DD)
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define PREFIX "MPLAPACK_TD"
#else
#define PREFIX "MPLAPACK_QD"
#endif
#endif
#if !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
#include <cuda_runtime.h>
#endif
#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
#if defined(MPLAPACK_BUILD_WITH_DD)
#include "../../dd/cuda/Rgemm_gpu_cuda.cu"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include "../../td/cuda/Rgemm_gpu_cuda.cu"
#else
#include "../../qd/cuda/Rgemm_gpu_cuda.cu"
#endif
#endif
int main() {
#if !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
    int devices=0;if(cudaGetDeviceCount(&devices) || !devices) return 77;
#endif
    setenv(PREFIX "_BLAS_CUDA_MIN","0",1);setenv(PREFIX "_BLAS_CUDA_L12","1",1);
    setenv(PREFIX "_CUDA_MIN_MNK","0",1);setenv(PREFIX "_CUDA_L12","1",1);
    const REAL u=REAL(4.0)*mplapack::get_arithmetic_params<REAL>().eps;
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    if(mpfrxx::default_precision_bits()>4096) {
        setenv(PREFIX "_MPBLAS_FAST","1",1);
        for(long n:{1L,7L,257L}) {
            std::vector<REAL> x(n),y(n);
            for(long i=0;i<n;++i) {x[i]=REAL(i%13-6)/REAL(7.0);y[i]=REAL(i%11-5)/REAL(9.0);}
            REAL r,again,tol=REAL(32*n*n)*u;INTEGER index=0;
            check(MPLAPACK_CUDA_ENTRY(Rdot)(n,x.data(),-1,y.data(),1,r),"dynamic dot status");
            check(abs(r-test_Rdot_reference(n,x.data(),-1,y.data(),1))<=tol,"dynamic dot bound");
            check(MPLAPACK_CUDA_ENTRY(Rdot)(n,x.data(),-1,y.data(),1,again) && same(r,again),"dynamic repeat");
            check(MPLAPACK_CUDA_ENTRY(Rasum)(n,x.data(),1,r) && abs(r-test_Rasum_reference(n,x.data(),1))<=tol,"dynamic asum");
            REAL norm_ref=test_Rnrm2_reference(n,x.data(),-1);
            REAL norm_tol=REAL(16+8*std::ceil(std::log2(double(std::max(2L,2*n)))))*u*norm_ref;
            check(MPLAPACK_CUDA_ENTRY(Rnrm2)(n,x.data(),-1,r) && abs(r-norm_ref)<=norm_tol,"dynamic norm");
            check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),1,index) && index==test_iRamax_reference(n,x.data(),1),"dynamic iamax");
            auto expected=y;
            for(long i=0;i<n;++i) expected[i]+=REAL(0.7)*x[n-1-i];
            check(MPLAPACK_CUDA_ENTRY(Raxpy)(n,REAL(0.7),x.data(),-1,y.data(),1),"dynamic axpy status");
            for(long i=0;i<n;++i) check(same(y[i],expected[i]),"dynamic axpy identity");
            check(MPLAPACK_CUDA_ENTRY(Rcopy)(n,x.data(),1,y.data(),1),"dynamic copy status");
            for(long i=0;i<n;++i) check(same(y[i],x[i]),"dynamic copy identity");
        }
        std::printf("BLAS1 dynamic real: %ld checks, %ld failures\n",checks,failures);
        return failures?1:0;
    }
#endif
    for(long n:{0L,1L,7L,64L,257L,1025L}) for(long inc:{-3L,1L,2L}) {
        size_t size=std::max(1L,1+(n-1)*std::labs(inc));
        std::vector<REAL> x(size),y(size);std::vector<COMPLEX> z(size),w(size);
        for(size_t i=0;i<size;++i) {
            x[i]=REAL(double(int(i%17)-8))/REAL(13.0);y[i]=REAL(double(int(i%11)-5))/REAL(7.0);
            z[i]=COMPLEX(x[i],y[i]);w[i]=COMPLEX(y[i],x[i]);
        }
        REAL r=17.0,again=19.0;COMPLEX c(17.0,19.0),cagain;INTEGER index=77;
        unsetenv(PREFIX "_MPBLAS_FAST");
        check(!MPLAPACK_CUDA_ENTRY(Rdot)(n,x.data(),inc,y.data(),inc,r) && r==REAL(17.0),"default false/result untouched");
        check(!MPLAPACK_CUDA_ENTRY(Cdotc)(n,z.data(),inc,w.data(),inc,c) && same(c,COMPLEX(17.0,19.0)),"Cdotc default untouched");
        check(!MPLAPACK_CUDA_ENTRY(Cdotu)(n,z.data(),inc,w.data(),inc,c) && same(c,COMPLEX(17.0,19.0)),"Cdotu default untouched");
        check(!MPLAPACK_CUDA_ENTRY(Rasum)(n,x.data(),inc,r) && r==REAL(17.0),"Rasum default untouched");
        check(!MPLAPACK_CUDA_ENTRY(RCasum)(n,z.data(),inc,r) && r==REAL(17.0),"RCasum default untouched");
        check(!MPLAPACK_CUDA_ENTRY(Rnrm2)(n,x.data(),inc,r) && r==REAL(17.0),"Rnrm2 default untouched");
        check(!MPLAPACK_CUDA_ENTRY(RCnrm2)(n,z.data(),inc,r) && r==REAL(17.0),"RCnrm2 default untouched");
        setenv(PREFIX "_MPBLAS_FAST","1",1);
        REAL tol=REAL(32.0)*REAL(std::max(1L,n))*REAL(std::max(1L,n))*u;
#define REAL_DOT_CHECK(NAME,REF) \
        check(MPLAPACK_CUDA_ENTRY(NAME)(n,x.data(),inc,y.data(),inc,r)==(n>0),#NAME " status"); \
        if(n>0) { check(MPLAPACK_CUDA_ENTRY(NAME)(n,x.data(),inc,y.data(),inc,again) && same(r,again),#NAME " repeat"); \
            check(abs(r-REF(n,x.data(),inc,y.data(),inc))<=tol,#NAME " bound"); }
        REAL_DOT_CHECK(Rdot,test_Rdot_reference)
#undef REAL_DOT_CHECK
#define COMPLEX_DOT_CHECK(NAME,REF) \
        check(MPLAPACK_CUDA_ENTRY(NAME)(n,z.data(),inc,w.data(),inc,c)==(n>0),#NAME " status"); \
        if(n>0) { check(MPLAPACK_CUDA_ENTRY(NAME)(n,z.data(),inc,w.data(),inc,cagain) && same(c,cagain),#NAME " repeat"); \
            check(abs(c-REF(n,z.data(),inc,w.data(),inc))<=tol,#NAME " bound"); }
        COMPLEX_DOT_CHECK(Cdotc,test_Cdotc_reference)
        COMPLEX_DOT_CHECK(Cdotu,test_Cdotu_reference)
#undef COMPLEX_DOT_CHECK
        REAL norm_units=REAL(16+8*std::ceil(std::log2(double(std::max(2L,2*n)))))*u;
#define REDUCE_CHECK(NAME,REF,X,OK,BOUND) \
        check(MPLAPACK_CUDA_ENTRY(NAME)(n,X.data(),inc,r)==(OK),#NAME " status"); \
        if(OK) { check(MPLAPACK_CUDA_ENTRY(NAME)(n,X.data(),inc,again) && same(r,again),#NAME " repeat"); \
            check(abs(r-REF(n,X.data(),inc))<=BOUND,#NAME " bound"); }
        REDUCE_CHECK(Rasum,test_Rasum_reference,x,n>0 && inc>0,tol)
        REDUCE_CHECK(RCasum,test_RCasum_reference,z,n>0 && inc>0,tol)
        REDUCE_CHECK(Rnrm2,test_Rnrm2_reference,x,n>0,norm_units*test_Rnrm2_reference(n,x.data(),inc))
        REDUCE_CHECK(RCnrm2,test_RCnrm2_reference,z,n>0,norm_units*test_RCnrm2_reference(n,z.data(),inc))
#undef REDUCE_CHECK
        check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),inc,index)==(n>0 && inc>0),"iRamax status");
        if(n>0 && inc>0) check(index==test_iRamax_reference(n,x.data(),inc),"iRamax first tie");
        check(MPLAPACK_CUDA_ENTRY(iCamax)(n,z.data(),inc,index)==(n>0 && inc>0),"iCamax status");
        if(n>0 && inc>0) check(index==test_iCamax_reference(n,z.data(),inc),"iCamax first tie");
        auto expected=y;long ix=inc<0?-(n-1)*inc:0;
        for(long i=0;i<n;++i) expected[ix+i*inc]+=REAL(0.7)*x[ix+i*inc];
        check(MPLAPACK_CUDA_ENTRY(Raxpy)(n,REAL(0.7),x.data(),inc,y.data(),inc)==(n>0),"Raxpy status");
        for(size_t i=0;i<size;++i) check(same(y[i],expected[i]),"Raxpy padding/identity");
        for(long i=0;i<n;++i) expected[ix+i*inc]=x[ix+i*inc];
        check(MPLAPACK_CUDA_ENTRY(Rcopy)(n,x.data(),inc,y.data(),inc)==(n>0),"Rcopy status");
        for(size_t i=0;i<size;++i) check(same(y[i],expected[i]),"Rcopy padding/identity");
        r=17.0;check(!MPLAPACK_CUDA_ENTRY(Rdot)(n,x.data(),0,y.data(),inc,r) && r==REAL(17.0),"invalid untouched");
    }
    for(double magnitude:{0.0,1e300,1e-300}) {
        REAL x[]={REAL(magnitude),REAL(magnitude),REAL(0.0)},r=17.0,again;
        COMPLEX z[]={COMPLEX(x[0],x[1]),COMPLEX(x[2],x[0])};
        REAL floor=0.0;
#if !defined(MPLAPACK_BUILD_WITH_MPFR)
        floor=REAL(0x1p-1074)*REAL(32.0);
#endif
        REAL ref=test_Rnrm2_reference(3,x,1);
        check(MPLAPACK_CUDA_ENTRY(Rnrm2)(3,x,1,r),"scaled Rnrm2 status");
        check(abs(r-ref)<=REAL(64.0)*u*ref+floor,"scaled Rnrm2 bound");
        check(MPLAPACK_CUDA_ENTRY(Rnrm2)(3,x,1,again) && same(r,again),"scaled Rnrm2 repeat");
        ref=test_RCnrm2_reference(2,z,1);
        check(MPLAPACK_CUDA_ENTRY(RCnrm2)(2,z,1,r),"scaled RCnrm2 status");
        check(abs(r-ref)<=REAL(64.0)*u*ref+floor,"scaled RCnrm2 bound");
        check(MPLAPACK_CUDA_ENTRY(RCnrm2)(2,z,1,again) && same(r,again),"scaled RCnrm2 repeat");
    }
    REAL negative_zero=-REAL(0.0),zero_output=REAL(1.0),zero_norm=REAL(1.0);
    check(MPLAPACK_CUDA_ENTRY(Rcopy)(1,&negative_zero,1,&zero_output,1) && same(zero_output,negative_zero),"copy signed zero");
    REAL expected_zero=negative_zero;expected_zero+=REAL(1.0)*negative_zero;
    check(MPLAPACK_CUDA_ENTRY(Raxpy)(1,REAL(1.0),&negative_zero,1,&zero_output,1) && same(zero_output,expected_zero),"axpy signed zero");
    check(MPLAPACK_CUDA_ENTRY(Rnrm2)(1,&negative_zero,1,zero_norm) && same(zero_norm,REAL(0.0)),"norm positive zero");
    REAL ties[]={REAL(0.0),REAL(-7.0),REAL(7.0)};
    COMPLEX cties[]={COMPLEX(0.0),COMPLEX(-3.0,4.0),COMPLEX(7.0)};
    INTEGER index=77;
    check(MPLAPACK_CUDA_ENTRY(iRamax)(3,ties,1,index) && index==2,"explicit real tie");
    check(MPLAPACK_CUDA_ENTRY(iCamax)(3,cties,1,index) && index==2,"explicit complex tie");
    setenv(PREFIX "_BLAS_CUDA_L12","0",1);setenv(PREFIX "_CUDA_L12","0",1);
    REAL r=17.0;
    check(!MPLAPACK_CUDA_ENTRY(Rdot)(3,ties,1,ties,1,r) && r==REAL(17.0),"per-call disabled");
    check(!MPLAPACK_CUDA_ENTRY(Raxpy)(3,REAL(1.0),ties,1,ties,1) && ties[1]==REAL(-7.0),"disabled output untouched");
    // Goal 13: one-pass GPU iamax. Chunk-boundary ties, all-equal vectors,
    // the maximum at the first and the last index (the scores are exact, so the
    // index must equal the reference scan in every launch shape).
    {
        setenv(PREFIX "_BLAS_CUDA_L12","1",1);setenv(PREFIX "_CUDA_L12","1",1);
        const long C=4096;
        std::mt19937_64 r3(20261110);
        for(long n:{C-1,C,C+1,2*C-1,2*C,2*C+1,3*C}) {
            std::vector<REAL> x((size_t)n);
            for(auto &v:x) v=REAL(double(int(r3()%2001)-1000))/REAL(997.0);
            INTEGER index=77;
            check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),1,index) && index==test_iRamax_reference(n,x.data(),1),"GPU iamax chunk random");
            std::vector<COMPLEX> z((size_t)n);
            for(long i=0;i<n;++i) z[i]=COMPLEX(x[i],x[(i+7)%n]);
            index=77;
            check(MPLAPACK_CUDA_ENTRY(iCamax)(n,z.data(),1,index) && index==test_iCamax_reference(n,z.data(),1),"GPU iCamax chunk random");
        }
        {   // all equal and first/last maxima
            const long n=2*C+3;
            std::vector<REAL> x((size_t)n,REAL(1.5));
            std::vector<COMPLEX> z((size_t)n,COMPLEX(REAL(0.5),REAL(0.5)));
            INTEGER index=77;
            check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),1,index) && index==1,"GPU iamax all-equal");
            check(MPLAPACK_CUDA_ENTRY(iCamax)(n,z.data(),1,index) && index==1,"GPU iCamax all-equal");
            x[0]=REAL(8.0);index=77;
            check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),1,index) && index==1,"GPU iamax max first");
            x[0]=REAL(1.5);x[n-1]=REAL(8.0);index=77;
            check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),1,index) && index==n,"GPU iamax max last");
            x[n-1]=REAL(1.5);x[C-1]=REAL(6.0);x[C]=REAL(6.0);index=77;
            check(MPLAPACK_CUDA_ENTRY(iRamax)(n,x.data(),1,index) && index==C,"GPU iamax boundary tie");
            z[C-1]=COMPLEX(REAL(2.0),REAL(2.0));z[C]=COMPLEX(REAL(1.0),REAL(3.0));index=77;
            check(MPLAPACK_CUDA_ENTRY(iCamax)(n,z.data(),1,index) && index==C,"GPU iCamax boundary tie");
        }
    }
    std::printf("BLAS1 GPU: %ld checks, %ld failures\n",checks,failures);
    unsetenv(PREFIX "_MPBLAS_FAST");
    return blas1_cpu_checks();
}
