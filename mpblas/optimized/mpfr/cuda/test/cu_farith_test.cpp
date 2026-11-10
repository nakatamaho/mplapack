// Goal12 three-way oracle: candidate, frozen canonical-zero baseline, MPFR RNDN.
#include <mpc_cuda/cu_freal.cuh>
#include "cu_freal_before_fast.cuh"
#include <mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define CUDA(x) do{auto e=(x);if(e!=cudaSuccess){fprintf(stderr,"%s: %s\n",#x,cudaGetErrorString(e));exit(1);}}while(0)
#endif
template<class A,class B> __host__ __device__ bool same(const A&a,const B&b){
    if(a.sign!=b.sign||a.exp!=b.exp)return false;
    for(int i=0;i<A::N;++i)if(a.m[i]!=b.m[i])return false;return true;
}
template<int P> using F=cu_fp::cu_freal<P>;
template<int P> __host__ __device__ void evaluate(const F<P>*a,const F<P>*b,F<P>*out,int*bad,int i){
    auto x=a[i],y=b[i];cu_fp_before_fast::cu_freal<P>u,v;
    u.sign=x.sign;u.exp=x.exp;v.sign=y.sign;v.exp=y.exp;
    for(int j=0;j<u.N;++j){u.m[j]=x.m[j];v.m[j]=y.m[j];}
    auto r0=x+y,r1=x-y,r2=x*y;out[3*i]=r0;out[3*i+1]=r1;out[3*i+2]=r2;
    bad[i]=!same(r0,u+v)||!same(r1,u-v)||!same(r2,u*v);
}
#ifdef __CUDACC__
template<int P> __global__ void kernel(const F<P>*a,const F<P>*b,F<P>*o,int*bad,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)evaluate(a,b,o,bad,i);}
#endif
template<int P> void to_mpfr(mpfr_ptr x,const F<P>&a){
    if(a.is_zero()){mpfr_set_zero(x,a.sign);return;}
    mp_limb_t m[F<P>::N];for(int i=0;i<F<P>::N;++i)m[i]=a.m[i];
    mpfr_t v;mpfr_custom_init_set(v,a.sign*MPFR_REGULAR_KIND,a.exp,P,m);mpfr_set(x,v,MPFR_RNDN);
}
template<int P> F<P> from_mpfr(mpfr_srcptr x){
    F<P> a;if(mpfr_zero_p(x)){a.set_zero();return a;}a.sign=mpfr_signbit(x)?-1:1;a.exp=mpfr_get_exp(x);
    auto m=static_cast<const mp_limb_t*>(mpfr_custom_get_significand(x));for(int i=0;i<a.N;++i)a.m[i]=m[i];return a;
}
static unsigned long long random64(unsigned long long&s){s^=s<<13;s^=s>>7;s^=s<<17;return s;}
template<int P> F<P> random_value(unsigned long long&s,long e){
    F<P>a;a.sign=(random64(s)&1)?1:-1;a.exp=e;for(auto&v:a.m)v=random64(s);
    a.m[a.N-1]|=1ULL<<63;a.m[0]&=~((1ULL<<a.SB)-1);return a;
}
template<int P> int run(long count){
    constexpr int batch=16384;std::vector<F<P>>a(batch),b(batch),out(3*batch);std::vector<int>bad(batch);
    std::vector<F<P>>hard_a,hard_b;mpfr_t x,y,z,t;mpfr_inits2(P,x,y,z,t,(mpfr_ptr)0);
    unsigned long long seed=0xf1337babaULL;
    auto push=[&](){hard_a.push_back(from_mpfr<P>(x));hard_b.push_back(from_mpfr<P>(y));};
    for(int k=1;k<=P;++k){mpfr_set_ui(x,1,MPFR_RNDN);mpfr_set_ui_2exp(t,1,-k,MPFR_RNDN);mpfr_sub(y,x,t,MPFR_RNDN);push();mpfr_neg(y,y,MPFR_RNDN);push();}
    for(int sign:{1,-1})for(int shift:{-1000000,0,1000000}){
        mpfr_set_si_2exp(x,sign,shift,MPFR_RNDN);mpfr_set(y,x,MPFR_RNDN);push();mpfr_neg(y,y,MPFR_RNDN);push();
        mpfr_set_zero(y,1);push();mpfr_set_zero(x,1);push();
    }
    // Products exactly at and immediately beside a midpoint: q+r=P.
    for(int q:{1,P/2,P-1})for(int direction:{-1,0,1}){
        mpfr_set_ui(x,1,MPFR_RNDN);mpfr_set_ui_2exp(t,1,-q,MPFR_RNDN);mpfr_add(x,x,t,MPFR_RNDN);
        mpfr_set_ui(y,1,MPFR_RNDN);mpfr_set_ui_2exp(t,1,-(P-q),MPFR_RNDN);mpfr_add(y,y,t,MPFR_RNDN);
        if(direction<0)mpfr_nextbelow(y);if(direction>0)mpfr_nextabove(y);push();
    }
    // d>=1 near-cancellation: A just above 2^e, B just below 2^e.
    for(int k=0;k<200;++k){
        long e=long(random64(seed)%41)-20;
        cu_fp::cu_freal<P> va;va.sign=1;va.exp=e;for(auto&v:va.m)v=0;va.m[va.N-1]=1ULL<<63;
        for(int t=0;t<3&&t<va.N;++t)va.m[va.N-1-t]=random64(seed)|(t?0:(1ULL<<63));
        va.m[0]&=~((1ULL<<va.SB)-1);
        cu_fp::cu_freal<P> vb;vb.sign=(k&1)?-1:1;vb.exp=e;for(auto&v:vb.m)v=0;
        for(auto&v:vb.m)v=~v;vb.m[vb.N-1]=~(1ULL<<(63-k%63));vb.m[0]|=((1ULL<<vb.SB)-1);
        vb.m[0]&=~((1ULL<<vb.SB)-1);vb.m[0]=~((1ULL<<vb.SB)-1);
        hard_a.push_back(va);hard_b.push_back(vb);
    }
    // Carry chains and halfway additions, both parity choices, binade crossing.
    for(int parity=0;parity<2;++parity){mpfr_set_ui(x,1,MPFR_RNDN);mpfr_nextbelow(x);if(parity)mpfr_nextbelow(x);
        mpfr_set_ui_2exp(y,1,-P-1,MPFR_RNDN);push();mpfr_set(y,x,MPFR_RNDN);push();}
    // Matched small exponent gaps (d=0..3) with full-width random significands,
    // and near-equal magnitudes for deep cancellation at d=0.
    for(int g=0;g<4;++g)for(int k=0;k<400;++k){
        long e=long(random64(seed)%41)-20;
        F<P> va=random_value<P>(seed,e),vb=random_value<P>(seed,e-g);
        if(k%2)std::swap(va,vb);
        hard_a.push_back(va);hard_b.push_back(vb);
        if(g==0&&k<200){ // B = A - tiny: deep cancellation, every result length
            F<P> vc=va;unsigned long long borrow=1+(random64(seed)%(unsigned long long)(k+1));
            borrow<<=F<P>::SB; // keep the unused low bits zero
            for(int t=0;t<F<P>::N&&borrow;++t){unsigned long long old=vc.m[t];vc.m[t]-=borrow;borrow=old<borrow;}
            hard_a.push_back(va);hard_b.push_back(vc);
            hard_a.push_back(va);hard_b.push_back(vc);hard_b.back().sign=-vc.sign;
        }
    }
#ifdef __CUDACC__
    F<P>*da,*db,*dr;int*de;CUDA(cudaMalloc(&da,batch*sizeof(F<P>)));CUDA(cudaMalloc(&db,batch*sizeof(F<P>)));
    CUDA(cudaMalloc(&dr,3*batch*sizeof(F<P>)));CUDA(cudaMalloc(&de,batch*sizeof(int)));
#endif
    const int gaps[]={0,1,2,62,63,64,65,127,128,P-2,P-1,P,P+1,P+2,P+129,100000};
    long total=count+hard_a.size(),errors=0;
    for(long start=0;start<total;start+=batch){int n=int(total-start<batch?total-start:batch);
        for(int i=0;i<n;++i){long idx=start+i;
            if(idx<(long)hard_a.size()){a[i]=hard_a[idx];b[i]=hard_b[idx];}
            else{long e=long(random64(seed)%2000001)-1000000;a[i]=random_value<P>(seed,e);b[i]=random_value<P>(seed,e-gaps[random64(seed)%16]);if(random64(seed)&1)std::swap(a[i],b[i]);}}
#ifdef __CUDACC__
        CUDA(cudaMemcpy(da,a.data(),n*sizeof(F<P>),cudaMemcpyHostToDevice));CUDA(cudaMemcpy(db,b.data(),n*sizeof(F<P>),cudaMemcpyHostToDevice));
        kernel<P><<<(n+63)/64,64>>>(da,db,dr,de,n);CUDA(cudaGetLastError());
        CUDA(cudaMemcpy(out.data(),dr,3*n*sizeof(F<P>),cudaMemcpyDeviceToHost));CUDA(cudaMemcpy(bad.data(),de,n*sizeof(int),cudaMemcpyDeviceToHost));
#else
        for(int i=0;i<n;++i)evaluate(a.data(),b.data(),out.data(),bad.data(),i);
#endif
        for(int i=0;i<n;++i){to_mpfr(x,a[i]);to_mpfr(y,b[i]);
            for(int op=0;op<3;++op){if(op==0)mpfr_add(z,x,y,MPFR_RNDN);else if(op==1)mpfr_sub(z,x,y,MPFR_RNDN);else mpfr_mul(z,x,y,MPFR_RNDN);
                if(bad[i]||!same(out[3*i+op],from_mpfr<P>(z))){
                    if(errors==0){fprintf(stderr,"A s=%d e=%ld m:",a[i].sign,(long)a[i].exp);for(int t=0;t<F<P>::N;++t)fprintf(stderr," %016llx",(unsigned long long)a[i].m[t]);
                    fprintf(stderr,"\nB s=%d e=%ld m:",b[i].sign,(long)b[i].exp);for(int t=0;t<F<P>::N;++t)fprintf(stderr," %016llx",(unsigned long long)b[i].m[t]);
                    fprintf(stderr,"\nR s=%d e=%ld m:",out[3*i+op].sign,(long)out[3*i+op].exp);for(int t=0;t<F<P>::N;++t)fprintf(stderr," %016llx",(unsigned long long)out[3*i+op].m[t]);
                    fprintf(stderr,"\nzexp=%ld zero=%d\n",(long)mpfr_get_exp(z),mpfr_zero_p(z));}
                    if(errors++<8)fprintf(stderr,"PB=%d index=%ld op=%d baseline_bad=%d\n",P,start+i,op,bad[i]);}}
        }
        if(errors)break;
    }
#ifdef __CUDACC__
    CUDA(cudaFree(da));CUDA(cudaFree(db));CUDA(cudaFree(dr));CUDA(cudaFree(de));
#endif
    printf("PB=%d random_pairs_per_op=%ld hard_pairs_per_op=%zu errors=%ld\n",P,count,hard_a.size(),errors);fflush(stdout);
    mpfr_clears(x,y,z,t,(mpfr_ptr)0);return errors?1:0;
}
int main(int argc,char**argv){
#ifdef __CUDACC__
    int n=0;if(cudaGetDeviceCount(&n)!=cudaSuccess||!n)return 77;long count=10000000;
#else
    long count=1000000;
#endif
    if(argc>1)count=strtol(argv[1],nullptr,10);if(count<1)return 2;
    return run<224>(count)|run<256>(count)|run<512>(count)|run<768>(count)|run<992>(count)|run<1024>(count)|run<2048>(count);
}
