// Canonical finite-zero contract, with a frozen abd9407 nonzero oracle.
// Standalone: c++ -O2 -std=c++17 -I<mpc_cuda>/include cu_fzero_test.cpp -lmpfr -lgmp
#include <mpc_cuda/cu_freal.cuh>
#include "cu_freal_before_zero.cuh"
#include <mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define CUDA(call) do { auto e=(call); if(e!=cudaSuccess){std::fprintf(stderr,"%s: %s\n",#call,cudaGetErrorString(e));std::exit(1);} }while(0)
#endif
using cu_fp::cu_freal;
template<class A,class B> __host__ __device__ bool equal(const A&a,const B&b){
    if(a.sign!=b.sign||a.exp!=b.exp)return false;
    for(int i=0;i<A::N;++i)if(a.m[i]!=b.m[i])return false;
    return true;
}
template<int P> __host__ __device__ cu_freal<P> apply(cu_freal<P>a,cu_freal<P>b,int op){
    switch(op){
    case 0:return a+b; case 1:return a-b; case 2:return a*b;
    case 3:return cu_fp::cu_neg(a); case 4:return cu_fp::cu_abs(a);
    case 5:return cu_fp::cu_scale2(a,3);
    case 6:return cu_fp::cu_convert<P>(cu_fp::cu_convert<P+32>(a));
    default:return cu_freal<P>::from_double(a.to_double());
    }
}
template<int P> __host__ __device__ bool old_equal(cu_freal<P>a,cu_freal<P>b,int op,const cu_freal<P>&r){
    cu_fp_before_zero::cu_freal<P> x,y;
    x.sign=a.sign;x.exp=a.exp;y.sign=b.sign;y.exp=b.exp;
    for(int i=0;i<x.N;++i){x.m[i]=a.m[i];y.m[i]=b.m[i];}
    switch(op){
    case 0:return equal(r,x+y);case 1:return equal(r,x-y);case 2:return equal(r,x*y);
    case 3:return equal(r,cu_fp_before_zero::cu_neg(x));
    case 4:return equal(r,cu_fp_before_zero::cu_abs(x));
    case 5:return equal(r,cu_fp_before_zero::cu_scale2(x,3));
    case 6:return equal(r,cu_fp_before_zero::cu_convert<P>(cu_fp_before_zero::cu_convert<P+32>(x)));
    default:return equal(r,cu_fp_before_zero::cu_freal<P>::from_double(x.to_double()));
    }
}
template<int P> __host__ __device__ void evaluate(const cu_freal<P>*a,const cu_freal<P>*b,cu_freal<P>*out,int*oldbad,int i){
    for(int op=0;op<8;++op){auto r=apply(a[i],b[i],op);out[8*i+op]=r;
        oldbad[8*i+op]=!r.is_zero()&&!old_equal(a[i],b[i],op,r);
    }
}
#ifdef __CUDACC__
template<int P> __global__ void kernel(const cu_freal<P>*a,const cu_freal<P>*b,cu_freal<P>*out,int*bad,int n){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)evaluate(a,b,out,bad,i);
}
#endif
template<int P> void to_mpfr(mpfr_ptr x,const cu_freal<P>&a){
    if(a.is_zero()){mpfr_set_zero(x,a.sign);return;}
    mp_limb_t limbs[cu_freal<P>::N];for(int i=0;i<a.N;++i)limbs[i]=a.m[i];
    mpfr_t view;mpfr_custom_init_set(view,a.sign*MPFR_REGULAR_KIND,a.exp,P,limbs);mpfr_set(x,view,MPFR_RNDN);
}
template<int P> bool matches(const cu_freal<P>&r,mpfr_srcptr x){
    if(mpfr_zero_p(x)){cu_freal<P> z;z.set_zero();return equal(r,z);}
    if(r.is_zero()||r.sign!=(mpfr_signbit(x)?-1:1)||r.exp!=mpfr_get_exp(x))return false;
    auto p=static_cast<const mp_limb_t*>(mpfr_custom_get_significand(x));
    for(int i=0;i<r.N;++i)if(r.m[i]!=p[i])return false;return true;
}
static unsigned long long rng(unsigned long long&s){s^=s<<13;s^=s>>7;s^=s<<17;return s;}
template<int P> int run(){
    using F=cu_freal<P>;std::vector<F> values;
    for(double d:{0.0,-0.0,1.0,-1.0,0.5,-0.5,2.0,-2.0})values.emplace_back(d);
    values[1].sign=-1;
    // Poisoned zero limbs are ignored whenever exp is the zero sentinel.
    F poison=values[1];for(auto &x:poison.m)x=~0ULL;values.push_back(poison);
    unsigned long long seed=0x123456789abcdefULL;
    for(int k=0;k<64;++k){F x;x.sign=(k&1)?-1:1;x.exp=k%21-10;
        for(auto &v:x.m)v=rng(seed);x.m[F::N-1]|=1ULL<<63;x.m[0]&=~((1ULL<<F::SB)-1);values.push_back(x);}
    std::vector<F>a,b;
    for(const auto&x:values)for(const auto&y:values){a.push_back(x);b.push_back(y);}
    int n=a.size();std::vector<F>out(8*n);std::vector<int>bad(8*n);
#ifdef __CUDACC__
    F *da,*db,*dr;int*de;CUDA(cudaMalloc(&da,n*sizeof(F)));CUDA(cudaMalloc(&db,n*sizeof(F)));
    CUDA(cudaMalloc(&dr,out.size()*sizeof(F)));CUDA(cudaMalloc(&de,bad.size()*sizeof(int)));
    CUDA(cudaMemcpy(da,a.data(),n*sizeof(F),cudaMemcpyHostToDevice));CUDA(cudaMemcpy(db,b.data(),n*sizeof(F),cudaMemcpyHostToDevice));
    kernel<P><<<(n+63)/64,64>>>(da,db,dr,de,n);CUDA(cudaGetLastError());
    CUDA(cudaMemcpy(out.data(),dr,out.size()*sizeof(F),cudaMemcpyDeviceToHost));CUDA(cudaMemcpy(bad.data(),de,bad.size()*sizeof(int),cudaMemcpyDeviceToHost));
    CUDA(cudaFree(da));CUDA(cudaFree(db));CUDA(cudaFree(dr));CUDA(cudaFree(de));
#else
    for(int i=0;i<n;++i)evaluate(a.data(),b.data(),out.data(),bad.data(),i);
#endif
    mpfr_t x,y,z;mpfr_inits2(P,x,y,z,(mpfr_ptr)0);int errors=0,nonzero_old=0;
    for(int i=0;i<n;++i){to_mpfr(x,a[i]);to_mpfr(y,b[i]);
        for(int op=0;op<8;++op){switch(op){
        case 0:mpfr_add(z,x,y,MPFR_RNDN);break;case 1:mpfr_sub(z,x,y,MPFR_RNDN);break;case 2:mpfr_mul(z,x,y,MPFR_RNDN);break;
        case 3:mpfr_neg(z,x,MPFR_RNDN);break;case 4:mpfr_abs(z,x,MPFR_RNDN);break;case 5:mpfr_mul_2si(z,x,3,MPFR_RNDN);break;
        case 6:mpfr_set(z,x,MPFR_RNDN);break;default:mpfr_set_d(z,mpfr_get_d(x,MPFR_RNDN),MPFR_RNDN);break;}
        const auto&r=out[8*i+op];nonzero_old+=!r.is_zero();
        if(!matches(r,z)||bad[8*i+op]){if(errors<8)std::fprintf(stderr,"PB=%d pair=%d op=%d oldbad=%d\n",P,i,op,bad[8*i+op]);++errors;}
    }}
    mpfr_clears(x,y,z,(mpfr_ptr)0);
    std::printf("PB=%d pairs=%d MPFR_checks=%d nonzero_old_checks=%d errors=%d\n",P,n,8*n,nonzero_old,errors);return errors;
}
int main(){
#ifdef __CUDACC__
    int n=0;if(cudaGetDeviceCount(&n)!=cudaSuccess||n==0)return 77;
#endif
    return (run<224>()+run<256>()+run<512>()+run<768>()+run<992>()+run<1024>()+run<2048>())?1:0;
}
