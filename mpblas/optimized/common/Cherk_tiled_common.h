// Ordered HERK NoTrans tiles. Same load/step phases in CUDA and host emulation.
#ifndef MPLAPACK_CHERK_TILED_COMMON_H
#define MPLAPACK_CHERK_TILED_COMMON_H
namespace mplapack_blas_dev {
template<class Job> struct cgpu_is_herk { static constexpr bool value=false; };
template<class Ops> struct cgpu_is_herk<job_Cherk<Ops>> {
    static constexpr bool value=sizeof(typename job_Cherk<Ops>::DC)<=64;
};
template<class Job> struct herk_shared {
    using V=typename Job::DC;
    static constexpr int T=sizeof(V)<=64?16:sizeof(V)<=256?8:4;
    V a[T][T+1], temp[T][T+1];
    bool active[T][T+1];
};
template<class Job> struct herk_constants {
#if defined(__CUDACC__)
    static __forceinline__ bool zero_alpha(const Job &j) { return j.alpha==typename Job::DR(0); }
#else
    static inline bool zero_alpha(const Job &j) { return j.alpha==typename Job::DR(0); }
#endif
};
template<class Job> static inline bool herk_tiled_eligible(const Job &j) {
    return lsame(j.tr,"N") && j.k>0 && !herk_constants<Job>::zero_alpha(j);
}
template<class Job> MPLAPACK_DEV_HD bool herk_mine(const Job &j,long i,long c) {
    return i<j.n && c<j.n && (lsame(j.uplo,"U")?i<=c:i>=c);
}
template<class Job> MPLAPACK_DEV_HD typename Job::DC herk_init(const Job &j,long i,long c) {
    using R=typename Job::DR; using V=typename Job::DC;
    if(!herk_mine(j,i,c) || j.beta==R(0)) return V(0);
    const V &old=j.C[i+c*j.ldc];
    if(i==c) return j.beta==R(1)?V(old.real()):V(j.beta*old.real());
    return j.beta==R(1)?old:j.beta*old;
}
template<class Job> MPLAPACK_DEV_HD void herk_load(const Job &j,herk_shared<Job> &s,long i0,long j0,long p,int x,int y) {
    if(p+y>=j.k) return;
    if(i0+x<j.n) s.a[y][x]=j.A[i0+x+(p+y)*j.lda];
    if(j0+x<j.n) {
        const auto v=j.A[j0+x+(p+y)*j.lda];
        s.active[y][x]=!(v==typename Job::DC(0));
        // The skip belongs to the unscaled value, not the rounded temp.
        if(s.active[y][x]) s.temp[y][x]=j.alpha*conj(v);
    }
}
template<class Job> MPLAPACK_DEV_HD void herk_step(const herk_shared<Job> &s,typename Job::DC &v,long count,int x,int y,bool diagonal) {
    for(long k=0;k<count;++k) if(s.active[k][y]) {
        const auto product=s.temp[k][y]*s.a[k][x];
        if(diagonal) v=typename Job::DC(v.real()+product.real());
        else v+=product;
    }
}
#if MPLAPACK_CBLAS_CUDA_DEVICE
template<class Job> __device__ void herk_tiled_body(const Job &j) {
    constexpr int T=herk_shared<Job>::T;
    const int x=threadIdx.x,y=threadIdx.y;
    const long i0=(long)blockIdx.x*T,j0=(long)blockIdx.y*T;
    if(lsame(j.uplo,"U")?i0>j0:j0>i0) return;
    __shared__ __align__(16) unsigned char bytes[sizeof(herk_shared<Job>)];
    auto &s=*reinterpret_cast<herk_shared<Job>*>(bytes);
    const long i=i0+x,c=j0+y;
    auto value=herk_init(j,i,c);
    const bool mine=herk_mine(j,i,c);
    for(long p=0;p<j.k;p+=T) {
        herk_load(j,s,i0,j0,p,x,y); __syncthreads();
        if(mine) herk_step(s,value,j.k-p<T?j.k-p:T,x,y,i==c);
        __syncthreads();
    }
    if(mine) j.C[i+c*j.ldc]=value;
}
template<class Job> __global__ void k_herk_tiled(const Job *arg) { herk_tiled_body(*arg); }
template<class Job> __global__ void k_herk_tiled_value(Job arg) { herk_tiled_body(arg); }
template<class Job> static int herk_tiled_stream(const Job &j,cudaStream_t stream) {
    constexpr int T=herk_shared<Job>::T;
    const unsigned blocks=(unsigned)((j.n+T-1)/T);
    if(blocks) k_herk_tiled_value<<<dim3(blocks,blocks),dim3(T,T),0,stream>>>(j);
    return cudaGetLastError();
}
template<class Job> static int herk_tiled_device(const Job &j,const Job *scratch) {
    constexpr int T=herk_shared<Job>::T;
    const unsigned blocks=(unsigned)((j.n+T-1)/T);
    if(blocks) k_herk_tiled<<<dim3(blocks,blocks),dim3(T,T)>>>(scratch);
    return cudaGetLastError();
}
#else
template<class Job> static int herk_tiled_host(const Job &j) {
    constexpr int T=herk_shared<Job>::T;
    for(long j0=0;j0<j.n;j0+=T) for(long i0=0;i0<j.n;i0+=T) {
        if(lsame(j.uplo,"U")?i0>j0:j0>i0) continue;
        herk_shared<Job> s;
        typename Job::DC value[T][T];
        for(int y=0;y<T;++y) for(int x=0;x<T;++x) value[y][x]=herk_init(j,i0+x,j0+y);
        for(long p=0;p<j.k;p+=T) {
            for(int y=0;y<T;++y) for(int x=0;x<T;++x) herk_load(j,s,i0,j0,p,x,y);
            for(int y=0;y<T;++y) for(int x=0;x<T;++x)
                if(herk_mine(j,i0+x,j0+y)) herk_step(s,value[y][x],j.k-p<T?j.k-p:T,x,y,i0+x==j0+y);
        }
        for(int y=0;y<T;++y) for(int x=0;x<T;++x)
            if(herk_mine(j,i0+x,j0+y)) j.C[i0+x+(j0+y)*j.ldc]=value[y][x];
    }
    return 0;
}
#endif
}
#endif
