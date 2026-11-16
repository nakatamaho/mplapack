// Private fixed-complex Upper Cholesky dot tiles; no arena-backed arithmetic.
#ifndef MPLAPACK_MPFR_POTRF_UPPER_CUDA_H
#define MPLAPACK_MPFR_POTRF_UPPER_CUDA_H
namespace mplapack_mpfr_cuda {
namespace {
template<class E> struct upper_update_job {
    using V=typename E::V; using M=typename E::M;
    static constexpr int T=8;
    struct shared { V a[T][T+1],b[T][T+1]; };
    E e; M a,b,c,constants;
    long m,n,k; bool rank;
    L3_HD bool mine(long i,long j) const { return i<m && j<n && (!rank || i<=j); }
    L3_HD void load(shared &sh,long i0,long j0,long q,int x,int y) {
        // k-contiguous loading for both operands of A**H * B.
        if(q+x>=k) return;
        if(i0+y<m) {
            sh.a[x][y]=conj(e.at(a,q+x,i0+y));
        }
        if(j0+y<n) e.set(sh.b[x][y],e.at(b,q+x,j0+y));
    }
    L3_HD void step(const shared &sh,V &v,long i,long j,int count,int x,int y) {
        if(!mine(i,j)) return;
        for(int l=0;l<count;++l) {
            auto product=e.reg(1); e.mul(product,sh.a[l][x],sh.b[l][y]);
            if(rank && i==j) v=V(v.real()+product.real());
            else e.add(v,v,product);
        }
    }
    L3_HD void finish(const V &v,long i,long j) {
        if(!mine(i,j)) return;
        auto out=e.reg(1);
        using R=typename E::DR;
        if(rank) {
            if(i==j) out=R(-1)*v.real()+R(1)*e.at(c,i,j).real();
            else out=R(-1)*v+R(1)*e.at(c,i,j);
        } else out=(-V(1))*v+V(1)*e.at(c,i,j); // reference -cone, including -0
        e.put(c,i,j,out);
    }
};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
template<class J> __global__ void upper_update_kernel(J j,int *failed) {
    __shared__ __align__(16) unsigned char storage[sizeof(typename J::shared)];
    auto &sh=*reinterpret_cast<typename J::shared*>(storage);
    const int x=threadIdx.x,y=threadIdx.y;
    const long i0=(long)blockIdx.x*J::T,j0=(long)blockIdx.y*J::T;
    if(j.rank && i0>j0) return;
    auto value=j.e.reg(0); j.e.zero(value);
    for(long q=0;q<j.k;q+=J::T) {
        j.load(sh,i0,j0,q,x,y); __syncthreads();
        j.step(sh,value,i0+x,j0+y,j.k-q<J::T?j.k-q:J::T,x,y); __syncthreads();
    }
    j.finish(value,i0+x,j0+y);
    if(!j.e.ok) atomicExch(failed,1);
}
template<class J> static int upper_update_device(J j,cudaStream_t stream,int *failed) {
    upper_update_kernel<<<dim3((unsigned)((j.m+J::T-1)/J::T),(unsigned)((j.n+J::T-1)/J::T)),dim3(J::T,J::T),0,stream>>>(j,failed);
    return cudaGetLastError();
}
#else
template<class J> static int upper_update_host(J j) {
    for(long j0=0;j0<j.n;j0+=J::T) for(long i0=0;i0<j.m;i0+=J::T) {
        if(j.rank && i0>j0) continue;
        typename J::shared sh;
        typename J::V values[J::T][J::T];
        for(int y=0;y<J::T;++y) for(int x=0;x<J::T;++x) j.e.zero(values[y][x]);
        for(long q=0;q<j.k;q+=J::T) {
            for(int y=0;y<J::T;++y) for(int x=0;x<J::T;++x) j.load(sh,i0,j0,q,x,y);
            for(int y=0;y<J::T;++y) for(int x=0;x<J::T;++x)
                j.step(sh,values[y][x],i0+x,j0+y,j.k-q<J::T?j.k-q:J::T,x,y);
        }
        for(int y=0;y<J::T;++y) for(int x=0;x<J::T;++x) j.finish(values[y][x],i0+x,j0+y);
    }
    return j.e.ok?0:1;
}
#endif
}
}
#endif
