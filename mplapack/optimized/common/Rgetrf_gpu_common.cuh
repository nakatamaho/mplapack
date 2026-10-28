// Internal reference-order LU, with next-panel/update stream lookahead.
#ifndef MPLAPACK_GETRF_GPU_COMMON_CUH
#define MPLAPACK_GETRF_GPU_COMMON_CUH
#include "Rpotrf_gpu_common.cuh"
namespace mplapack_lapack_gpu {
template<class E> struct lu_base {
    E *a; long *piv,*info; bool single;
    MPLAPACK_DEV_HD void operator()(long) const {
        if(single) *piv=1;
        long p=*piv-1;
        *info=a[p]==E(0)?1:0;
        if(!*info && p) { E t=a[0]; a[0]=a[p]; a[p]=t; }
    }
};
template<class Ops> MPLAPACK_DEV_HD dreal<Ops> lu_abs(const dreal<Ops> &x) { return abs(x); }
template<class Ops> MPLAPACK_DEV_HD dreal<Ops> lu_abs(const dcomplex<dreal<Ops>> &z) {
    using R=dreal<Ops>;
    R a=abs(z.real()),b=abs(z.imag());
    if(a<b) { R t=a; a=b; b=t; }
    if(a==R(0)) return R(0);
    R t=b/a;
    return a*R(square_root<Ops>((R(1)+t*t).v));
}
template<class E,class R> struct lu_reciprocal {
    E *a,*reciprocal; R sfmin; long *multiply;
    MPLAPACK_DEV_HD void operator()(long) const {
        *multiply=lu_abs(a[0])>=sfmin;
        if(*multiply) *reciprocal=E(1)/a[0];
    }
};
template<class E,class R> struct lu_scale {
    E *a; const E *reciprocal; const long *multiply;
    MPLAPACK_DEV_HD void operator()(long i) const {
        a[i+1]=*multiply?*reciprocal*a[i+1]:a[i+1]/a[0];
    }
};
template<class Ops> class real_lu : public chol_schedule {
    using F=typename Ops::F;
    buffer<F> work[2],constants;
    size_t capacity=0;
  public:
    using E=dreal<Ops>;
    real_lu(long n,long nb) {
        capacity=(size_t)n*nb;
        if(!error) error=work[0].allocate(capacity);
        if(!error) error=work[1].allocate(capacity);
        if(!error) error=constants.allocate(2);
        F c[2]={Ops::from_d(1),Ops::from_d(-1)};
        if(!error) error=transfer(constants.data(),c,2,false);
    }
    ~real_lu() { finish(); }
    void solve(long m,long n,E *a,E *b,long ld) {
        if(error || !m || !n) return;
        mplapack_blas::gpu_lane_args<Ops> g{};
        g.s.op=mplapack_blas::L3_TRSM; g.s.left=1; g.s.upper=0;
        g.s.trans=0; g.s.nounit=0; g.s.alpha_one=1; g.s.m=m; g.s.n=n;
        g.a=reinterpret_cast<F*>(a); g.c=reinterpret_cast<F*>(b);
        g.al=constants.data(); g.lda=ld; g.ldc=ld;
#if defined(__CUDACC__)
        long lanes=mplapack_blas::blas_lanes(g.s);
        mplapack_blas::k_blas_lanes<Ops><<<(unsigned)((lanes+63)/64),64,0,stream()>>>(g,lanes);
        error=cudaGetLastError();
#else
        error=mplapack_blas::blas_device<Ops>(g);
#endif
    }
    void update(long m,long n,long k,E *a,E *b,E *c,long ld) {
        if(error || !m || !n || !k) return;
        using namespace mplapack_gemm;
        // iMlaenv may request nb<=1 (whole-matrix recursion). Bound the
        // alpha*B scratch by splitting output columns, never the dot product.
        if((size_t)k>capacity) { error=1; return; }
        long width=(long)std::min((size_t)n,capacity/(size_t)k);
        for(long col=0;col<n && !error;col+=width) {
            shape s{}; s.m=m; s.n=std::min(width,n-col); s.k=k; s.lda=ld; s.ldb=ld; s.ldc=ld; s.nota=true; s.notb=true;
            gpu_scalars<F> sc{Ops::from_d(-1),Ops::from_d(1),false,true};
            F *aa=reinterpret_cast<F*>(a),*bb=reinterpret_cast<F*>(b+col*ld),*cc=reinterpret_cast<F*>(c+col*ld),*scratch=work[far?1:0].data();
#if defined(__CUDACC__)
            k_alpha_opB<Ops><<<gpu_grid2(s.k,s.n),dim3(GPU_TILE,GPU_TILE),0,stream()>>>(s,sc.alpha,bb,scratch);
            error=cudaGetLastError();
            if(!error) {
                k_tiled<Ops><<<gpu_grid2(s.m,s.n),dim3(GPU_TILE,GPU_TILE),0,stream()>>>(s,sc,aa,scratch,cc,0);
                error=cudaGetLastError();
            }
#else
            error=gemm_device<Ops>(s,sc,aa,bb,cc,scratch);
#endif
        }
    }
};
template<class Ops> class complex_lu : public chol_schedule {
  public:
    using R=dreal<Ops>; using E=dcomplex<R>;
    complex_lu(long,long) {}
    ~complex_lu() { finish(); }
    void solve(long m,long n,E *a,E *b,long ld) {
        job_Ctrsm<Ops> j{};
        j.side[0]='L'; j.uplo[0]='L'; j.ta[0]='N'; j.diag[0]='U';
        j.m=m; j.n=n; j.lda=ld; j.ldb=ld; j.alpha=E(1,0); j.A=a; j.B=b;
        launch(j,j.units());
    }
    void update(long m,long n,long k,E *a,E *b,E *c,long ld) {
        if(!m || !n || !k) return;
        job_Cgemm<Ops> j{}; j.ta[0]='N'; j.tb[0]='N';
        j.m=m; j.n=n; j.k=k; j.lda=ld; j.ldb=ld; j.ldc=ld;
        j.alpha=E(-1,0); j.beta=E(1,0); j.A=a; j.B=b; j.C=c;
        launch(j,j.units());
    }
};
template<class Blas,class R> class lu_context : public Blas {
  public:
    using E=typename Blas::E;
  private:
    buffer<long> status;
    buffer<E> reciprocal;
    R sfmin;
  public:
    lu_context(long n,long nb,R safe):Blas(n,nb),sfmin(safe) {
        if(!this->error) this->error=status.allocate(1);
        if(!this->error) this->error=reciprocal.allocate(1);
    }
    ~lu_context() { this->finish(); }
    void swap(long n,E *a,long ld,long k1,long k2,const long *piv) {
        this->launch(swaps<E>{a,ld,k1,k2,1,piv},n);
    }
    long recursive(long m,long n,E *a,long ld,long *piv) {
        if(this->error || !m || !n) return 0;
        if(m==1 || n==1) {
            if(m!=1) {
#if defined(__CUDACC__)
                k_iamax<<<1,256,0,this->stream()>>>(a,m,piv);
                this->error=cudaGetLastError();
#else
                this->error=iamax(a,m,piv);
#endif
            }
            this->launch(lu_base<E>{a,piv,status.data(),m==1},1);
            long info=this->feedback(status.data());
            if(m>1 && !info) {
                this->launch(lu_reciprocal<E,R>{a,reciprocal.data(),sfmin,status.data()},1);
                this->launch(lu_scale<E,R>{a,reciprocal.data(),status.data()},m-1);
            }
            return info;
        }
        long n1=std::min(m,n)/2,n2=n-n1;
        long info=recursive(m,n1,a,ld,piv);
        swap(n2,a+n1*ld,ld,1,n1,piv);
        this->solve(n1,n2,a,a+n1*ld,ld);
        this->update(m-n1,n2,n1,a+n1,a+n1*ld,a+n1+n1*ld,ld);
        long i=recursive(m-n1,n2,a+n1+n1*ld,ld,piv+n1);
        if(!info && i) info=i+n1;
        this->launch(adjust_pivots{piv+n1,n1},std::min(m,n)-n1);
        swap(n1,a,ld,n1+1,std::min(m,n),piv);
        return info;
    }
};
template<class Context> static long lu(Context &ctx,long m,long n,long nb,typename Context::E *a,long ld,long *piv) {
    long mn=std::min(m,n),info=0;
    if(nb<=1 || nb>=mn) {
        info=ctx.recursive(m,n,a,ld,piv); ctx.finish(); return info;
    }
    for(long j=0;j<mn;j+=nb) {
        long jb=std::min(nb,mn-j);
        long i=ctx.recursive(m-j,jb,a+j+j*ld,ld,piv+j);
        if(!info && i) info=j+i;
        // The previous far update reads the current panel's left block.
        // Join before its row interchanges; panel arithmetic is independent.
        ctx.join_update();
        ctx.launch(adjust_pivots{piv+j,j},jb);
        ctx.swap(j,a,ld,j+1,j+jb,piv);
        for(int part=0;part<2;++part) {
            long c=j+jb+(part?std::min(nb,n-j-jb):0);
            long cols=part?n-c:std::min(nb,n-c);
            ctx.far=part!=0;
            if(cols) {
                ctx.swap(cols,a+c*ld,ld,j+1,j+jb,piv);
                ctx.solve(jb,cols,a+j+j*ld,a+j+c*ld,ld);
                ctx.update(m-j-jb,cols,jb,a+j+jb+j*ld,a+j+c*ld,a+j+jb+c*ld,ld);
            }
            if(!part) ctx.release_update();
        }
        ctx.far=false;
        if(ctx.error) break;
    }
    ctx.finish(); return info;
}
}
#endif
