// MPFR oracle for the additive mpc_cuda cu_fround patch. Plain C++ host
// build, or included by cu_fround_device_test.cu for the identical GPU test.
#include <mpc_cuda/cu_fround.cuh>
#include <mpfr.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#ifdef CU_FROUND_DEVICE_TEST
#include <cuda_runtime.h>
#define CUDA_OK(x) do { cudaError_t e=(x); if(e!=cudaSuccess) { \
    std::fprintf(stderr,"%s: %s\n",#x,cudaGetErrorString(e)); std::exit(1); } } while(0)
#endif
using namespace cu_fp;

template<int P> struct sample { cu_freal<P> a,b; };
template<int P> struct result { cu_freal<P> q,r; int div_status,sqrt_status; bool de,se; int ties; };

template<int P> __host__ __device__ static int midpoint_tests()
{
    // Same-P finite binary operands cannot produce exact div/sqrt midpoints:
    // an odd (P+1)-bit midpoint numerator (or its square) needs >P bits.
    // Test those branches independently using exact wider integers, not
    // operands rounded to P before the operation (which destroys the tie).
    const int W=P+64, N=cu_freal<W>::N;
    int bad=0;
    for(int odd=0;odd<2;++odd) {
        cu_freal<W> t=cu_convert<W>(cu_freal<P>(1.0));
        // t = 1 + (odd+1/2)*ulp(P), hence both tie parities.
        t.m[(cu_freal<W>::SB+63)/64] |= 1ULL<<((cu_freal<W>::SB+63)%64);
        if(odd) t.m[(cu_freal<W>::SB+64)/64] |= 1ULL<<((cu_freal<W>::SB+64)%64);
        cu_freal<W> one(1.0);
        bad += fround_detail::residual(t,t,one)!=0;
        for(int dir=-1;dir<=1;++dir) {
            cu_freal<W> numerator=t;
            if(dir) numerator.m[0]+=dir>0 ? (1ULL<<cu_freal<W>::SB) : -(1ULL<<cu_freal<W>::SB);
            // With SB=32 the midpoint's low limb is zero: propagate borrow.
            if(dir<0 && cu_freal<W>::SB) --numerator.m[1];
            bad += fround_detail::residual(numerator,t,one)!=dir;
            cu_freal<P> got=fround_detail::finish<P>(t,dir);
            cu_freal<W> lo=t,hi=t;
            const int bit=cu_freal<W>::SB+63;
            lo.m[bit/64] &= ~(1ULL<<(bit%64));
            hi=cu_fadd<W>(lo,cu_scale2<W>(one,1-P));
            cu_freal<P> want=cu_convert<P>(dir<0 || (dir==0 && !odd)?lo:hi);
            bad += cu_cmp<P>(got,want)!=0;
            // a = t^2 +/- one integer unit: compare exact 2N-limb squares.
            cu_limb square[2*N], a[2*N]; cu_mul_full<W>(square,t.m,t.m);
            for(int k=0;k<2*N;++k) a[k]=square[k];
            if(dir) { cu_limb carry=1;
                for(int k=0;k<2*N && carry;++k) {
                    cu_limb old=a[k]; a[k]=dir>0?old+carry:old-carry;
                    carry=dir>0?(a[k]<old):(old<carry);
                }
            }
            bad += fround_detail::compare(a,2*N,0,square,2*N,0)!=dir;
        }
    }
    return bad;
}

template<int P> __host__ __device__ static int contracts()
{
    cu_freal<P> x(-1.0), y(0.0), v;
    int bad=cu_sqrt_rn<P>(v,x)==0;
    bad+=cu_fdiv_rn<P>(v,x,y)==0;
    bad+=cu_fdiv_rn<P>(x,x,x)!=0 || cu_cmp<P>(x,cu_freal<P>(1.0))!=0;
    bad+=cu_sqrt_rn<P>(x,x)!=0 || cu_cmp<P>(x,cu_freal<P>(1.0))!=0;
    return bad;
}

template<int P> __host__ __device__ static result<P> evaluate(const sample<P>& s)
{
    result<P> o;
    o.q.set_zero(); o.r.set_zero(); // domain errors need not initialize output
    o.div_status=fround_detail::div(o.q,s.a,s.b,o.de);
    cu_freal<P> a=cu_abs<P>(s.a);
    o.sqrt_status=fround_detail::sqrt(o.r,a,o.se);
    o.ties=0;
    return o;
}
#ifdef CU_FROUND_DEVICE_TEST
template<int P> __global__ void kernel(const sample<P>* s,result<P>* o,int n)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n) { o[i]=evaluate(s[i]); if(i==0) o[i].ties=midpoint_tests<P>()+contracts<P>(); }
}
#endif

template<int P> static void to_mpfr(mpfr_ptr x,const cu_freal<P>& a)
{
    if(a.is_zero()) { mpfr_set_zero(x,1); return; }
    mp_limb_t m[cu_freal<P>::N];
    for(int k=0;k<cu_freal<P>::N;++k) m[k]=a.m[k];
    mpfr_t v;
    mpfr_custom_init_set(v,a.sign*MPFR_REGULAR_KIND,a.exp,P,m);
    mpfr_set(x,v,MPFR_RNDN);
}
template<int P> static cu_freal<P> from_mpfr(mpfr_srcptr a)
{
    cu_freal<P> x;
    if(mpfr_zero_p(a)) { x.set_zero(); return x; }
    x.sign=mpfr_signbit(a)?-1:1; x.exp=mpfr_get_exp(a);
    auto m=static_cast<const mp_limb_t*>(mpfr_custom_get_significand(a));
    for(int k=0;k<cu_freal<P>::N;++k) x.m[k]=m[k];
    return x;
}
template<int P> static cu_freal<P> random_value(std::mt19937_64& rng)
{
    cu_freal<P> a;
    a.sign=(rng()&1)?1:-1; a.exp=long(rng()%2000001)-1000000;
    for(int k=0;k<a.N;++k) a.m[k]=rng();
    a.m[a.N-1]|=1ULL<<63;
    if(a.SB) a.m[0]&=~0xffffffffULL;
    return a;
}
template<int P> static int run(size_t cases)
{
    std::mt19937_64 rng(0x07a00000+P);
    const int batch=4096;
    std::vector<sample<P>> s(batch);
    std::vector<result<P>> o(batch);
    size_t bad=0, de=0,se=0, checked=0;
    mpfr_t a,b,q,r,g;
    mpfr_inits2(P,a,b,q,r,g,(mpfr_ptr)0);
    mpz_t z; mpz_init(z);
    // 256 hard pairs, then random pairs. Near-midpoint quotients on BOTH
    // sides: b=floor/ceil(2^P/3), a=b+1, q near 1+3*2^-P (a PB midpoint).
    // sqrt(1+2^(1-P)) is just below 1+2^-P by about 2^(-2P-1).
#ifdef CU_FROUND_DEVICE_TEST
    sample<P>* ds; result<P>* dout;
    CUDA_OK(cudaMalloc(&ds,batch*sizeof(sample<P>)));
    CUDA_OK(cudaMalloc(&dout,batch*sizeof(result<P>)));
#else
    bad+=midpoint_tests<P>()+contracts<P>();
#endif
    for(size_t first=0;first<cases+256;first+=batch) {
        int n=int(std::min(size_t(batch),cases+256-first));
        for(int i=0;i<n;++i) {
            size_t id=first+i;
            s[i].a=random_value<P>(rng); s[i].b=random_value<P>(rng);
            if(id<256) {
                int kind=id%16;
                mpfr_set_ui(a,1,MPFR_RNDN); mpfr_set_ui(b,1,MPFR_RNDN);
                if(kind<2) {
                    mpz_set_ui(z,1); mpz_mul_2exp(z,z,P); mpz_fdiv_q_ui(z,z,3);
                    if(kind&1) mpz_add_ui(z,z,1);
                    mpfr_set_z(b,z,MPFR_RNDN); mpfr_add_ui(a,b,1,MPFR_RNDN);
                } else if(kind==2 || kind==3) {
                    // Random odd denominator d with P-1 bits. An inverse
                    // modulo 2^P constructs midpoint k/2^P and residual +/-1
                    // exactly, with a numerator that still fits P bits.
                    mpz_t d,mod,k,num; mpz_inits(d,mod,k,num,(mpz_ptr)0);
                    std::vector<unsigned long long> limbs((P+63)/64);
                    for(auto& limb:limbs) limb=rng();
                    mpz_import(d,limbs.size(),-1,sizeof(limbs[0]),0,0,limbs.data());
                    mpz_fdiv_r_2exp(d,d,P-1); mpz_setbit(d,P-2); mpz_setbit(d,0);
                    mpz_set_ui(mod,1); mpz_mul_2exp(mod,mod,P);
                    mpz_invert(k,d,mod);
                    if(kind==3) mpz_sub(k,mod,k);
                    mpz_add(k,k,mod); mpz_mul(num,k,d);
                    if(kind==3) mpz_add_ui(num,num,1); else mpz_sub_ui(num,num,1);
                    mpz_fdiv_q_2exp(num,num,P);
                    mpfr_set_z(a,num,MPFR_RNDN); mpfr_set_z(b,d,MPFR_RNDN);
                    mpz_clears(d,mod,k,num,(mpz_ptr)0);
                } else if(kind==4) {
                    mpfr_set_ui_2exp(b,1,1-P,MPFR_RNDN);
                    mpfr_add(a,a,b,MPFR_RNDN); mpfr_set_ui(b,1,MPFR_RNDN);
                } else if(kind==5) {
                    // Hensel-lift k^2 == -7 (mod 2^(P+2)). Choose k with
                    // P+1 bits: t=k/2^P is a PB midpoint, and
                    // a=t^2+7*2^(-2P) fits PB exactly, just ABOVE the boundary.
                    mpz_t k,tmp,mod; mpz_inits(k,tmp,mod,(mpz_ptr)0);
                    mpz_set_ui(k,1);
                    for(int bit=3;bit<P+2;++bit) {
                        mpz_mul(tmp,k,k); mpz_add_ui(tmp,tmp,7);
                        if(mpz_tstbit(tmp,bit)) mpz_setbit(k,bit-1);
                    }
                    mpz_set_ui(mod,1); mpz_mul_2exp(mod,mod,P+1);
                    if(!mpz_tstbit(k,P)) mpz_sub(k,mod,k);
                    mpz_mul(tmp,k,k); mpz_add_ui(tmp,tmp,7);
                    mpfr_set_z(a,tmp,MPFR_RNDN); mpfr_mul_2si(a,a,-2*P,MPFR_RNDN);
                    mpz_clears(k,tmp,mod,(mpz_ptr)0);
                } else if(kind==6) mpfr_set_zero(a,1);
                else if(kind==7) mpfr_set_zero(b,1);
                else if(kind==8) { mpfr_nextbelow(a); mpfr_set(b,a,MPFR_RNDN); }
                else if(kind==9) { mpfr_set_ui(a,9,MPFR_RNDN); mpfr_set_ui(b,3,MPFR_RNDN); }
                else if(kind==10) { mpfr_set_ui(a,25,MPFR_RNDN); mpfr_set_ui(b,5,MPFR_RNDN); }
                else if(kind==11) { mpfr_nextabove(a); mpfr_nextbelow(b); }
                else if(kind==12) { mpfr_set_ui_2exp(a,1,mpfr_get_emax()-16,MPFR_RNDN); mpfr_set(b,a,MPFR_RNDN); }
                else if(kind==13) { mpfr_set_ui_2exp(a,1,mpfr_get_emin()+16,MPFR_RNDN); mpfr_set(b,a,MPFR_RNDN); }
                else if(kind==14) { mpfr_nextbelow(a); mpfr_nextbelow(b); }
                else if(kind==15) {
                    std::vector<unsigned long long> limbs((P/2+63)/64);
                    for(auto& limb:limbs) limb=rng();
                    mpz_import(z,limbs.size(),-1,sizeof(limbs[0]),0,0,limbs.data());
                    mpz_fdiv_r_2exp(z,z,P/2-2); mpz_setbit(z,P/2-3);
                    mpfr_set_z(b,z,MPFR_RNDN); mpfr_mul(a,b,b,MPFR_RNDN); // exact square, a/b=b
                }
                // Exact powers of two, exact quotients, signs, exponent parity.
                if(id&16) mpfr_neg(a,a,MPFR_RNDN);
                if(id&32) mpfr_neg(b,b,MPFR_RNDN);
                if(id&64 && kind<12) { mpfr_mul_2si(a,a,-127,MPFR_RNDN); mpfr_mul_2si(b,b,-127,MPFR_RNDN); }
                s[i].a=from_mpfr<P>(a); s[i].b=from_mpfr<P>(b);
            }
        }
#ifdef CU_FROUND_DEVICE_TEST
        CUDA_OK(cudaMemcpy(ds,s.data(),n*sizeof(sample<P>),cudaMemcpyHostToDevice));
        kernel<P><<<(n+63)/64,64>>>(ds,dout,n);
        CUDA_OK(cudaGetLastError());
        CUDA_OK(cudaMemcpy(o.data(),dout,n*sizeof(result<P>),cudaMemcpyDeviceToHost));
#else
        for(int i=0;i<n;++i) o[i]=evaluate(s[i]);
#endif
        for(int i=0;i<n;++i) {
            auto& v=o[i]; bad+=v.ties; de+=v.de; se+=v.se;
            size_t id=first+i;
            if(id<256 && !(id&64)) {
                if(id%16<4 && !v.de) ++bad;
                if((id%16==4 || id%16==5) && !v.se) ++bad;
            }
            to_mpfr<P>(a,s[i].a); to_mpfr<P>(b,s[i].b);
            if(mpfr_zero_p(b)) bad+=v.div_status==0;
            else {
                mpfr_div(q,a,b,MPFR_RNDN); to_mpfr<P>(g,v.q);
                if(v.div_status || !mpfr_equal_p(q,g)) {
                    if(bad<5) std::fprintf(stderr,"PB=%d div case=%zu\n",P,first+i);
                    ++bad;
                }
                ++checked;
            }
            mpfr_abs(a,a,MPFR_RNDN); mpfr_sqrt(r,a,MPFR_RNDN); to_mpfr<P>(g,v.r);
            if(v.sqrt_status || !mpfr_equal_p(r,g)) {
                if(bad<5) std::fprintf(stderr,"PB=%d sqrt case=%zu\n",P,first+i);
                ++bad;
            }
            ++checked;
        }
    }
    bad+=de==0 || se==0;
    std::printf("PB=%d pairs=%zu checked=%zu mismatches=%zu exact_div=%zu exact_sqrt=%zu\n",P,cases+256,checked,bad,de,se);
    std::fflush(stdout);
#ifdef CU_FROUND_DEVICE_TEST
    CUDA_OK(cudaFree(ds)); CUDA_OK(cudaFree(dout));
#endif
    mpz_clear(z); mpfr_clears(a,b,q,r,g,(mpfr_ptr)0);
    return bad?1:0;
}
int main(int argc,char** argv)
{
#ifdef CU_FROUND_DEVICE_TEST
    int n=0; if(cudaGetDeviceCount(&n)!=cudaSuccess || !n) return 77;
    CUDA_OK(cudaDeviceSetLimit(cudaLimitStackSize,65536));
    size_t cases=1000000;
#else
    size_t cases=100000;
#endif
    if(argc==2) cases=std::strtoull(argv[1],nullptr,10);
    int bad=0;
    bad|=run<224>(cases); bad|=run<256>(cases); bad|=run<512>(cases);
    bad|=run<768>(cases); bad|=run<1024>(cases); bad|=run<2048>(cases);
    return bad;
}
