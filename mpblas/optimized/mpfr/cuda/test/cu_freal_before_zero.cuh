/* cu_freal.cuh -- fixed-precision, register-resident binary floating-point real
 * for CUDA (host + device).  The precision PB (in bits) is a COMPILE-TIME
 * template parameter, a multiple of 32, so the significand lives in registers
 * (no global-memory arena, no runtime-precision dispatch).  Arithmetic is
 * bit-exact with MPFR round-to-nearest-even (RNDN).
 *
 *   cu_freal<256> a, b, c;          // 256-bit mantissa
 *   a = 1.5; b = some_double;
 *   c = a * b + a;                  // operators, or cu_fmul/cu_fadd
 *   double d = (double) c;          // or c.to_double()
 *
 * Storage (MPFR-compatible, left-justified):
 *   value = sign * M * 2^(exp - N*64),   M = sum m[i]*2^(64 i),  m[N-1] MSB set,
 *   N = ceil(PB/64) limbs, the low SB = N*64-PB bits of M are zero (SB in {0,32}).
 *   exp follows the MPFR convention (EXP of the value); zero is a special flag.
 */
// Frozen abd9407 oracle; intentionally preserves the upstream zero-sign bug.
#ifndef MPLAPACK_TEST_FREAL_BEFORE_ZERO_CUH
#define MPLAPACK_TEST_FREAL_BEFORE_ZERO_CUH

#include <cmath>
#include <cstdint>

/* Usable from plain host C++ too (no nvcc): make the execution-space tags
 * no-ops when not compiling CUDA.  Under nvcc they are already defined. */
#if !defined(__CUDACC__)
#  ifndef __host__
#    define __host__
#  endif
#  ifndef __device__
#    define __device__
#  endif
#endif

namespace cu_fp_before_zero {

typedef unsigned long long cu_limb;
#define CU_FREAL_ZERO_EXP (-0x7fffffffffffffLL)

/* ---------------- 64x64 -> 128 ---------------- */
__host__ __device__ static inline void
cu_umul (cu_limb &hi, cu_limb &lo, cu_limb a, cu_limb b)
{
#if defined(__CUDA_ARCH__)
  asm ("mul.lo.u64 %0,%2,%3;\n\tmul.hi.u64 %1,%2,%3;"
       : "=l"(lo),"=l"(hi) : "l"(a),"l"(b));
#else
  __uint128_t p = (__uint128_t)a*b; lo=(cu_limb)p; hi=(cu_limb)(p>>64);
#endif
}

/* ---------------- generic bit helpers over a limb array ---------------- */
__host__ __device__ static inline cu_limb
cu_wbit (const cu_limb *W, int len, long i)
{ if (i<0||i>=(long)len*64) return 0; return (W[i>>6]>>(i&63))&1; }

__host__ __device__ static inline cu_limb
cu_wsticky_below (const cu_limb *W, int len, long i)      /* OR of bits [0..i-1] */
{
  if (i<=0) return 0;
  long full=i>>6; cu_limb s=0;
  for (long j=0;j<full && j<len;j++) s|=W[j];
  int rem=(int)(i&63); if (rem && full<len) s|= W[full] & (((cu_limb)1<<rem)-1);
  return s!=0;
}

/* ============================ the type ============================ */
template<int PB>
struct cu_freal
{
  static_assert (PB>0 && (PB%32)==0, "PB must be a positive multiple of 32");
  static const int P  = PB;
  static const int N  = (PB+63)/64;         /* limbs                       */
  static const int SB = N*64 - PB;          /* unused low bits, 0 or 32     */

  int     sign;                             /* +1 / -1                      */
  long    exp;                              /* MPFR exponent; ZERO flag     */
  cu_limb m[N];                             /* m[N-1] MSB set when nonzero  */

  __host__ __device__ bool is_zero () const { return exp==CU_FREAL_ZERO_EXP; }
  __host__ __device__ void set_zero () { sign=1; exp=CU_FREAL_ZERO_EXP;
    for(int i=0;i<N;i++) m[i]=0; }

  /* ---- finalize: normalize+round a result into this type ----
   * D[len]    : little-endian magnitude buffer (the significant bits)
   * msb       : index of the result MSB within D (>=0; <0 means zero)
   * Ehi       : value-weight exponent of bit `msb`  (value bit msb has weight 2^Ehi)
   * brnd,bsti : round/sticky material strictly BELOW D bit 0
   * Rounds to P bits (low SB bits forced zero), RNDN.                    */
  __host__ __device__ static cu_freal
  finalize (const cu_limb *D, int len, int msb, long Ehi,
            cu_limb brnd, cu_limb bsti, int sgn)
  {
    cu_freal r;
    if (msb<0){ r.set_zero(); return r; }
    long Eh = Ehi;
    long lo = (long)msb - P + 1;           /* keep bits [lo .. msb]        */
    cu_limb V[N];                          /* P-bit value (top SB bits 0)  */
    cu_limb rbit, sticky;
    if (lo>0)
      {
        rbit   = cu_wbit (D,len,lo-1);
        sticky = brnd | bsti | cu_wsticky_below (D,len,lo-1);
        long ws=lo>>6; int bs=(int)(lo&63);
        for (int i=0;i<N;i++){ long s=ws+i;
          cu_limb a=(s>=0&&s<len)?D[s]:0, b=(s+1>=0&&s+1<len)?D[s+1]:0;
          V[i]= bs? ((a>>bs)|(b<<(64-bs))) : a; }
      }
    else if (lo==0)
      {
        rbit = brnd; sticky = bsti;
        for (int i=0;i<N;i++) V[i]=(i<len)?D[i]:0;
      }
    else /* lo<0: fewer than P significant bits -> exact, left-justify */
      {
        rbit=0; sticky=0;
        int ls=(int)(-lo); int ws=ls>>6, bs=ls&63;
        for (int i=N-1;i>=0;i--){ int s=i-ws;
          cu_limb a=(s>=0&&s<len)?D[s]:0, b=(s-1>=0&&s-1<len)?D[s-1]:0;
          V[i]= bs? ((a<<bs)|(b>>(64-bs))) : a; }
      }
    /* round to nearest even at V bit 0 */
    if (rbit && (sticky || (V[0]&1)))
      {
        cu_limb c=1;
        for (int i=0;i<N && c;i++){ V[i]+=c; c=(V[i]==0); }
        bool ov = SB ? (((V[(P)>>6]>>(P&63))&1)!=0) : (c!=0);
        if (ov){ for(int i=0;i<N;i++) V[i]=0;
                 V[(P-1)>>6]=(cu_limb)1<<((P-1)&63); Eh++; }
      }
    /* M = V << SB  (left-justify; low SB bits zero) */
    if (SB){
      for (int i=N-1;i>=0;i--){
        cu_limb a=V[i], b=(i-1>=0)?V[i-1]:0;
        r.m[i]=(a<<SB)|(b>>(64-SB)); }
    } else { for (int i=0;i<N;i++) r.m[i]=V[i]; }
    r.sign=sgn; r.exp=Eh+1;
    return r;
  }

  /* ---------------- double <-> ---------------- */
  __host__ __device__ static cu_freal from_double (double d)
  {
    cu_freal r;
    if (d==0.0){ r.set_zero(); return r; }
    int sgn = d<0?-1:1; d = d<0?-d:d;
    int k; double f = frexp (d,&k);              /* d=f*2^k, f in [0.5,1)   */
    cu_limb F = (cu_limb) ldexp (f,53);          /* 53-bit, MSB at bit 52   */
    return finalize (&F, 1, 52, (long)k-1, 0,0, sgn);
  }
  __host__ __device__ double to_double () const
  {
    if (is_zero()) return 0.0;
    /* top 53 bits of M (MSB at N*64-1) */
    cu_limb F, rbit, sticky;
    int hb = N*64-53;                            /* bit index of F's LSB    */
    int w=hb>>6, b=hb&63;
    if (b){ F = (m[w]>>b) | (w+1<N? m[w+1]<<(64-b):0); }
    else    F = m[w];
    long rp = hb-1;                              /* round bit index         */
    rbit = cu_wbit (m,N,rp);
    sticky = cu_wsticky_below (m,N,rp);
    long e = exp;
    if (rbit && (sticky || (F&1))){ F++; if (F==(1ULL<<53)){ F>>=1; e++; } }
    return (double)sign * ldexp ((double)F, e-53);
  }

  /* ---- ergonomic ctors / conversions ---- */
  __host__ __device__ cu_freal () {}
  __host__ __device__ cu_freal (double d) { *this = from_double(d); }
  __host__ __device__ explicit operator double () const { return to_double(); }
};

/* ---------------- N x N -> 2N schoolbook (register-resident) ---------- */
template<int PB> __host__ __device__ static void
cu_mul_full (cu_limb r[2*cu_freal<PB>::N],
             const cu_limb a[cu_freal<PB>::N], const cu_limb b[cu_freal<PB>::N])
{
  const int N = cu_freal<PB>::N;
  for (int i=0;i<2*N;i++) r[i]=0;
  for (int j=0;j<N;j++)
    {
      cu_limb cl=0;
      for (int i=0;i<N;i++)
        {
          cu_limb hi,lo,rl=r[i+j];
#if defined(__CUDA_ARCH__)
          asm ("mad.lo.cc.u64 %0,%2,%3,%4;\n\t"
               "madc.hi.u64   %1,%2,%3,0;\n\t"
               "add.cc.u64    %0,%0,%5;\n\t"
               "addc.u64      %1,%1,0;"
               : "=&l"(lo),"=&l"(hi) : "l"(a[i]),"l"(b[j]),"l"(cl),"l"(rl));
#else
          cu_umul(hi,lo,a[i],b[j]);
          cu_limb t=lo+cl; hi+=(t<lo); lo=t; t=lo+rl; hi+=(t<lo); lo=t;
#endif
          r[i+j]=lo; cl=hi;
        }
      r[N+j]=cl;
    }
}

/* ---------------- multiply ---------------- */
template<int PB> __host__ __device__ static cu_freal<PB>
cu_fmul (const cu_freal<PB> &a, const cu_freal<PB> &b)
{
  typedef cu_freal<PB> F; const int N=F::N;
  if (a.is_zero()||b.is_zero()){ F r; r.set_zero(); return r; }
  cu_limb P2[2*N];
  cu_mul_full<PB> (P2, a.m, b.m);
  int msb = (P2[2*N-1] & (1ULL<<63)) ? 2*N*64-1 : 2*N*64-2;
  long Ehi = a.exp + b.exp - 2*N*64 + msb;
  return F::finalize (P2, 2*N, msb, Ehi, 0,0, a.sign*b.sign);
}

/* ---------------- magnitude compare ---------------- */
template<int PB> __host__ __device__ static int
cu_cmpmag (const cu_freal<PB> &a, const cu_freal<PB> &b)
{
  if (a.exp!=b.exp) return a.exp>b.exp?1:-1;
  for (int i=cu_freal<PB>::N-1;i>=0;i--) if (a.m[i]!=b.m[i]) return a.m[i]>b.m[i]?1:-1;
  return 0;
}

/* ---------------- magnitude add: |a|+|b| ---------------- */
template<int PB> __host__ __device__ static cu_freal<PB>
cu_addmag (cu_freal<PB> a, cu_freal<PB> b, int sgn)
{
  typedef cu_freal<PB> F; const int N=F::N;
  if (b.exp>a.exp){ F t=a; a=b; b=t; }
  long d = a.exp-b.exp;
  cu_limb Bsh[N]; cu_limb brnd, bsti;
  if (d>=N*64){ for(int i=0;i<N;i++) Bsh[i]=0; }
  else { long ws=d>>6; int bs=(int)(d&63);
    for (int i=0;i<N;i++){ long s=ws+i;
      cu_limb x=(s<N)?b.m[s]:0, y=(s+1<N)?b.m[s+1]:0;
      Bsh[i]= bs? ((x>>bs)|(y<<(64-bs))) : x; } }
  brnd = cu_wbit (b.m,N,d-1); bsti = cu_wsticky_below (b.m,N,d-1);
  /* sum (N+1 limbs) */
  cu_limb D[N+1]; cu_limb c=0;
  for (int i=0;i<N;i++){ cu_limb s=a.m[i]+c; c=(s<c); s+=Bsh[i]; c+=(s<Bsh[i]); D[i]=s; }
  D[N]=c;
  int msb = c? N*64 : N*64-1;
  long Ehi = a.exp - N*64 + msb;
  return F::finalize (D, N+1, msb, Ehi, brnd, bsti, sgn);
}

/* ---------------- magnitude sub: |a|-|b|, |a|>|b| ---------------- */
template<int PB> __host__ __device__ static cu_freal<PB>
cu_submag (const cu_freal<PB> &a, const cu_freal<PB> &b, int sgn)
{
  typedef cu_freal<PB> F; const int N=F::N; const int WL=N+2;   /* 128-bit guard */
  long d = a.exp-b.exp;
  cu_limb A[WL], B[WL]; cu_limb ext=0;
  for (int i=0;i<WL;i++){ A[i]=0; B[i]=0; }
  for (int i=0;i<N;i++) A[i+2]=a.m[i];
  if (d < WL*64){
    long ws=d>>6; int bs=(int)(d&63);
    cu_limb T[WL]; for(int i=0;i<WL;i++) T[i]=0;
    for (int i=0;i<N;i++) T[i+2]=b.m[i];
    for (long i=0;i<ws && i<WL;i++) ext |= (T[i]!=0);
    if (bs && ws<WL) ext |= (T[ws] & (((cu_limb)1<<bs)-1))!=0;
    for (int i=0;i<WL;i++){ long s=ws+i;
      cu_limb x=(s<WL)?T[s]:0, y=(s+1<WL)?T[s+1]:0;
      B[i]= bs? ((x>>bs)|(y<<(64-bs))) : x; }
  } else ext=1;
  /* D = A - B  (- borrow if ext, see note in 1024 proto) */
  cu_limb D[WL]; cu_limb br=0;
  for (int i=0;i<WL;i++){ cu_limb ai=A[i], bi=B[i]; cu_limb s=ai-bi; cu_limb bb=(ai<bi);
    cu_limb s2=s-br; bb|=(s<br); D[i]=s2; br=bb; }
  if (ext){ cu_limb bb=1; for(int i=0;i<WL && bb;i++){ cu_limb s=D[i]-bb; bb=(D[i]<bb); D[i]=s; } }
  int hl=-1; for (int i=WL-1;i>=0;i--){ if (D[i]){ hl=i; break; } }
  if (hl<0){ F r; r.set_zero(); return r; }
  int top=63; { cu_limb v=D[hl]; while(!((v>>top)&1)) top--; }
  int msb = hl*64+top;
  long Ehi = a.exp - WL*64 + msb;
  return F::finalize (D, WL, msb, Ehi, 0, ext, sgn);
}

/* ---------------- add / sub ---------------- */
template<int PB> __host__ __device__ static cu_freal<PB>
cu_fadd (const cu_freal<PB> &a, const cu_freal<PB> &b)
{
  if (a.is_zero()) return b;
  if (b.is_zero()) return a;
  if (a.sign==b.sign) return cu_addmag<PB>(a,b,a.sign);
  int c = cu_cmpmag<PB>(a,b);
  if (c==0){ cu_freal<PB> r; r.set_zero(); return r; }
  return (c>0)? cu_submag<PB>(a,b,a.sign) : cu_submag<PB>(b,a,b.sign);
}
template<int PB> __host__ __device__ static cu_freal<PB>
cu_fsub (const cu_freal<PB> &a, cu_freal<PB> b)
{ b.sign = -b.sign; return cu_fadd<PB>(a,b); }

/* ---------------- utilities (sign / scale / convert / int) ---------------- */
template<int PB> __host__ __device__ static inline cu_freal<PB>
cu_neg (cu_freal<PB> a){ if(!a.is_zero()) a.sign=-a.sign; return a; }
template<int PB> __host__ __device__ static inline cu_freal<PB>
cu_abs (cu_freal<PB> a){ if(!a.is_zero()) a.sign=1; return a; }
template<int PB> __host__ __device__ static inline cu_freal<PB>
cu_scale2 (cu_freal<PB> a, long k){ if(!a.is_zero()) a.exp+=k; return a; }  /* a*2^k */
template<int PB> __host__ __device__ static inline int
cu_sgn (const cu_freal<PB>&a){ return a.is_zero()?0:a.sign; }

/* full signed compare: -1 if a<b, 0 if equal, 1 if a>b */
template<int PB> __host__ __device__ static int
cu_cmp (const cu_freal<PB>&a, const cu_freal<PB>&b)
{
  if (a.is_zero()&&b.is_zero()) return 0;
  if (a.is_zero()) return b.sign>0?-1:1;
  if (b.is_zero()) return a.sign>0?1:-1;
  if (a.sign!=b.sign) return a.sign>0?1:-1;
  int m = cu_cmpmag<PB>(a,b);                 /* magnitude */
  return a.sign>0? m : -m;
}

/* re-round a value of one precision to another (exact if widening) */
template<int TO, int FROM> __host__ __device__ static inline cu_freal<TO>
cu_convert (const cu_freal<FROM>&x)
{
  if (x.is_zero()){ cu_freal<TO> r; r.set_zero(); return r; }
  return cu_freal<TO>::finalize (x.m, cu_freal<FROM>::N,
                                 cu_freal<FROM>::N*64-1, x.exp-1, 0,0, x.sign);
}

/* exact small-integer -> cu_freal */
template<int PB> __host__ __device__ static cu_freal<PB>
cu_from_si (long v)
{
  cu_freal<PB> r;
  if (v==0){ r.set_zero(); return r; }
  int sgn = v<0?-1:1; unsigned long long a = v<0? (unsigned long long)(-(v+1))+1ULL : (unsigned long long)v;
  int nb=64; while(!((a>>(nb-1))&1)) nb--;       /* bit length of a */
  return cu_freal<PB>::finalize (&a, 1, nb-1, (long)nb-1, 0,0, sgn);
}

/* round to nearest integer, return as long (assumes |round(x)| < 2^62) */
template<int PB> __host__ __device__ static long
cu_lround (const cu_freal<PB>&x)
{
  if (x.is_zero()) return 0;
  const int N=cu_freal<PB>::N; long e=x.exp;
  if (e<=0){                                     /* |x|<1 */
    if (e==0){                                   /* [0.5,1): rounds to 1 unless exactly .5 (->0) */
      bool half=(x.m[N-1]==(1ULL<<63)); for(int i=0;i<N-1;i++) half&=(x.m[i]==0);
      long r=half?0:1; return x.sign<0?-r:r;
    }
    return 0;
  }
  unsigned long long ip=0;
  for (long b=0;b<e;b++){ long bit=(long)N*64-e+b; ip |= ((x.m[bit>>6]>>(bit&63))&1ULL)<<b; }
  long rb=(long)N*64-e-1;
  unsigned long long rbit = rb>=0? ((x.m[rb>>6]>>(rb&63))&1ULL):0;
  unsigned long long st = cu_wsticky_below (x.m, N, rb);
  if (rbit && (st || (ip&1))) ip++;
  long r=(long)ip; return x.sign<0?-r:r;
}

/* ---------------- operator sugar ---------------- */
template<int PB> __host__ __device__ static inline cu_freal<PB>
operator* (const cu_freal<PB>&a,const cu_freal<PB>&b){ return cu_fmul<PB>(a,b); }
template<int PB> __host__ __device__ static inline cu_freal<PB>
operator+ (const cu_freal<PB>&a,const cu_freal<PB>&b){ return cu_fadd<PB>(a,b); }
template<int PB> __host__ __device__ static inline cu_freal<PB>
operator- (const cu_freal<PB>&a,const cu_freal<PB>&b){ return cu_fsub<PB>(a,b); }

} /* namespace cu_fp */
#endif /* CU_MPC_CUDA_FREAL_CUH */
