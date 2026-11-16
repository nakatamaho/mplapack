#include <mpblas.h>
#include <mpblas_mpfr_cuda.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
#define Cherk upper_reference_herk
#include "mpblas/reference/Cherk.cpp"
#undef Cherk
int main() {
 setenv("MPLAPACK_MPFR_CUDA","1",1);setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME","0",1);
 gmp_randstate_t rng;gmp_randinit_default(rng);gmp_randseed_ui(rng,29);
 long checked=0,zeros=0,calls=0;
 for(long p:{256,512,768,1024,2048}) {
  mpfrxx::set_default_precision_bits(p);
  for(long n:{1,7,63,64,65}) for(const char *u:{"U","L"}) {
   const long k=17,lda=std::max(n,k)+3,ldc=n+3;
   std::vector<mpc_class>a(lda*(n+1)),ref(ldc*(n+1));
   long serial=0;
   for(auto *v:{&a,&ref})for(auto &z:*v){mpfr_class r,i;mpfr_urandomb(r.get_mpfr_t(),rng);mpfr_urandomb(i.get_mpfr_t(),rng);mpfr_mul_2si(r.get_mpfr_t(),r.get_mpfr_t(),serial%5-2,MPFR_RNDN);if(serial&1)r=-r;if(serial&2)i=-i;z=mpc_class(r,i);++serial;}
   for(size_t j=0;j<a.size();j+=17)a[j]=mpc_class(mpfr_class(-0.0),mpfr_class(0.0));
   mpfr_prec_t pr=ref[0].real_precision(),pi=ref[0].imag_precision();
   if(pr!=p || pi!=p){std::fprintf(stderr,"precision mismatch requested=%ld actual=%ld/%ld\n",p,long(pr),long(pi));return 3;}
   auto got=ref;mpfr_class alpha=-.75,beta=.625;
   upper_reference_herk(u,"C",n,k,alpha,a.data(),lda,beta,ref.data(),ldc);
   if(!Cherk_mpfr_cuda(u,"C",n,k,alpha,a.data(),lda,beta,got.data(),ldc)){std::fprintf(stderr,"declined p=%ld n=%ld %s\n",p,n,u);return 1;}
   for(size_t j=0;j<ref.size();++j)for(int part=0;part<2;++part){
    mpfr_class x=part?ref[j].imag():ref[j].real(),y=part?got[j].imag():got[j].real();
    if(!mpfr_equal_p(x.get_mpfr_t(),y.get_mpfr_t()) || mpfr_get_prec(x.get_mpfr_t())!=mpfr_get_prec(y.get_mpfr_t())){std::fprintf(stderr,"mismatch p=%ld n=%ld %s index=%zu component=%d\n",p,n,u,j,part);return 2;}
    if(mpfr_zero_p(x.get_mpfr_t()) && mpfr_signbit(x.get_mpfr_t())!=mpfr_signbit(y.get_mpfr_t()))++zeros;
    ++checked;
   }
   ++calls;
  }
  std::fprintf(stderr,"p=%ld passed\n",p);
 }
 std::printf("%ld direct HERK calls, %ld components/padding compared, 0 ordinary mismatches, %ld legacy public zero signs\n",calls,checked,zeros);
 gmp_randclear(rng);
}
