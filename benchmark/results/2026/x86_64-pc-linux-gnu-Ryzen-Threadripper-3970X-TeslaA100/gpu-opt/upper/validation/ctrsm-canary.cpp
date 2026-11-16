#include <mpblas_mpfr_cuda.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
int main() {
 setenv("MPLAPACK_MPFR_CUDA","1",1);setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME","0",1);
 for(long p:{256,512,768,1024,2048}) for(long n:{1,13,64}) {
  mpfrxx::set_default_precision_bits(p);mpc_class diag(mpfr_class(1.25),mpfr_class(.125)),al(1),expected=al/conj(diag);
  std::vector<mpc_class>a(n*n,mpc_class(0)),b(n*3,mpc_class(1));
  for(long i=0;i<n;++i)a[i+i*n]=diag;
  mpfr_prec_t pr=b[0].real_precision(),pi=b[0].imag_precision();if(pr!=p || pi!=p)return 3;
  std::fprintf(stderr,"begin p=%ld n=%ld\n",p,n);
  if(!Ctrsm_mpfr_cuda("L","U","C","N",n,3,al,a.data(),n,b.data(),n)) return 1;
  for(const auto &v:b)if(v!=expected)return 2;
  std::fprintf(stderr,"pass p=%ld n=%ld\n",p,n);
 }
}
