#include <mpblas_mpfr_cuda.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
int main() {
 setenv("MPLAPACK_MPFR_CUDA","1",1);setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME","0",1);
 for(long p:{256,512,768,1024,2048}) for(long n:{1,7,13}) {
  mpfrxx::set_default_precision_bits(p);std::vector<mpfr_class>a(n*n,mpfr_class(0)),b(n*3,mpfr_class(1));mpfr_class al=1;
  for(long i=0;i<n;++i)a[i+i*n]=2;
  std::fprintf(stderr,"begin p=%ld n=%ld\n",p,n);
  if(!Rtrsm_mpfr_cuda("L","U","T","N",n,3,al,a.data(),n,b.data(),n)) return 1;
  for(const auto &v:b)if(v!=mpfr_class(0.5))return 2;
  std::fprintf(stderr,"pass p=%ld n=%ld\n",p,n);
 }
}
