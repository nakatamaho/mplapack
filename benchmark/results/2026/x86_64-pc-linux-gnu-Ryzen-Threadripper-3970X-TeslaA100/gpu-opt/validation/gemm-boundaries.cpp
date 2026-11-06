// Supplement the regular suite with large/edge dimensions at changed precisions.
// Reuse its exact reference, random operands, padding and bit comparison logic.
// Never use a shared library's shadowing Rgemm dispatcher as the oracle.
#define Rgemm goal10_reference_gemm
#include "../../../../../../mpblas/reference/Rgemm.cpp"
#define main regular_gemm_test_main
#include "../../../../../../mpblas/optimized/mpfr/cuda/test/Rgemm_cuda_test.cpp"
#undef main
#undef Rgemm
int main() {
    if(mplapack_mpfr_cuda::device_count()==0) return 77;
    setenv("MPLAPACK_MPFR_CUDA","1",1);
    setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);
    gmp_randinit_default(rng);gmp_randseed_ui(rng,42);
    const mplapackint shapes[][3]={{7,7,7},{63,63,63},{64,64,64},{65,65,65},{200,17,9},{1000,17,7}};
    for(int prec:{256,512,1024}) {
        mpfrxx::set_default_precision_bits(prec);
        // m=0, n=0, and k=0 with beta=1 are true quick returns: false,
        // with every supplied output (including negative zero) untouched.
        for(int zero_dim=0;zero_dim<3;++zero_dim) {
            std::vector<mpfr_class> a(9,1),b(9,2),c(9,3);
            mpfr_set_zero(c[0].get_mpfr_t(),-1);auto original=c;
            bool used=Rgemm_mpfr_cuda("N","N",zero_dim==0?0:3,zero_dim==1?0:3,zero_dim==2?0:3,
                mpfr_class(1),a.data(),3,b.data(),3,mpfr_class(1),c.data(),3);
            if(used) ++failures;
            for(int i=0;i<9;++i) { ++checked;if(!same_bits(original[i],c[i])) ++failures; }
        }
        for(auto ta:{"N","T","C"}) for(auto tb:{"N","T","C"})
            for(auto &s:shapes) for(double beta:{0.,1.,.75,-1.25}) for(int alpha=0;alpha<3;++alpha)
            {
                long old=signed_zero_diffs;
                check_eligible(prec,ta,tb,s[0],s[1],s[2],beta,alpha);
                if(signed_zero_diffs!=old) std::printf("legacy zero signs: %d %s%s %ld,%ld,%ld beta=%g alpha=%d count=%ld\n",prec,ta,tb,(long)s[0],(long)s[1],(long)s[2],beta,alpha,signed_zero_diffs-old);
            }
    }
    gmp_randclear(rng);
    std::printf("%ld elements, %ld zero-sign differences, %d failures\n",checked,signed_zero_diffs,failures);
    return failures?1:0;
}
