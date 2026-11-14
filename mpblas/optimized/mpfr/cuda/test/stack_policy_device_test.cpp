// Autotools device-side counterparts of the seven CTest restore families.
#include "stack_policy_test.h"
using namespace mpfr_stack_test;
int main(int argc,char **argv) {
    int devices=0; if(cudaGetDeviceCount(&devices) || !devices) return 77;
    if(argc!=2) return 2;
    setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    setenv("MPLAPACK_MPFR_LAPACK_SOLVE_CUDA_MIN","0",1);
    setenv("MPLAPACK_MPFR_LAPACK_DRIVER_CUDA_MIN","0",1);
    setenv("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN","0",1);
    std::string family=argv[1];
    long p=family=="Rgemm"?200:512;
    mpfrxx::set_default_precision_bits(p);mpfrxx::set_default_mpc_precision_bits(p);
    if(family=="Cfactor" || family=="Cblas") {
        policies(family.c_str(),family=="Cfactor"?factor_bytes(p,true):size_t(16384)*((p+511)/512),[&] {
            result<mpc_class> r; r.a.assign(9,mpc_class(.125));r.b.assign(3,mpc_class(.25));r.piv.assign(3,-7);
            for(int i=0;i<3;++i) r.a[i+i*3]=mpc_class(4);
            if(family=="Cfactor") r.computed=Cgetrf_mpfr_cuda(3,3,r.a.data(),3,r.piv.data(),r.info);
            else r.computed=Ctrsm_mpfr_cuda("L","U","N","N",3,1,mpc_class(1),r.a.data(),3,r.b.data(),3);
            return r;
        });
    } else if(family=="solver" || family=="Rfactor" || family=="inverse" || family=="Rgemm" || family=="blas3") {
        size_t bytes=family=="inverse"?inverse_bytes(p):family=="Rgemm"?65536:factor_bytes(p,false);
        policies(family.c_str(),bytes,[&] {
            result<mpfr_class> r;r.a.assign(9,mpfr_class(.125));r.b.assign(9,mpfr_class(.25));r.work.assign(9,mpfr_class(.5));r.piv.assign(3,-7);
            for(int i=0;i<3;++i) r.a[i+i*3]=4;
            if(family=="Rfactor") r.computed=Rgetrf_mpfr_cuda(3,3,r.a.data(),3,r.piv.data(),r.info);
            if(family=="solver") r.computed=Rgesv_mpfr_cuda(3,1,r.a.data(),3,r.piv.data(),r.b.data(),3,r.info);
            if(family=="inverse") r.computed=Rtrtri_mpfr_cuda("U","N",3,r.a.data(),3,r.info);
            if(family=="Rgemm") r.computed=Rgemm_mpfr_cuda("N","N",3,3,3,mpfr_class(1),r.a.data(),3,r.b.data(),3,mpfr_class(.25),r.work.data(),3);
            if(family=="blas3") r.computed=Rtrsm_mpfr_cuda("L","U","N","N",3,1,mpfr_class(1),r.a.data(),3,r.b.data(),3);
            return r;
        });
    } else return 2;
    return 0;
}
