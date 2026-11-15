// Autotools-only adapters: one self-contained CUDA library, no CPU archive.
// Preinclude headers globally before placing the unchanged factor/inverse
// test sources (including their reference oracles) in a private namespace.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <cuda_runtime_api.h>
#include "stack_policy_test.h"

#if defined(MPLAPACK_STACK_DEVICE_RGEMM) || defined(MPLAPACK_STACK_DEVICE_BLAS3)
// BLAS eligibility is cached, so do not toggle it around an oracle call.
// Compile the unchanged reference bodies under internal names instead. These
// six routines only use scalar arithmetic and CPU argument-checking helpers.
namespace {
namespace stack_cpu_oracle {
#if defined(MPLAPACK_STACK_DEVICE_RGEMM)
#include "../../../../reference/Rgemm.cpp"
#else
#include "../../../../reference/Rsyrk.cpp"
#include "../../../../reference/Rsyr2k.cpp"
#include "../../../../reference/Rsymm.cpp"
#include "../../../../reference/Rtrmm.cpp"
#include "../../../../reference/Rtrsm.cpp"
#endif
}
}
#if defined(MPLAPACK_STACK_DEVICE_RGEMM)
#define Rgemm stack_cpu_oracle::Rgemm
#include "Rgemm_cuda_test.cpp"
#undef Rgemm
#else
#define Rsyrk stack_cpu_oracle::Rsyrk
#define Rsyr2k stack_cpu_oracle::Rsyr2k
#define Rsymm stack_cpu_oracle::Rsymm
#define Rtrmm stack_cpu_oracle::Rtrmm
#define Rtrsm stack_cpu_oracle::Rtrsm
#include "blas3_cuda_test.cpp"
#undef Rsyrk
#undef Rsyr2k
#undef Rsymm
#undef Rtrmm
#undef Rtrsm
#endif
#else
namespace {
void initialize_cpu_blas() {
    // Before the first cached public BLAS gate read. Private device-resident
    // LAPACK kernels do not use these host BLAS gates and must still compute.
    static const bool initialized = [] {
        const std::string maximum = std::to_string(std::numeric_limits<long>::max());
        mpfr_stack_test::require(!setenv("MPLAPACK_MPFR_CUDA_MIN_MNK",maximum.c_str(),1) &&
                                !setenv("MPLAPACK_MPFR_CUDA_L3","0",1) &&
                                !setenv("MPLAPACK_MPFR_CUDA_L12","0",1),
                                "select public BLAS CPU paths for LAPACK oracles");
        std::puts("CPU oracles: CUDA library CPU paths; direct LAPACK GPU entries remain enabled");
        return true;
    }();
    (void)initialized;
}
#if defined(MPLAPACK_STACK_DEVICE_SOLVER)
namespace stack_cpu_oracle {
// The shared library's LAPACK enable flag is read per call. Disable only that
// flag during a CPU oracle; never cache MPLAPACK_MPFR_CUDA=0, which would also
// disable later direct BLAS entries. Restore the exact caller environment.
class cpu_lapack_scope {
    const char *key = "MPLAPACK_MPFR_LAPACK_CUDA";
    bool had;
    std::string saved;
public:
    cpu_lapack_scope() : had(std::getenv(key) != nullptr),
        saved(had ? std::getenv(key) : "") {
        initialize_cpu_blas();
        mpfr_stack_test::require(!setenv(key,"0",1),"disable LAPACK GPU for CPU oracle");
    }
    ~cpu_lapack_scope() {
        mpfr_stack_test::require(!(had ? setenv(key,saved.c_str(),1) : unsetenv(key)),
                                "restore LAPACK GPU eligibility after CPU oracle");
    }
    cpu_lapack_scope(const cpu_lapack_scope &) = delete;
    cpu_lapack_scope &operator=(const cpu_lapack_scope &) = delete;
};
#define CPU_ORACLE(R) \
    template<class... Args> void R(Args&&... args) { \
        cpu_lapack_scope cpu; ::R(std::forward<Args>(args)...); \
    }
CPU_ORACLE(Rgetrf)
CPU_ORACLE(Cgetrf)
CPU_ORACLE(Rpotrf)
CPU_ORACLE(Cpotrf)
CPU_ORACLE(Rgetrs)
CPU_ORACLE(Cgetrs)
CPU_ORACLE(Rpotrs)
CPU_ORACLE(Cpotrs)
CPU_ORACLE(Rgesv)
CPU_ORACLE(Cgesv)
CPU_ORACLE(Rposv)
CPU_ORACLE(Cposv)
#undef CPU_ORACLE
}
// The solver's main comes from Rsolver_gpu_test.cpp, so select the CPU BLAS
// gates before main runs: the gates are cached on their first read, and a
// public BLAS call ahead of the first oracle would otherwise cache the GPU
// defaults and silently move the oracles' BLAS onto the GPU.
const bool solver_cpu_blas_selected = (initialize_cpu_blas(), true);
}
// Keep the solver's real main: unlike an ordinary renamed function, main is
// allowed to reach its closing brace without an explicit return statement.
// These aliases affect CPU calls only; ENTRY's token-pasted *_mpfr_cuda names
// are not expanded to the CPU aliases.
#define Rgetrf stack_cpu_oracle::Rgetrf
#define Cgetrf stack_cpu_oracle::Cgetrf
#define Rpotrf stack_cpu_oracle::Rpotrf
#define Cpotrf stack_cpu_oracle::Cpotrf
#define Rgetrs stack_cpu_oracle::Rgetrs
#define Cgetrs stack_cpu_oracle::Cgetrs
#define Rpotrs stack_cpu_oracle::Rpotrs
#define Cpotrs stack_cpu_oracle::Cpotrs
#define Rgesv stack_cpu_oracle::Rgesv
#define Cgesv stack_cpu_oracle::Cgesv
#define Rposv stack_cpu_oracle::Rposv
#define Cposv stack_cpu_oracle::Cposv
#include "../../../../../mplapack/optimized/common/test/Rsolver_gpu_test.cpp"
#undef Rgetrf
#undef Cgetrf
#undef Rpotrf
#undef Cpotrf
#undef Rgetrs
#undef Cgetrs
#undef Rpotrs
#undef Cpotrs
#undef Rgesv
#undef Cgesv
#undef Rposv
#undef Cposv
#else
namespace stack_suite {
#define main stack_suite_main
#if defined(MPLAPACK_STACK_DEVICE_RFACTOR)
#include "Rfactor_cuda_test.cpp"
#elif defined(MPLAPACK_STACK_DEVICE_CFACTOR)
#include "Cfactor_cuda_test.cpp"
#elif defined(MPLAPACK_STACK_DEVICE_INVERSE)
#include "Rinverse_cuda_test.cpp"
#else
#error Select an autotools stack-policy device test family
#endif
#undef main
}
}
int main(int argc,char **argv) {
    initialize_cpu_blas();
    return stack_suite::stack_suite_main(argc,argv);
}
#endif
#endif
