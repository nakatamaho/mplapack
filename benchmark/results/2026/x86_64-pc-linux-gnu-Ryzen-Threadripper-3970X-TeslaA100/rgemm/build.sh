#!/bin/bash
# Build the Rgemm comparison (benchmark/Rgemm_gpu_compare.cpp -> bench_<type>)
# and the MPFR sweep (benchmark/Rgemm_mpfr_sweep.cpp -> mpfr_sweep) against an
# autotools build of this tree with --enable-{dd,td,qd,mpfr}-cuda.
#   usage: build.sh [OUTDIR=./bin]
set -eu
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
O=${1:-./bin}; mkdir -p "$O"
T=$R/mpblas/optimized
I="-I$R/include -I$R/external/i/QD/include -I$R/external/i/GMPFRXX_MKII/include -I$R/external/i/MPFR/include -I$R/external/i/GMP/include -I$R/external/i/MPC/include -I/usr/local/cuda/include"
CXXF="-O3 -std=gnu++17 -ffp-contract=off -fopenmp -Wl,--disable-new-dtags"
RP=""; for d in $T/dd/cuda/.libs $T/td/cuda/.libs $T/qd/cuda/.libs $T/mpfr/cuda/.libs $R/external/i/QD/lib $R/external/i/MPFR/lib $R/external/i/GMP/lib $R/external/i/MPC/lib /usr/local/cuda/lib64; do RP="$RP -Wl,-rpath,$d"; done
LC="-L$R/external/i/QD/lib -L$R/external/i/MPFR/lib -L$R/external/i/GMP/lib -L$R/external/i/MPC/lib -L/usr/local/cuda/lib64 -lcudart -lgomp $RP"
for b in dd td qd; do
  g++ $CXXF -DMPLAPACK_BUILD_WITH_${b^^} $I $R/benchmark/Rgemm_gpu_compare.cpp -o $O/bench_$b \
      -L$T/$b/cuda/.libs -lmplapack_${b}_opt_cuda -Wl,--no-as-needed -lqd -Wl,--as-needed $LC
done
g++ $CXXF -DMPLAPACK_BUILD_WITH_MPFR $I $R/benchmark/Rgemm_gpu_compare.cpp -o $O/bench_mpfr \
    -L$T/mpfr/cuda/.libs -lmplapack_mpfr_opt_cuda -lmpc -lmpfr -lgmp $LC
g++ $CXXF -DMPLAPACK_BUILD_WITH_MPFR $I $R/benchmark/Rgemm_mpfr_sweep.cpp -o $O/mpfr_sweep \
    -L$T/mpfr/cuda/.libs -lmplapack_mpfr_opt_cuda -lmpc -lmpfr -lgmp $LC
echo "built in $O"
