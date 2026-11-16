#!/bin/bash
# GPU kernel-only Rgemm(TRANSA,TRANSB,n,n,n) sweep with benchmark/Rkernel_bench.<type>_cuda_kernel
# (data on the GPU, CUDA events, no transfers).  n = multiples of 7, multiples
# of 31 and powers of two up to NMAX; best of 3 after one warm-up for
# n <= 2000, one run without warm-up above (as the 2026-10 MPFR sweep).
# Output: OUTDIR/<TYPE>_<FORM>.csv (n,seconds,gflops).
#   usage: run_gemm_sweep.sh OUTDIR TYPE NMAX [FORM=NN]   TYPE = dd|td|qd|mpfr<bits>, FORM = NN|NT|TN|TT
set -eu
out=$1; t=$2; nmax=$3; form=${4:-NN}
R=$(cd "$(dirname "$0")/.." && pwd)
export LD_LIBRARY_PATH=$R/external/i/QD/lib:$R/external/i/MPFR/lib:$R/external/i/GMP/lib:$R/external/i/MPC/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
case $t in mpfr*) bin=$R/benchmark/Rkernel_bench.mpfr_cuda_kernel; export MPFRXX_DEFAULT_PRECISION_BITS=${t#mpfr} ;;
           *) bin=$R/benchmark/Rkernel_bench.${t}_cuda_kernel ;; esac
sizes=$( (seq 7 7 $nmax; seq 31 31 $nmax; p=1; while [ $p -le $nmax ]; do echo $p; p=$((p*2)); done) | sort -nu)
small=$(echo "$sizes" | awk '$1<=2000' | paste -sd,)
large=$(echo "$sizes" | awk '$1>2000' | paste -sd,)
mkdir -p $out
{ echo "n,seconds,gflops"
  $bin gemm $form $small --label $t --reps 3 --warmup 1 | grep -v '^#' | cut -d, -f4,6,7
  if [ -n "$large" ]; then $bin gemm $form $large --label $t --reps 1 --warmup 0 | grep -v '^#' | cut -d, -f4,6,7; fi
} > $out/${t}_$form.csv.tmp && mv $out/${t}_$form.csv.tmp $out/${t}_$form.csv
