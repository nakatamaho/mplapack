#!/bin/bash
# MPFR 512 on the V100, following ../../x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/mpfr512-sweep-2:
# Rgemm NN/NT/TN/TT kernel-only sweeps (as ../../*TeslaA100/gemm-sweep), then
# Rgetrf / Rpotrf L / Rpotrf U through the GPU entries (whole entry incl. packing
# and transfers, warmup 1, best of 3 for n <= 2000, 1 run above), all to 4096.
set -u
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
O=$(cd "$(dirname "$0")" && pwd)
# HPC SDK 24.11: CUDA 12.6 runtime and cuBLAS (math_libs) live in separate trees
H=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11
export LD_LIBRARY_PATH=$R/external/i/QD/lib:$R/external/i/MPFR/lib:$R/external/i/GMP/lib:$R/external/i/MPC/lib:$H/cuda/12.6/lib64:$H/math_libs/12.6/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export MPFRXX_DEFAULT_PRECISION_BITS=512
# every n through the GPU entry (the default threshold rejects small n and the entry returns false)
export MPLAPACK_MPFR_LAPACK_CUDA_MIN=1
# GEMM_ONLY=1 skips the factorizations, FACTOR_ONLY=1 the GEMM sweeps
log() { echo "$(date +%FT%T) $* load=$(cut -d' ' -f1-3 /proc/loadavg) gpu_procs=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" >> $O/run.log; }
git -C $R rev-parse HEAD > $O/commit.txt
git -C $R diff --stat >> $O/commit.txt
nvidia-smi --query-gpu=index,name,driver_version,clocks.max.sm,memory.total --format=csv > $O/gpu-info.txt
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-start.txt
[ "${FACTOR_ONLY:-0}" = 1 ] || for f in NN NT TN TT; do
  log "gemm $f start"; $R/benchmark/run_gemm_sweep.sh $O mpfr512 4096 $f; log "gemm $f end rc=$?"
done
if [ "${GEMM_ONLY:-0}" != 1 ]; then
  F=$R/benchmark/Rfactor.mpfr_cuda_total
  log "getrf start"; $F --entry --warmup 1 --max 4096 > $O/getrf.csv 2> $O/getrf.err; log "getrf end rc=$?"
  log "potrf L start"; $F --potrf --entry --warmup 1 --max 4096 > $O/potrf_L.csv 2> $O/potrf_L.err; log "potrf L end rc=$?"
  log "potrf U start"; $F --potrf --upper --entry --warmup 1 --max 4096 > $O/potrf_U.csv 2> $O/potrf_U.err; log "potrf U end rc=$?"
fi
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-end.txt
log done
