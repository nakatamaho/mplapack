#!/bin/bash
# MPFR 512 re-measurement at HEAD (after goals 12/13): Rgemm NN/NT/TN/TT kernel-only
# sweep (as ../gemm-sweep, all forms to 4096) and Rgetrf / Rpotrf L / Rpotrf U
# through the GPU entries (Rgetrf_mpfr_cuda, Rpotrf_mpfr_cuda; whole entry
# incl. packing and transfers, warmup 1, best of 3 for n <= 2000, 1 run above).
set -u
R=/home/docker/mplapack-acceralation
O=$(cd "$(dirname "$0")" && pwd)
export LD_LIBRARY_PATH=$R/external/i/QD/lib:$R/external/i/MPFR/lib:$R/external/i/GMP/lib:$R/external/i/MPC/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export MPFRXX_DEFAULT_PRECISION_BITS=512
# every n through the GPU entry (the default threshold rejects small n and the entry returns false)
export MPLAPACK_MPFR_LAPACK_CUDA_MIN=1
# FACTOR_ONLY=1 skips the GEMM sweeps
log() { echo "$(date +%FT%T) $* load=$(cut -d' ' -f1-3 /proc/loadavg) gpu_procs=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" >> $O/run.log; }
git -C $R rev-parse HEAD > $O/commit.txt
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-start.txt
[ "${FACTOR_ONLY:-0}" = 1 ] || for f in NN NT TN TT; do
  log "gemm $f start"; $R/benchmark/run_gemm_sweep.sh $O mpfr512 4096 $f; log "gemm $f end rc=$?"
done
F=$R/benchmark/Rfactor.mpfr_cuda_total
log "getrf start"; $F --entry --warmup 1 --max 4096 > $O/getrf.csv 2> $O/getrf.err; log "getrf end rc=$?"
log "potrf L start"; $F --potrf --entry --warmup 1 --max 4096 > $O/potrf_L.csv 2> $O/potrf_L.err; log "potrf L end rc=$?"
log "potrf U start"; $F --potrf --upper --entry --warmup 1 --max 4096 > $O/potrf_U.csv 2> $O/potrf_U.err; log "potrf U end rc=$?"
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-end.txt
log done
