#!/bin/bash
# dd/td/qd Rgemm on the V100 with the sizes of the MPFR 512 sweep
# (../mpfr512-sweep, ../../*TeslaA100/gemm-sweep): kernel only
# (benchmark/run_gemm_sweep.sh: multiples of 7 and 31 and powers of two to
# 4096, best of 3 after one warm-up for n <= 2000, one run above), NN/NT/TN/TT.
set -u
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
O=$(cd "$(dirname "$0")" && pwd)
H=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11
export LD_LIBRARY_PATH=$H/cuda/12.6/lib64:$H/math_libs/12.6/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
log() { echo "$(date +%FT%T) $* load=$(cut -d' ' -f1-3 /proc/loadavg) gpu_procs=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" >> $O/run.log; }
git -C $R rev-parse HEAD > $O/commit.txt
git -C $R diff --stat >> $O/commit.txt
nvidia-smi --query-gpu=index,name,driver_version,clocks.max.sm,memory.total --format=csv > $O/gpu-info.txt
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-start.txt
for t in ${TYPES:-dd td qd}; do for f in NN NT TN TT; do
  log "$t $f start"; $R/benchmark/run_gemm_sweep.sh $O $t 4096 $f; log "$t $f end rc=$?"
done; done
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-end.txt
log done
