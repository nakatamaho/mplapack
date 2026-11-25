#!/bin/bash
# dd/td/qd public Rgemm (whole call: host<->device transfers + kernels) on one
# and two V100s: MPLAPACK_{DD,TD,QD}_CUDA_NGPUS=1 vs 2, NN, square n.
# Best/mean of LOOPS runs after one warm-up, no reference check (-NOCHECK).
set -u
R=$(cd "$(dirname "$0")/../../../../.." && pwd)
O=$(cd "$(dirname "$0")" && pwd)
H=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11
export LD_LIBRARY_PATH=$R/external/i/QD/lib:$H/cuda/12.6/lib64:$H/math_libs/12.6/targets/x86_64-linux/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
SIZES=${SIZES:-1024 2048 4096 8192}
log() { echo "$(date +%FT%T) $* load=$(cut -d' ' -f1-3 /proc/loadavg) gpu_procs=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)" >> $O/run.log; }
git -C $R rev-parse HEAD > $O/commit.txt
git -C $R diff --stat -- mpblas mplapack benchmark/Rgemm.cpp >> $O/commit.txt
nvidia-smi --query-gpu=index,name,driver_version,clocks.max.sm,memory.total,pcie.link.gen.max,pcie.link.width.max --format=csv > $O/gpu-info.txt
nvidia-smi topo -m > $O/gpu-topo.txt
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-start.txt
for t in ${TYPES:-dd td qd}; do
  T=$(echo $t | tr a-z A-Z)
  for g in 1 2; do
    f=$O/${t}_ngpus$g.csv
    echo "m,n,k,transa,transb,loops,mean,median,min,max,stddev,cv" > $f
    for n in $SIZES; do
      log "$t ngpus=$g n=$n start"
      env MPLAPACK_${T}_CUDA_NGPUS=$g $R/benchmark/Rgemm.${t}_cuda_total -N $n -M $n -K $n \
        -TOTALSTEPS 1 -LOOPS ${LOOPS:-3} -WARMUP 1 -NOCHECK -NN -CSV -STATS | grep -E '^[0-9]' >> $f; rc=${PIPESTATUS[0]}
      log "$t ngpus=$g n=$n end rc=$rc"
    done
  done
done
nvidia-smi --query-compute-apps=pid,process_name --format=csv > $O/gpu-processes-at-end.txt
log done
