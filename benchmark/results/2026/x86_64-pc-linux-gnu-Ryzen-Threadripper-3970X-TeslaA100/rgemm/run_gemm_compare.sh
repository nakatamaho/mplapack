#!/bin/bash
# usage: run_gemm_compare.sh dd|td|qd|mpfr [BINDIR=./bin]   (stdout -> gemm_compare_<type>.txt)
# Binaries from build.sh.  Matrices for the accuracy comparison go to ./dump-<type>.
t=$1; X=${2:-./bin}/bench_$t; D=./dump-$t; mkdir -p $D
case $t in
  dd) P=MPLAPACK_DD; sizes="1024 2048"; acc=1024 ;;
  td) P=MPLAPACK_TD; sizes="1024"; acc=1024 ;;
  qd) P=MPLAPACK_QD; sizes="512 1024"; acc=512 ;;
  mpfr) P=MPLAPACK_MPFR; sizes="512 1024"; acc=512 ;;
esac
if [ $t = mpfr ]; then
  variants="cpu default wino256 wino512 runtime"
else
  variants="cpu default naive wino256 wino512 ozaki"
fi
envof() {
  case $1 in
    cpu) if [ $t = mpfr ]; then echo "${P}_CUDA=0"; else echo "${P}_CUDA=0 ${P}_BLAS_CUDA=0"; fi ;;
    default) echo "X=1" ;;
    naive) echo "${P}_CUDA_KERNEL=naive" ;;
    wino256) echo "${P}_CUDA_WINOGRAD_CUTOFF=256" ;;
    wino512) echo "${P}_CUDA_WINOGRAD_CUTOFF=512" ;;
    ozaki) echo "${P}_CUDA_OZAKI=1" ;;
    runtime) echo "${P}_CUDA_FORCE_RUNTIME=1" ;;
  esac
}
for n in $sizes; do for R in R C; do for v in $variants; do
  dump=""; [ $n = $acc ] && dump=$D/${R}_${n}_$v.bin
  out=$(env $(envof $v) $X $R $n 3 $dump 2>&1 | tail -1)
  echo "$t $R n=$n $v : $out"
done; done; done
echo "== accuracy (n=$acc, vs cpu)"
for R in R C; do for v in $variants; do [ $v = cpu ] && continue
  echo "$t $R n=$acc $v : $($X cmp $R $acc $D/${R}_${acc}_cpu.bin $D/${R}_${acc}_$v.bin)"
done; done
echo "== verbose (default, n=${sizes%% *})"
for R in R C; do env ${P}_CUDA_VERBOSE=1 ${P}_BLAS_CUDA_VERBOSE=1 $X $R ${sizes%% *} 1 2>&1 | sed "s/^/$t $R verbose: /"; done
rm -f $D/*.bin
