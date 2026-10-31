#!/bin/sh
# Full solver cases at every precision, isolated so host emulation needs no GPU
# and different MPFR/MPC defaults cannot leak between precision jobs.
set -u
work=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-solver-host.XXXXXX") || exit 1
trap 'rm -rf "$work"' EXIT HUP INT TERM
pids=""
for precision in 256 512 768 1024 2048 200; do
    ./solver_gpu_host_test_mpfr "$precision" > "$work/$precision.log" 2>&1 &
    pids="$pids $!"
done
status=0
for pid in $pids; do
    wait "$pid" || status=1
done
for precision in 256 512 768 1024 2048 200; do
    cat "$work/$precision.log"
done
exit "$status"
