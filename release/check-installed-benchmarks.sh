#!/bin/bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <prefix>" >&2
    exit 2
fi

prefix="$1"
case "${prefix}" in
    /*) ;;
    *)
        echo "ERROR: prefix must be an absolute path: ${prefix}" >&2
        exit 1
        ;;
esac

cpu_detection="${prefix}/misc/cpu_os_detection.sh"
if [ ! -f "${cpu_detection}" ]; then
    echo "ERROR: installed CPU detection helper not found: ${cpu_detection}" >&2
    exit 1
fi
if [ ! -x "${cpu_detection}" ]; then
    echo "ERROR: installed CPU detection helper is not executable: ${cpu_detection}" >&2
    exit 1
fi

echo "=== Checking installed CPU detection helper: ${cpu_detection} ==="
detection_output="$(sh "${cpu_detection}")"
echo "cpu_os_detection.sh: ${detection_output}"
case "${detection_output}" in
    *"|"*"|"*) ;;
    *)
        echo "ERROR: unexpected cpu_os_detection.sh output: ${detection_output}" >&2
        exit 1
        ;;
esac

benchmark_dirs=()
while IFS= read -r benchmark_dir; do
    benchmark_dirs+=("${benchmark_dir}")
done < <(find "${prefix}/lib" -path '*/mplapack/benchmark' -type d -print 2>/dev/null | sort)

if [ "${#benchmark_dirs[@]}" -eq 0 ]; then
    echo "SKIP: no installed benchmark directory found under ${prefix}/lib"
    exit 0
fi

placeholders='%%PREFIX%%|%%CPU_OS_DETECTION_SH%%|%%CC_TAG%%'
representative_bins=(
    daxpy_ref
    dgemm_ref
    Raxpy.double
    Rgemm.double
    Rgetrf.double
    Rpotrf.double
)

for benchmark_dir in "${benchmark_dirs[@]}"; do
    echo "=== Checking installed benchmark scripts: ${benchmark_dir} ==="

    shopt -s nullglob
    go_scripts=("${benchmark_dir}"/go.*.sh)
    shopt -u nullglob
    if [ "${#go_scripts[@]}" -eq 0 ]; then
        echo "ERROR: no installed go.*.sh scripts found in ${benchmark_dir}" >&2
        exit 1
    fi

    if grep -E "${placeholders}" "${go_scripts[@]}"; then
        echo "ERROR: unresolved benchmark placeholders found in ${benchmark_dir}" >&2
        exit 1
    fi

    for go_script in "${go_scripts[@]}"; do
        sh -n "${go_script}"
    done

    missing_representative=0
    for bin in "${representative_bins[@]}"; do
        if [ -x "${benchmark_dir}/${bin}" ]; then
            echo "found benchmark binary: ${benchmark_dir}/${bin}"
        elif [ -x "${benchmark_dir}/${bin}.exe" ]; then
            echo "found benchmark binary: ${benchmark_dir}/${bin}.exe"
        else
            echo "missing representative benchmark binary: ${benchmark_dir}/${bin}"
            missing_representative=1
        fi
    done
    if [ "${missing_representative}" -ne 0 ]; then
        echo "ERROR: missing representative benchmark binaries in ${benchmark_dir}" >&2
        exit 1
    fi
done
