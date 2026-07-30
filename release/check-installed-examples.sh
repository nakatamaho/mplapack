#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: $0 <prefix> <makefile-name[,fallback...]> [jobs] [make-var ...]" >&2
    exit 2
fi

prefix="$1"
makefile_spec="$2"
jobs="${3:-}"
shift 2
if [ "$#" -gt 0 ]; then
    shift
fi
make_args=("$@")
IFS=',' read -r -a makefile_names <<< "${makefile_spec}"

case "${prefix}" in
    /*) ;;
    *)
        echo "ERROR: prefix must be an absolute path: ${prefix}" >&2
        exit 1
        ;;
esac

example_dirs=(
    "${prefix}/share/examples/mpblas"
    "${prefix}/share/examples/mplapack/00_GeneralLinearEquations"
    "${prefix}/share/examples/mplapack/01_PositiveDefiniteLinearEquations"
    "${prefix}/share/examples/mplapack/02_LeastSquares"
    "${prefix}/share/examples/mplapack/03_SymmetricEigenproblems"
    "${prefix}/share/examples/mplapack/04_NonsymmetricEigenproblems"
    "${prefix}/share/examples/mplapack/05_SingularValueDecomposition"
    "${prefix}/share/examples/mplapack/06_SymmetricIndefiniteLinearEquations"
    "${prefix}/share/examples/mplapack/07_GeneralizedSymmetricDefiniteEigenproblems"
    "${prefix}/share/examples/mplapack/08_GeneralizedNonsymmetricEigenproblems"
    "${prefix}/share/examples/mplapack/09_GeneralizedSingularValueDecomposition"
    "${prefix}/share/examples/mplapack/90_PrecisionComparison"
)

for example_dir in "${example_dirs[@]}"; do
    if [ ! -d "${example_dir}" ]; then
        echo "ERROR: installed example directory not found: ${example_dir}" >&2
        exit 1
    fi
    selected_makefile=""
    for candidate in "${makefile_names[@]}"; do
        if [ -f "${example_dir}/${candidate}" ]; then
            selected_makefile="${candidate}"
            break
        fi
    done
    if [ -z "${selected_makefile}" ]; then
        echo "ERROR: installed example Makefile not found in ${example_dir}: ${makefile_spec}" >&2
        exit 1
    fi

    echo "=== Checking installed examples: ${example_dir} (${selected_makefile}) ==="
    make_cmd=(make -f "${selected_makefile}")
    if [ -n "${jobs}" ]; then
        make_cmd+=("-j${jobs}")
    fi
    (
        cd "${example_dir}"
        "${make_cmd[@]}" "${make_args[@]}" clean
        "${make_cmd[@]}" "${make_args[@]}"
    )
done
