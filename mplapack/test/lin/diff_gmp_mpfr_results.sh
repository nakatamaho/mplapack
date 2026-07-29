#!/bin/sh

set -eu

baseline_version=${1:-2.2.1}
platform=${2:-Ryzen_Threadripper_3970X_ubuntu26_04_gcc-15_2_0}
baseline_root=${3:-results/${baseline_version}/${platform}}
current_root=${4:-results/${platform}}

status=0

for backend in gmp mpfr; do
    baseline_dir=${baseline_root}/${backend}
    current_dir=${current_root}/${backend}

    if [ ! -d "${baseline_dir}" ]; then
        printf '%s\n' "missing baseline directory: ${baseline_dir}" >&2
        status=1
        continue
    fi
    if [ ! -d "${current_dir}" ]; then
        printf '%s\n' "missing current directory: ${current_dir}" >&2
        status=1
        continue
    fi

    for baseline_file in "${baseline_dir}"/*.out; do
        [ -e "${baseline_file}" ] || continue
        name=${baseline_file##*/}
        current_file=${current_dir}/${name}
        if [ ! -f "${current_file}" ]; then
            printf '%s\n' "missing current file: ${current_file}" >&2
            status=1
            continue
        fi
        diff -u "${baseline_file}" "${current_file}" || status=1
    done

    for current_file in "${current_dir}"/*.out; do
        [ -e "${current_file}" ] || continue
        name=${current_file##*/}
        baseline_file=${baseline_dir}/${name}
        if [ ! -f "${baseline_file}" ]; then
            printf '%s\n' "missing baseline file: ${baseline_file}" >&2
            status=1
        fi
    done
done

exit "${status}"
