#!/bin/bash
fable_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fable_repo_root="$(cd "${fable_script_dir}/.." && pwd)"
. "${fable_script_dir}/clang_format_common.sh"

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

cd "${fable_repo_root}/mpblas/reference"

FILES=`ls *cpp | grep -v mplapackinit.cpp`

for _file in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${_file} |  tr ':' ' ' | sed -e 's/^typename //' > ${_file%.*}.hpp
done
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' Mxerbla.cpp |  tr ':' ' ' | sed -e 's/^typename //' > Mxerbla.hpp
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' Mlsame.cpp |  tr ':' ' ' | sed -e 's/^typename //' > Mlsame.hpp

cat *hpp | grep -v abssq > header_all
rm *hpp

# Backend names and their REAL/COMPLEX types come from backends.txt.
. "${fable_repo_root}/misc/backends.sh"
for mplib in $(mplapack_backend_names); do
    real_type=$(mplapack_backend_field "$mplib" real)
    complex_type=$(mplapack_backend_field "$mplib" complex)
    cp header_all mpblas_${mplib}.h
    sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
    sed -i -e "s/COMPLEX/${complex_type}/g" mpblas_${mplib}.h
    sed -i -e "s/REAL/${real_type}/g" mpblas_${mplib}.h
    sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
    sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h

    fable_clang_format_stdout mpblas_${mplib}.h | fable_sort_prototypes "$mplib" > l ; mv l mpblas_${mplib}.h
    cat "${fable_repo_root}/mpblas/reference/mpblas_${mplib}.h.in" mpblas_${mplib}.h > "${fable_repo_root}/include/mpblas_${mplib}.h"
    rm mpblas_${mplib}.h
    echo "#endif" >> "${fable_repo_root}/include/mpblas_${mplib}.h"

done
mv header_all mpblas_generic.h

for f in mpblas_generic.h; do
fable_clang_format_inplace "$f"
done
