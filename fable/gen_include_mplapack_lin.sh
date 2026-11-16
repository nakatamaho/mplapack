#!/bin/bash
fable_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fable_repo_root="$(cd "${fable_script_dir}/.." && pwd)"
. "${fable_script_dir}/clang_format_common.sh"

cd "${fable_repo_root}/mplapack/test/lin/common"

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp`
for filename in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' "${filename}" | sed -E 's/^typename[[:space:]]*:[[:space:]]*//' > "${filename%.*}.hpp"
done

printf "REAL Rlamch(const char *cmach);" > Rlamch.hpp
printf "INTEGER iMlaenv2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);" > iMlaenv2stage.hpp
printf "INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);" > iMlaenv.hpp

cat *hpp \
  | grep -v abs1 \
  | grep -vE '^[[:space:]]*-[[:space:]]+' \
  | grep -v main \
  | grep -v program_ \
  | LC_ALL=C sort | uniq > header_all

rm *hpp

# Backend names and their REAL/COMPLEX types come from backends.txt.
. "${fable_repo_root}/misc/backends.sh"
for mplib in $(mplapack_backend_names); do
    real_type=$(mplapack_backend_field "$mplib" real)
    complex_type=$(mplapack_backend_field "$mplib" complex)
    case "$mplib" in
        gmp) cat header_all | grep -v mpfr ;;
        mpfr) cat header_all | grep -v gmp ;;
        *) cat header_all | grep -v gmp | grep -v mpfr ;;
    esac > mplapack_lin_${mplib}.h
    sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h
    sed -i -e "s/COMPLEX/${complex_type}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/REAL/${real_type}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h
    sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h

    fable_clang_format_stdout mplapack_lin_${mplib}.h | fable_sort_prototypes "$mplib" > l ; mv l mplapack_lin_${mplib}.h
    cat "${fable_repo_root}/mplapack/test/lin/common/mplapack_lin_${mplib}.h.in" mplapack_lin_${mplib}.h > "${fable_repo_root}/include/mplapack_lin_${mplib}.h"
    rm mplapack_lin_${mplib}.h
    echo "#endif" >> "${fable_repo_root}/include/mplapack_lin_${mplib}.h"

done

mv header_all mplapack_lin_generic.h

for f in mplapack_lin_generic.h; do
fable_clang_format_inplace "$f"
done
