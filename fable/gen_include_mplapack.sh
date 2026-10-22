#!/bin/bash
fable_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fable_repo_root="$(cd "${fable_script_dir}/.." && pwd)"
. "${fable_script_dir}/clang_format_common.sh"

cd "${fable_repo_root}/mplapack/reference"

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp | grep -v Rlamch`
for filename in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${filename} |  tr ':' ' ' | sed -e 's/^typename //' >  ${filename%.*}.hpp
done

printf "REAL Rlamch(const char *cmach);\nREAL Rlamc3(REAL a, REAL b);\n" > Rlamch.hpp

cat *hpp | LC_ALL=C sort | grep -v iseed_is_all_minus_one | grep -v rlaruv_print_nondet_banner_once | grep -v fixed_point | grep -v Mxerbla | grep -v abs1 | grep -v abs2 | grep -v ___mplapack_ | grep -v nondeterministic | grep -v advance_iseed | grep -v iseed_to_seed64 | grep -v abssq | grep -v ___random_mplapack_gmp > header_all

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
    esac > mplapack_${mplib}.h
    sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h
    sed -i -e "s/COMPLEX/${complex_type}/g" mplapack_${mplib}.h
    sed -i -e "s/REAL/${real_type}/g" mplapack_${mplib}.h
    sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h
    sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
    case "$mplib" in
        gmp) printf "void mplapack_gmp_initialize(void);" >> mplapack_${mplib}.h ;;
        mpfr)
            printf "void ___mplapack_mpfr_initialize(void);" >> mplapack_${mplib}.h
            printf "void mplapack_mpfr_finalize(void);" >> mplapack_${mplib}.h
            ;;
    esac

    fable_clang_format_stdout mplapack_${mplib}.h | fable_sort_prototypes "$mplib" > l ; mv l mplapack_${mplib}.h
    {
        cat "${fable_repo_root}/mplapack/reference/mplapack_${mplib}.h.in"
        # Preserve the established MPFR header separator after the template.
        if [ "$mplib" = mpfr ]; then printf '\n'; fi
        cat mplapack_${mplib}.h
    } > "${fable_repo_root}/include/mplapack_${mplib}.h"
    rm mplapack_${mplib}.h
    echo "#endif" >> "${fable_repo_root}/include/mplapack_${mplib}.h"

done

mv header_all "${fable_repo_root}/mplapack/reference/mplapack_generic.h"

for f in mplapack_generic.h; do
fable_clang_format_inplace "$f"
done

{
cat <<'EOF'
/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_generic.h,v 1.18 2010/08/07 03:15:46 nakatamaho Exp $
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#ifndef MPLAPACK_GENERIC_H
#define MPLAPACK_GENERIC_H

/* MPLAPACK prototypes */

EOF
LC_ALL=C sort "${fable_repo_root}/mplapack/reference/mplapack_generic.h"
printf "\n#endif\n"
} > "${fable_repo_root}/include/mplapack_generic.h"
