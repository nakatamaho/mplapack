#!/bin/bash

echo "fix dylibs for dir: $1, prefix $2"
DIR=${1:-}
PREFIX=${2:-}
if [ ! -d "$DIR" ]; then
    echo "skip missing directory: $DIR"
    exit 0
fi
cd "$DIR" || exit 1
for _file in *; do
    [ -f "$_file" ] || continue
    echo "$_file"
    otool -L "$_file" | grep "external/mpc/../i/MPC/lib/libmpc"
    if [ $? = 0 ]; then
        ORG=`otool -L "$_file"  | grep MPC | awk '{print $1}'` ;\
        _FIX=`otool -L "$_file" | grep MPC | awk '{print $1}' | sed 's|/| |g' | awk '{print $NF}'` ; FIX="$PREFIX/lib/$_FIX" ;\
        echo "install_name_tool -change $ORG $FIX $_file" ;\
        install_name_tool -change "$ORG" "$FIX" "$_file"
    fi
    otool -L "$_file" | grep "external/mpfr/../i/MPFR/lib/libmpfr"
    if [ $? = 0 ]; then
        ORG=`otool -L "$_file"  | grep MPFR | awk '{print $1}'` ;\
        _FIX=`otool -L "$_file" | grep MPFR | awk '{print $1}' | sed 's|/| |g' | awk '{print $NF}'` ; FIX="$PREFIX/lib/$_FIX" ;\
        echo "install_name_tool -change $ORG $FIX $_file" ;\
        install_name_tool -change "$ORG" "$FIX" "$_file"
    fi
    otool -L "$_file" | grep "external/gmp/../i/GMP/lib/libgmpxx"
    if [ $? = 0 ]; then
        ORG=`otool -L "$_file"  | grep GMP | grep gmpxx | awk '{print $1}'` ;\
        _FIX=`otool -L "$_file" | grep GMP | grep gmpxx | awk '{print $1}' | sed 's|/| |g' | awk '{print $NF}'` ; FIX="$PREFIX/lib/$_FIX" ;\
        echo "install_name_tool -change $ORG $FIX $_file" ;\
        install_name_tool -change "$ORG" "$FIX" "$_file"
    fi
    otool -L "$_file" | grep "external/gmp/../i/GMP/lib/libgmp"
    if [ $? = 0 ]; then
        ORG=`otool -L "$_file"  | grep GMP | grep -v gmpxx | awk '{print $1}'` ;\
        _FIX=`otool -L "$_file" | grep GMP | grep -v gmpxx | awk '{print $1}' | sed 's|/| |g' | awk '{print $NF}'` ; FIX="$PREFIX/lib/$_FIX" ;\
        echo "install_name_tool -change $ORG $FIX $_file" ;\
        install_name_tool -change "$ORG" "$FIX" "$_file"
    fi
    otool -L "$_file" | grep "external/qd/../i/QD/lib/libqd"
    if [ $? = 0 ]; then
        ORG=`otool -L "$_file"  | grep QD | awk '{print $1}'` ;\
        _FIX=`otool -L "$_file" | grep QD | awk '{print $1}' | sed 's|/| |g' | awk '{print $NF}'` ; FIX="$PREFIX/lib/$_FIX" ;\
        echo "install_name_tool -change $ORG $FIX $_file" ;\
        install_name_tool -change "$ORG" "$FIX" "$_file"
    fi
done
