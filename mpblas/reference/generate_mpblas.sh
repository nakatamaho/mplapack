bash -x ~/mplapack/misc/conv_all_blas.sh
patch -p3 < patch
sed -i -e 's/poINTEGER/point/g' *cpp
bash -x ~/mplapack/misc/make_include_blas.sh
bash -x ~/mplapack/misc/indent.sh
