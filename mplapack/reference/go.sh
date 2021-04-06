bash -x ~/mplapack/misc/conv_all_lapack.sh 
bash -x ~/mplapack/misc/make_include_lapack.sh
cp hand/*.cpp  hand/*.hpp .
#patch -p0 < patch