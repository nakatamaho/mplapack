#!/bin/sh

#remove everything

rm -f /home/docker/mplapack/include/mplapack_{mpfr,gmp,qd,dd,_Float128,_Float64x,double}.h
rm -f /home/docker/mplapack/mpblas/mpblas_generic.h
rm -f /home/docker/mplapack/mplapack/mplapack_generic.h
rm -rf /home/docker/mplapack/external/lapack/work

# convert BLAS, LAPACK subroutine names to MPLAPACK function name
# typically dgemm -> Rgemm, zgemm -> Cgemm
# manual mappings can be found in MANUAL_MAPPINGS of gen_mplapack_name_map.sh
cd /home/docker/mplapack/external/lapack/ ; make extract
bash /home/docker/mplapack/fable/gen_mplapack_name_map.sh

# convert BLAS and LAPACK: 1st pass
cd /home/docker/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC
bash /home/docker/mplapack/fable/convert_blas_all.sh ; mv *cpp /home/docker/mplapack/mpblas/reference/

cd /home/docker/mplapack/external/lapack/work/internal/lapack-3.9.1/LAPACK/SRC
bash /home/docker/mplapack/fable/convert_lapack_all.sh ; mv *cpp /home/docker/mplapack/mplapack/mpblas/reference/

### make prototype header
bash /home/docker/mplapack/fable/make_include_blas.sh
bash /home/docker/mplapack/fable/make_include_lapack.sh

#2nd pass
cd /home/docker/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC
bash /home/docker/mplapack/fable/convert_blas_all.sh ; mv *cpp /home/docker/mplapack/mpblas/reference/

cd /home/docker/mplapack/external/lapack/work/internal/lapack-3.9.1/LAPACK/SRC
bash /home/docker/mplapack/fable/convert_lapack_all.sh ; mv *cpp /home/docker/mplapack/mplapack/mpblas/reference/

### make prototype header
bash /home/docker/mplapack/fable/make_include_blas.sh
bash /home/docker/mplapack/fable/make_include_lapack.sh
