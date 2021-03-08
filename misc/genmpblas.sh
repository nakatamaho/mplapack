rm mpblas_gmp.h
rm mpblas_mpfr.h
rm mpblas_qd.h
rm mpblas_dd.h

cp mpblas_generic.h mpblas_mpfr.h 
sed -i.bak "s/MPBLAS_GENERIC_H/MPBLAS_MPFR_H/g" mpblas_mpfr.h
sed -i.bak "s/INTEGER/mplapackint/g" mpblas_mpfr.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mpblas_mpfr.h
sed -i.bak "s/REAL/mpreal/g" mpblas_mpfr.h
sed -i.bak "s/COMPLEX/mpcomplex/g" mpblas_mpfr.h
sed -i.bak "s/Mlsame/Mlsame_mpfr/g" mpblas_mpfr.h
sed -i.bak "s/Mxerbla/Mxerbla_mpfr/g" mpblas_mpfr.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mplapack_config.h>\\
#include <mpreal.h>\\
#include <mpcomplex.h>\\
#include <mutils_mpfr.h>/" mpblas_mpfr.h

cp mpblas_generic.h mpblas_gmp.h 
sed -i.bak "s/MPBLAS_GENERIC_H/MPBLAS_GMP_H/g" mpblas_gmp.h
sed -i.bak "s/INTEGER/mplapackint/g" mpblas_gmp.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mpblas_gmp.h
sed -i.bak "s/REAL/mpf_class/g" mpblas_gmp.h
sed -i.bak "s/COMPLEX/mpc_class/g" mpblas_gmp.h
sed -i.bak "s/Mlsame/Mlsame_gmp/g" mpblas_gmp.h
sed -i.bak "s/Mxerbla/Mxerbla_gmp/g" mpblas_gmp.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mplapack_config.h>\\
#include <gmpxx.h>\\
#include <mpc_class.h>\\
#include <mutils_gmp.h>/" mpblas_gmp.h

cp mpblas_generic.h mpblas_qd.h 
sed -i.bak "s/MPBLAS_GENERIC_H/MPBLAS_QD_H/g" mpblas_qd.h
sed -i.bak "s/INTEGER/mplapackint/g" mpblas_qd.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mpblas_qd.h
sed -i.bak "s/REAL/qd_real/g" mpblas_qd.h
sed -i.bak "s/COMPLEX/qd_complex/g" mpblas_qd.h
sed -i.bak "s/Mlsame/Mlsame_qd/g" mpblas_qd.h
sed -i.bak "s/Mxerbla/Mxerbla_qd/g" mpblas_qd.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mplapack_config.h>\\
#include <qd\/qd_real.h>\\
#include <qd_complex.h>\\
#include <mutils_qd.h>/" mpblas_qd.h

cp mpblas_generic.h mpblas_dd.h 
sed -i.bak "s/MPBLAS_GENERIC_H/MPBLAS_DD_H/g" mpblas_dd.h
sed -i.bak "s/INTEGER/mplapackint/g" mpblas_dd.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mpblas_dd.h
sed -i.bak "s/REAL/dd_real/g" mpblas_dd.h
sed -i.bak "s/COMPLEX/dd_complex/g" mpblas_dd.h
sed -i.bak "s/Mlsame/Mlsame_dd/g" mpblas_dd.h
sed -i.bak "s/Mxerbla/Mxerbla_dd/g" mpblas_dd.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mplapack_config.h>\\
#include <qd\/dd_real.h>\\
#include <dd_complex.h>\\
#include <mutils_dd.h>/" mpblas_dd.h

cp mpblas_generic.h mpblas_double.h 
sed -i.bak "s/MPBLAS_GENERIC_H/MPBLAS_DOUBLE_H/g" mpblas_double.h
sed -i.bak "s/INTEGER/mplapackint/g" mpblas_double.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mpblas_double.h
sed -i.bak "s/REAL/double/g" mpblas_double.h
sed -i.bak "s/COMPLEX/std::complex<double>/g" mpblas_double.h
sed -i.bak "s/Mlsame/Mlsame_double/g" mpblas_double.h
sed -i.bak "s/Mxerbla/Mxerbla_double/g" mpblas_double.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mplapack_config.h>\\
#include <mutils_double.h>/" mpblas_double.h

