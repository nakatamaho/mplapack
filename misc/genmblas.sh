rm mblas_gmp.h
rm mblas_mpfr.h
rm mblas_qd.h
rm mblas_dd.h

cp mblas_generic.h mblas_mpfr.h 
sed -i.bak "s/MBLAS_GENERIC_H/MBLAS_MPFR_H/g" mblas_mpfr.h
sed -i.bak "s/INTEGER/mpackint/g" mblas_mpfr.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mblas_mpfr.h
sed -i.bak "s/REAL/mpreal/g" mblas_mpfr.h
sed -i.bak "s/COMPLEX/mpcomplex/g" mblas_mpfr.h
sed -i.bak "s/Mlsame/Mlsame_mpfr/g" mblas_mpfr.h
sed -i.bak "s/Mxerbla/Mxerbla_mpfr/g" mblas_mpfr.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mpack_config.h>\\
#include <mpreal.h>\\
#include <mpcomplex.h>\\
#include <mutils_mpfr.h>/" mblas_mpfr.h

cp mblas_generic.h mblas_gmp.h 
sed -i.bak "s/MBLAS_GENERIC_H/MBLAS_GMP_H/g" mblas_gmp.h
sed -i.bak "s/INTEGER/mpackint/g" mblas_gmp.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mblas_gmp.h
sed -i.bak "s/REAL/mpf_class/g" mblas_gmp.h
sed -i.bak "s/COMPLEX/mpc_class/g" mblas_gmp.h
sed -i.bak "s/Mlsame/Mlsame_gmp/g" mblas_gmp.h
sed -i.bak "s/Mxerbla/Mxerbla_gmp/g" mblas_gmp.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mpack_config.h>\\
#include <gmpxx.h>\\
#include <mpc_class.h>\\
#include <mutils_gmp.h>/" mblas_gmp.h

cp mblas_generic.h mblas_qd.h 
sed -i.bak "s/MBLAS_GENERIC_H/MBLAS_QD_H/g" mblas_qd.h
sed -i.bak "s/INTEGER/mpackint/g" mblas_qd.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mblas_qd.h
sed -i.bak "s/REAL/qd_real/g" mblas_qd.h
sed -i.bak "s/COMPLEX/qd_complex/g" mblas_qd.h
sed -i.bak "s/Mlsame/Mlsame_qd/g" mblas_qd.h
sed -i.bak "s/Mxerbla/Mxerbla_qd/g" mblas_qd.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mpack_config.h>\\
#include <qd\/qd_real.h>\\
#include <qd_complex.h>\\
#include <mutils_qd.h>/" mblas_qd.h

cp mblas_generic.h mblas_dd.h 
sed -i.bak "s/MBLAS_GENERIC_H/MBLAS_DD_H/g" mblas_dd.h
sed -i.bak "s/INTEGER/mpackint/g" mblas_dd.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mblas_dd.h
sed -i.bak "s/REAL/dd_real/g" mblas_dd.h
sed -i.bak "s/COMPLEX/dd_complex/g" mblas_dd.h
sed -i.bak "s/Mlsame/Mlsame_dd/g" mblas_dd.h
sed -i.bak "s/Mxerbla/Mxerbla_dd/g" mblas_dd.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mpack_config.h>\\
#include <qd\/dd_real.h>\\
#include <dd_complex.h>\\
#include <mutils_dd.h>/" mblas_dd.h

cp mblas_generic.h mblas_double.h 
sed -i.bak "s/MBLAS_GENERIC_H/MBLAS_DOUBLE_H/g" mblas_double.h
sed -i.bak "s/INTEGER/mpackint/g" mblas_double.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mblas_double.h
sed -i.bak "s/REAL/double/g" mblas_double.h
sed -i.bak "s/COMPLEX/std::complex<double>/g" mblas_double.h
sed -i.bak "s/Mlsame/Mlsame_double/g" mblas_double.h
sed -i.bak "s/Mxerbla/Mxerbla_double/g" mblas_double.h
sed -i.bak "s/%%REPLACEHERE%%/#include <mpack_config.h>\\
#include <mutils_double.h>/" mblas_double.h

