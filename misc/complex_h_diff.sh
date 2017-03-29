    sed 's/qd_complex/COMPLEX/g' qd_complex.h | sed 's/qd_real/REAL/g' > ll
    sed 's/dd_complex/COMPLEX/g' dd_complex.h | sed 's/dd_real/REAL/g' > ll2
    sed 's/mpc_class/COMPLEX/g'  mpc_class.h  | sed 's/mpf_class/REAL/g' > ll3
    sed 's/mpcomplex/COMPLEX/g'  ../mpfrc++/mpcomplex.h  | sed 's/mpreal/REAL/g' > ll4
    sed 's/double_complex/COMPLEX/g'  double_complex.h  | sed 's/double/REAL/g' > ll5
