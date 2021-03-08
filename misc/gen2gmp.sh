sed -i.bak "s/INTEGER/mplapackint/g" *.cpp
sed -i.bak "s/REAL/mpf_class/g" *.cpp
sed -i.bak "s/COMPLEX/mpc_class/g" *.cpp
sed -i.bak "s/mpblas.h/mpblas_gmp.h/g" *.cpp
