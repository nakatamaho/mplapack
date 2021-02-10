sed -i.bak "s/INTEGER/mpackint/g" *.cpp
sed -i.bak "s/REAL/mpf_class/g" *.cpp
sed -i.bak "s/COMPLEX/mpc_class/g" *.cpp
sed -i.bak "s/mblas.h/mblas_gmp.h/g" *.cpp
