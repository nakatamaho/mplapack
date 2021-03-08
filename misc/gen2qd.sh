sed -i.bak "s/INTEGER/mplapackint/g" *.cpp
sed -i.bak "s/REAL/qd_real/g" *.cpp
sed -i.bak "s/COMPLEX/qd_complex/g" *.cpp
sed -i.bak "s/mpblas.h/mpblas_qd.h/g" *.cpp
