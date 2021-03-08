sed -i.bak "s/INTEGER/mplapackint/g" *.cpp
sed -i.bak "s/REAL/dd_real/g" *.cpp
sed -i.bak "s/COMPLEX/dd_complex/g" *.cpp
sed -i.bak "s/mpblas.h/mpblas_dd.h/g" *.cpp
