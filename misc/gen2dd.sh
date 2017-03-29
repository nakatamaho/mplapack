sed -i.bak "s/INTEGER/mpackint/g" *.cpp
sed -i.bak "s/REAL/dd_real/g" *.cpp
sed -i.bak "s/COMPLEX/dd_complex/g" *.cpp
sed -i.bak "s/mblas.h/mblas_dd.h/g" *.cpp
