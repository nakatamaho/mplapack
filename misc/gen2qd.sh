sed -i.bak "s/INTEGER/mpackint/g" *.cpp
sed -i.bak "s/REAL/qd_real/g" *.cpp
sed -i.bak "s/COMPLEX/qd_complex/g" *.cpp
sed -i.bak "s/mblas.h/mblas_qd.h/g" *.cpp
