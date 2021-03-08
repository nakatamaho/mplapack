cd build
cd gcc/mplapack-0.8.0
/usr/bin/time make -k check >& ../../../log.checkall.gcc.sh
cd ../..
cd icc/mplapack-0.8.0
/usr/bin/time make -k check >& ../../../log.checkall.icc.sh
cd ../..
