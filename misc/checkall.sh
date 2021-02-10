cd build
cd gcc/mpack-0.8.0
/usr/bin/time make -k check >& ../../../log.checkall.gcc.sh
cd ../..
cd icc/mpack-0.8.0
/usr/bin/time make -k check >& ../../../log.checkall.icc.sh
cd ../..
