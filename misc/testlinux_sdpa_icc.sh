#!/bin/sh
#test on Linux
MPLAPACK_VER=0.7.0
HOST=sdpa.indsys.chuo-u.ac.jp
WORKDIR=mplapack-work-icc
MPLAPACK_CONFIG_SH=configlinux_sdpa_icc.sh

MPLAPACK_DIST_DIR=mplapack-$MPLAPACK_VER
MPLAPACK_DIST=mplapack-$MPLAPACK_VER.tar.gz

#install mplapack
ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPLAPACK_DIST $HOST:$WORKDIR
scp misc/$MPLAPACK_CONFIG_SH $HOST:$WORKDIR
ssh $HOST "cd $WORKDIR ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST ; cd $MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH ; make -j4 ; make check ; make install) 2>&1 | tee log.linux" | tee log.linux_sdpa_icc
