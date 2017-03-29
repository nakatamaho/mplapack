#!/bin/sh
#test on Linux
MPACK_VER=0.7.0
HOST=sdpa.indsys.chuo-u.ac.jp
WORKDIR=mpack-work-icc
MPACK_CONFIG_SH=configlinux_sdpa_icc.sh

MPACK_DIST_DIR=mpack-$MPACK_VER
MPACK_DIST=mpack-$MPACK_VER.tar.gz

#install mpack
ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPACK_DIST $HOST:$WORKDIR
scp misc/$MPACK_CONFIG_SH $HOST:$WORKDIR
ssh $HOST "cd $WORKDIR ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST ; cd $MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH ; make -j4 ; make check ; make install) 2>&1 | tee log.linux" | tee log.linux_sdpa_icc
