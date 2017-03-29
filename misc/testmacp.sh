#!/bin/sh
#test on MacOSX
MPACK_VER=0.6.7
HOST=172.27.42.30
WORKDIR=mpack-work
MPACK_CONFIG_SH=configmacp.sh

MPACK_DIST_DIR=mpack-$MPACK_VER
MPACK_DIST=mpack-$MPACK_VER.tar.gz

#install mpack
ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPACK_DIST $HOST:$WORKDIR
scp misc/$MPACK_CONFIG_SH $HOST:$WORKDIR
ssh $HOST "cd $WORKDIR ; md5 $MPACK_DIST ; tar xvfz $MPACK_DIST ; cd $MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH ; make ; make install) | tee log.macp" | tee log.macp
