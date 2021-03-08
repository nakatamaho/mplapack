#!/bin/sh
#test on MacOSX
MPLAPACK_VER=0.6.7
HOST=172.27.42.30
WORKDIR=mplapack-work
MPLAPACK_CONFIG_SH=configmacp.sh

MPLAPACK_DIST_DIR=mplapack-$MPLAPACK_VER
MPLAPACK_DIST=mplapack-$MPLAPACK_VER.tar.gz

#install mplapack
ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPLAPACK_DIST $HOST:$WORKDIR
scp misc/$MPLAPACK_CONFIG_SH $HOST:$WORKDIR
ssh $HOST "cd $WORKDIR ; md5 $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST ; cd $MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH ; make ; make install) | tee log.macp" | tee log.macp
