#!/bin/sh
#test on Linux
#ssh -N -L2003:tesla:22 maho@sdpa.indsys.chuo-u.ac.jp &
#add following lines to /home/maho/.ssh/config
#Host fedora13
#Hostname localhost
#User maho
#Port 2003

MPACK_VER=0.7.0
HOST=fedora13
WORKDIR=mpack-work-fedora13
MPACK_CONFIG_SH=configlinux.sh

MPACK_DIST_DIR=mpack-$MPACK_VER
MPACK_DIST=mpack-$MPACK_VER.tar.gz

#install mpack
ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPACK_DIST $HOST:$WORKDIR
scp misc/$MPACK_CONFIG_SH $HOST:$WORKDIR
ssh $HOST "cd $WORKDIR ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST ; cd $MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH ; make -j4 ; make install) 2>&1 | tee log.linux" | tee log.linux_fedora13

