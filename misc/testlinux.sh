#!/bin/sh
#test on Linux
MPLAPACK_VER=0.8.0
HOST=172.27.42.56

MPLAPACK_DIST=mplapack-$MPLAPACK_VER.tar.gz

WORKDIR=mplapack-work
MPLAPACK_CONFIG_SH=configlinux.sh
MPLAPACK_DIST_DIR=mplapack-$MPLAPACK_VER

WORKDIR_FULL=mplapack-work-full
MPLAPACK_CONFIG_SH_FULL=configlinux_full.sh
MPLAPACK_DIST_DIR_FULL=mplapack-$MPLAPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPLAPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH ; make -j20) 2>&1 | tee log.build.linux" | tee log.build.linux
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.linux" | tee log.check.linux

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPLAPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH_FULL ; make -j20) 2>&1 | tee log.build.full.linux" | tee log.build.full.linux
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.full.linux" | tee log.check.full.linux

