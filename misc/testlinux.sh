#!/bin/sh
#test on Linux
MPACK_VER=0.8.0
HOST=172.27.42.56

MPACK_DIST=mpack-$MPACK_VER.tar.gz

WORKDIR=mpack-work
MPACK_CONFIG_SH=configlinux.sh
MPACK_DIST_DIR=mpack-$MPACK_VER

WORKDIR_FULL=mpack-work-full
MPACK_CONFIG_SH_FULL=configlinux_full.sh
MPACK_DIST_DIR_FULL=mpack-$MPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH ; make -j20) 2>&1 | tee log.build.linux" | tee log.build.linux
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.linux" | tee log.check.linux

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH_FULL ; make -j20) 2>&1 | tee log.build.full.linux" | tee log.build.full.linux
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.full.linux" | tee log.check.full.linux

