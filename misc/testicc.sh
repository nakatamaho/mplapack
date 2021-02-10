#!/bin/sh
#test on Linux
MPACK_VER=0.7.0
HOST=sdpa.indsys.chuo-u.ac.jp

MPACK_DIST=mpack-$MPACK_VER.tar.gz

WORKDIR=mpack-icc-work
MPACK_CONFIG_SH=configicc.sh
MPACK_DIST_DIR=mpack-$MPACK_VER

WORKDIR_FULL=mpack-icc-work-full
MPACK_CONFIG_SH_FULL=configicc_full.sh
MPACK_DIST_DIR_FULL=mpack-$MPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH ; make -j20) 2>&1 | tee log.build.icc" | tee log.build.icc
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.icc" | tee log.check.icc

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH_FULL ; make -j20) 2>&1 | tee log.build.full.icc" | tee log.build.full.icc
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.full.icc" | tee log.check.full.icc

