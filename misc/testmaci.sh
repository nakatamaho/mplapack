#!/bin/sh
#test on MacOSX

MPLAPACK_VER=0.7.0
HOST=172.27.42.137

MPLAPACK_DIST=mplapack-$MPLAPACK_VER.tar.gz

WORKDIR=mplapack-work
MPLAPACK_CONFIG_SH=configmaci.sh
MPLAPACK_DIST_DIR=mplapack-$MPLAPACK_VER

WORKDIR_FULL=mplapack-work-full
MPLAPACK_CONFIG_SH_FULL=configmaci_full.sh
MPLAPACK_DIST_DIR_FULL=mplapack-$MPLAPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPLAPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH ; make -j8) 2>&1 | tee log.build.maci" | tee log.build.maci
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ;  make -k check 2>&1 | tee log.check.maci" | tee log.check.maci

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPLAPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH_FULL ; make -j8) 2>&1 | tee log.build.full.maci" | tee log.build.full.maci
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ; make -k check 2>&1 | tee log.check.full.maci" | tee log.check.full.maci

