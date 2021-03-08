#!/bin/sh
#test build on FBSD
MPLAPACK_VER=0.7.0
HOST=172.27.42.231

MPLAPACK_DIST=mplapack-$MPLAPACK_VER.tar.gz

WORKDIR=mplapack-work
MPLAPACK_CONFIG_SH=configfbsd.sh
MPLAPACK_DIST_DIR=mplapack-$MPLAPACK_VER

WORKDIR_FULL=mplapack-work-full
MPLAPACK_CONFIG_SH_FULL=configfbsd_full.sh
MPLAPACK_DIST_DIR_FULL=mplapack-$MPLAPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPLAPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH ; gmake -j8) 2>&1 | tee log.build.fbsd" | tee log.build.fbsd
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ;  gmake -k check 2>&1 | tee log.check.fbsd" | tee log.check.fbsd

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPLAPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH_FULL ; gmake -j8) 2>&1 | tee log.build.full.fbsd" | tee log.build.full.fbsd
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ;  gmake -k check 2>&1 | tee log.check.full.fbsd" | tee log.check.full.fbsd

