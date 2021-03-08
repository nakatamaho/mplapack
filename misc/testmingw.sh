#!/bin/sh
#test mingw build on FreeBSD
MPLAPACK_VER=0.7.0
HOST=172.27.42.154

MPLAPACK_DIST=mplapack-$MPLAPACK_VER.tar.gz

WORKDIR=mplapack-mingw-work
MPLAPACK_CONFIG_SH=configmingw.sh
MPLAPACK_DIST_DIR=mplapack-$MPLAPACK_VER

WORKDIR_FULL=mplapack-mingw-work-full
MPLAPACK_CONFIG_SH_FULL=configmingw_full.sh
MPLAPACK_DIST_DIR_FULL=mplapack-$MPLAPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPLAPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH ; gmake -j8) 2>&1 | tee log.build.mingw" | tee log.build.mingw
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPLAPACK_DIST_DIR ;  TESTS_ENVIRONMENT=wine gmake -k check 2>&1 | tee log.check.mingw" | tee log.check.mingw

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPLAPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPLAPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPLAPACK_DIST ; tar xvfz $MPLAPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ; (sh ../$MPLAPACK_CONFIG_SH_FULL ; gmake -j8) 2>&1 | tee log.build.full.mingw" | tee log.build.full.mingw
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPLAPACK_DIST_DIR ; TESTS_ENVIRONMENT=wine gmake -k check 2>&1 | tee log.check.full.mingw" | tee log.check.full.mingw

