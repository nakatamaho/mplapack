#!/bin/sh
#test mingw build on FreeBSD
MPACK_VER=0.7.0
HOST=172.27.42.154

MPACK_DIST=mpack-$MPACK_VER.tar.gz

WORKDIR=mpack-mingw-work
MPACK_CONFIG_SH=configmingw.sh
MPACK_DIST_DIR=mpack-$MPACK_VER

WORKDIR_FULL=mpack-mingw-work-full
MPACK_CONFIG_SH_FULL=configmingw_full.sh
MPACK_DIST_DIR_FULL=mpack-$MPACK_VER

ssh $HOST "rm -rf $WORKDIR"
ssh $HOST "mkdir $WORKDIR"
scp $MPACK_DIST $HOST:$WORKDIR
/usr/bin/time scp misc/$MPACK_CONFIG_SH $HOST:$WORKDIR
/usr/bin/time ssh $HOST "cd $WORKDIR ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH ; gmake -j8) 2>&1 | tee log.build.mingw" | tee log.build.mingw
/usr/bin/time ssh $HOST "cd $WORKDIR/$MPACK_DIST_DIR ;  TESTS_ENVIRONMENT=wine gmake -k check 2>&1 | tee log.check.mingw" | tee log.check.mingw

ssh $HOST "rm -rf $WORKDIR_FULL"
ssh $HOST "mkdir $WORKDIR_FULL"
scp $MPACK_DIST $HOST:$WORKDIR_FULL
/usr/bin/time scp misc/$MPACK_CONFIG_SH_FULL $HOST:$WORKDIR_FULL
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL ; md5sum $MPACK_DIST ; tar xvfz $MPACK_DIST"
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPACK_DIST_DIR ; (sh ../$MPACK_CONFIG_SH_FULL ; gmake -j8) 2>&1 | tee log.build.full.mingw" | tee log.build.full.mingw
/usr/bin/time ssh $HOST "cd $WORKDIR_FULL/$MPACK_DIST_DIR ; TESTS_ENVIRONMENT=wine gmake -k check 2>&1 | tee log.check.full.mingw" | tee log.check.full.mingw

