cd ~/
rm -rf mplapack
git clone --depth 1 --branch release/2.1 git@github.com:nakatamaho/mplapack.git
cd mplapack
(git log -1 ; LANG=C date ; bash misc/reconfig.ubuntu24.04.mingw64.sh ; make -j32 ; make install ; LANG=C date ; make LOG_COMPILER=wine check -j 32; LANG=C date)
