cd ~/
rm -rf mplapack
git clone --depth 1 --branch release/2.1 git@github.com:nakatamaho/mplapack.git
cd mplapack
(git --no-pager log -1 ; LANG=C date ; bash misc/reconfig.macOS.sh ; make -j4 ; make install ; LANG=C date ; make check -j 4; LANG=C date)
