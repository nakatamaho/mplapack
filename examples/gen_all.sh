TOP=`pwd`

cd $TOP
cd $TOP/mpblas
find . -maxdepth 1 ! -name '.' \
                   ! -name 'Makefile' \
                   ! -name 'Makefile.am' \
                   ! -name 'Makefile.in' \
                   ! -name 'generic' \
                   -exec rm -rf {} +
cd generic
cd mpblas/generic ; bash -x generate.sh

cd $TOP/mplapack
DIRS=`ls -d */ | sed 's|generic/||g'`
echo $DIRS
for _dir in $DIRS; do
    cd $TOP/mplapack/$_dir
    # Remove everything except Makefile, Makefile.am, Makefile.in, and generic/
    find . -maxdepth 1 ! -name '.' \
                       ! -name 'Makefile' \
                       ! -name 'Makefile.am' \
                       ! -name 'Makefile.in' \
                       ! -name 'Matrix*txt' \
                       ! -name 'generic' \
                       -exec rm -rf {} +
    cd generic
    bash -x ../../generic/generate.sh
done
