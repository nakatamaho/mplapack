TOP=`pwd`

cd $TOP
cd mpblas/generic ; bash -x generate.sh

cd $TOP
cd mplapack
DIRS=`ls -d */ | sed 's|generic/||g'`
echo $DIRS

for _dir in $DIRS; do
    cd $TOP
    cd mplapack/$_dir/generic ; bash -x ../../generic/generate.sh
done
