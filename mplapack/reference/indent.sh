FILES=`ls *cpp`

for _file in $FILES; do
    bash ~/mplapack/misc/indent.sh $_file
done