
FILES=`ls *mpfr`
for _file in $FILES; do
libtool --mode=execute valgrind --leak-check=full --show-leak-kinds=all --main-stacksize=18388608 $_file 2>&1| tee ${_file}.log
done