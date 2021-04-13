bash -x ~/mplapack/misc/conv_all_lapack.sh 
grep -e ldR -e ldC *cpp | awk '{print $NF}' | sort | uniq | sed 's/];//g' | sed -e 's/ldC/z/g' -e 's/ldR/d/g' | awk '{print "and prev_tok.value != \"" $1 "\"" }' 
cp hand/*.cpp  hand/*.hpp .
#patch -p0 < patch