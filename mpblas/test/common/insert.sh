FILES=`ls *test.cpp | grep -v iCamax`

for _file in $FILES; do
    echo $_file 
    LINE=`cat -n $_file | grep fail | grep ing | awk '{print $1}'`
    if [ x$LINE != x"" ]; then
      sed -i "$((LINE+2))d" $_file 
      sed -i "$((LINE+1))r insert.txt" $_file 
    fi
done