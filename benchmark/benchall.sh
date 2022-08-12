FILES=`ls go*.sh`

for _file in $FILES; do
    bash -x $_file
done

FILES=`ls *.plt`

for _file in $FILES; do
    gnuplot $_file > ${_file%.*}.pdf
done
