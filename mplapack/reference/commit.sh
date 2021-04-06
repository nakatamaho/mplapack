if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi
bash -x ~/mplapack/misc/indent.sh $1
git add $1
git commit $1 -m "."