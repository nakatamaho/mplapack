rm -f FILELIST

for _file in `cat DONE | awk '{print $3}'`; do 
echo "${_file} \\" >> FILELIST
done

sed -e "/%%insert here%%/e cat FILELIST" Makefile.am.in > Makefile.am
sed -i -e "s/%%insert here%%//g" Makefile.am
head -c -1 Makefile.am > Makefile.am_
mv Makefile.am_ Makefile.am
sed -i -e '$s/\\//' Makefile.am
