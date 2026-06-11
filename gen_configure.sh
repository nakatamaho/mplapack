#!/bin/bash
aclocal
cd mplapack/test/compare/
bash ./gen.Makefile.am.sh     
cd ../../..
autoheader
automake -a -v --add-missing
autoconf
