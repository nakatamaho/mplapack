# usage
# cd /home/docker/mplapack/examples/mpblas/generic ; bash -x generate.sh
FILES=`ls R*generic.cpp C*generic.cpp`
MPLIBS="mpfr gmp _Float128 _Float64x double dd qd"
if [ `uname` = "Darwin" ]; then
    SED=gsed
else
    SED=sed
fi

_FILE=`ls R*.cpp | head -1 | $SED 's/_/ /g' | awk '{print $1}'`
$SED -e "s|%%ROUTINE%%|$_FILE|g" Makefile.freebsd.in > ../Makefile.freebsd.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" Makefile.macos.in   > ../Makefile.macos.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" Makefile.linux.in   > ../Makefile.linux.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" Makefile.mingw.in   > ../Makefile.mingw.in

SOURCEFILES=""
rm -f source_files
for _file in $FILES; do
    for _mplib in $MPLIBS; do
        resultfilename=`echo $_file | $SED "s/generic/${_mplib}/g"`
        if echo $_file | grep ^R ; then
            cat header_${_mplib} ${_file} > ../$resultfilename
        elif echo $_file | grep ^C ; then
            cat header_${_mplib}_complex ${_file} > ../$resultfilename
        else
            echo "unknown type"
	    exit -1
        fi 
        SOURCEFILES=`echo $SOURCEFILES ${resultfilename}`
        if [ x"$_mplib" = x"gmp" ]; then
            $SED -i -e "s/REAL/mpf_class/g" -e "s/COMPLEX/mpc_class/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"mpfr" ]; then
            $SED -i -e "s/REAL/mpreal/g" -e "s/COMPLEX/mpcomplex/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"double" ]; then
            $SED -i -e "s/REAL/double/g" -e "s/COMPLEX/std::complex<double>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi 
        if [ x"$_mplib" = x"dd" ]; then
            $SED -i -e "s/REAL/dd_real/g" -e "s/COMPLEX/dd_complex/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi 
        if [ x"$_mplib" = x"qd" ]; then
            $SED -i -e "s/REAL/qd_real/g" -e "s/COMPLEX/qd_complex/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"_Float128" ]; then
            $SED -i -e "s/REAL/_Float128/g" -e "s/COMPLEX/std::complex<_Float128>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"_Float64x" ]; then
            $SED -i -e "s/REAL/_Float64x/g" -e "s/COMPLEX/std::complex<_Float64x>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
    done
done

echo "mpblasexamples_PROGRAMS =" > ../Makefile.am

for _mplib in $MPLIBS; do

    if [ x"$_mplib" = x"mpfr" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_MPFR" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(top_srcdir)/mpfrc++ -I\$(GMP_INCLUDEDIR) -I\$(MPFR_INCLUDEDIR) -I\$(MPC_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(MPC_LIBDIR) -L\$(MPFR_LIBDIR) -L\$(GMP_LIBDIR) -lmpfr -lmpc -lgmp"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"gmp" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_GMP" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(GMP_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(GMP_LIBDIR) -lgmp"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"_Float128" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE__FLOAT128" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "if WANT_QUADMATH" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -lquadmath"  >> ../Makefile.am
        echo "else" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib}"  >> ../Makefile.am   
        echo "endif" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"_Float64x" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE__FLOAT64X" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib}"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"double" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_DOUBLE" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib}"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"dd" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_DD" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(QD_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am   
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"qd" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_QD" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(QD_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am   
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi
done
echo ""               >> ../Makefile.am

path=`pwd` 
echo "mpblasexamplesdir=\$(prefix)/share/examples/mpblas/"   >> ../Makefile.am
echo ""               >> ../Makefile.am
echo "mpblasexamples_DATA = $SOURCEFILES \\" >> ../Makefile.am
echo "Makefile.freebsd Makefile.linux Makefile.macos Makefile.mingw" >> ../Makefile.am
echo ""               >> ../Makefile.am
cat >> ../Makefile.am << EOF
install-data-hook:
if IS_MACOS
	bash \$(top_builddir)/misc/fix_dylib_macOS.sh \$(mpblasexamplesdir) \$(prefix)
endif
EOF
