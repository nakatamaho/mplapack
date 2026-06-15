# usage
# cd /home/docker/mplapack/examples/mpblas/generic ; bash -x generate.sh
FILES=`ls R*generic.cpp C*generic.cpp`
MPLIBS="mpfr gmp binary128 binary80 double dd qd"
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
cp Makefile.linux_cuda.in ../Makefile.linux_cuda.in

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
        if [ x"$_mplib" = x"binary128" ]; then
            $SED -i -e "s/REAL/mplapack_binary128_t/g" -e "s/COMPLEX/std::complex<mplapack_binary128_t>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"binary80" ]; then
            $SED -i -e "s/REAL/mplapack_binary80_t/g" -e "s/COMPLEX/std::complex<mplapack_binary80_t>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" ../$resultfilename
        fi
    done
done

echo "mpblasexamples_PROGRAMS =" > ../Makefile.am

for _mplib in $MPLIBS; do

    if [ x"$_mplib" = x"mpfr" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_MPFR" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/${_mplib}/${_mplib}_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(top_srcdir)/mpfrc++ -I\$(GMP_INCLUDEDIR) -I\$(MPFR_INCLUDEDIR) -I\$(MPC_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(MPC_LIBDIR) -L\$(MPFR_LIBDIR) -L\$(GMP_LIBDIR) -lmpfr -lmpc -lgmp"  >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt -L\$(MPC_LIBDIR) -L\$(MPFR_LIBDIR) -L\$(GMP_LIBDIR) -lmpfr -lmpc -lgmp"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"gmp" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_GMP" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/${_mplib}/${_mplib}_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(GMP_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(GMP_LIBDIR) -lgmp"  >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt -L\$(GMP_LIBDIR) -lgmp"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"binary128" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_BINARY128" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/binary128/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/binary128/binary128_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib}" >> ../Makefile.am
        echo "if MPLAPACK_BINARY128_MODE_QUADMATH" >> ../Makefile.am
        echo "${_mplib}_libdepends += -lquadmath" >> ../Makefile.am
        echo "endif" >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt" >> ../Makefile.am
        echo "if MPLAPACK_BINARY128_MODE_QUADMATH" >> ../Makefile.am
        echo "${_mplib}_opt_libdepends += -lquadmath" >> ../Makefile.am
        echo "endif" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"binary80" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_BINARY80" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/binary80/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/binary80/binary80_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib}"  >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"double" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_DOUBLE" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/${_mplib}/${_mplib}_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib}"  >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"dd" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_DD" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/${_mplib}/${_mplib}_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(QD_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"qd" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_QD" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        executefilenames_opt=`echo $executefilenames | $SED "s/${_mplib}/${_mplib}_opt/g"`
        echo "mpblasexamples_PROGRAMS += $executefilenames $executefilenames_opt" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(QD_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -L\$(top_builddir)/mpblas/reference -lmpblas_${_mplib} -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am
        echo "${_mplib}_opt_libdepends = -L\$(top_builddir)/mpblas/optimized/${_mplib} -lmpblas_${_mplib}_opt -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        for _file in $FILES; do
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_opt_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_opt_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_opt_LDFLAGS = \$(${_mplib}_opt_libdepends)" >> ../Makefile.am
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
echo "Makefile.freebsd Makefile.linux Makefile.linux_cuda Makefile.macos Makefile.mingw" >> ../Makefile.am
echo ""               >> ../Makefile.am
cat >> ../Makefile.am << EOF
install-data-hook:
if IS_MACOS
	bash \$(top_srcdir)/misc/fix_dylib_macOS.sh \$(mpblasexamplesdir) \$(prefix)
endif
EOF

cat >> ../Makefile.am << 'EOF'

EXTRA_DIST = \
	Makefile.freebsd.in \
	Makefile.linux.in \
	Makefile.linux_cuda.in \
	Makefile.macos.in \
	Makefile.mingw.in

GENERATED_MAKEFILES = \
	Makefile.freebsd \
	Makefile.linux \
	Makefile.linux_cuda \
	Makefile.macos \
	Makefile.mingw

# Resolve %%LIBQUADMATH%% based on whether quadmath backend is selected.
if MPLAPACK_BINARY128_MODE_QUADMATH
LIBQUADMATH_SUBST = -lquadmath
else
LIBQUADMATH_SUBST =
endif

if MPLAPACK_BINARY128_MODE_QUADMATH
MPLAPACK_SED_SUBST = \
        -e 's|%%MPLAPACKDIR%%|$(prefix)|g' \
        -e 's|%%LIBQUADMATH%%| -lquadmath|g'
else
MPLAPACK_SED_SUBST = \
        -e 's|%%MPLAPACKDIR%%|$(prefix)|g' \
        -e 's|%%LIBQUADMATH%%||g'
endif

DISABLED_EXAMPLE_SUFFIXES =
if ENABLE_MPFR
else
DISABLED_EXAMPLE_SUFFIXES += mpfr
endif
if ENABLE_GMP
else
DISABLED_EXAMPLE_SUFFIXES += gmp
endif
if ENABLE_QD
else
DISABLED_EXAMPLE_SUFFIXES += qd
endif
if ENABLE_DD
else
DISABLED_EXAMPLE_SUFFIXES += dd
endif
if ENABLE_DOUBLE
else
DISABLED_EXAMPLE_SUFFIXES += double
endif
if ENABLE_BINARY128
else
DISABLED_EXAMPLE_SUFFIXES += binary128
endif
if ENABLE_BINARY80
else
DISABLED_EXAMPLE_SUFFIXES += binary80
endif

Makefile.freebsd: Makefile.freebsd.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.linux: Makefile.linux.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.linux_cuda: Makefile.linux_cuda.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.macos: Makefile.macos.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.mingw: Makefile.mingw.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

all-local: $(GENERATED_MAKEFILES)

clean-local:
	rm -f $(GENERATED_MAKEFILES)
EOF
