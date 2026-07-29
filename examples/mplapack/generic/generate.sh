# usage
# cd /home/docker/mplapack/examples/mplapack/00_LinearEquations/generic ; bash -x ../../generate.sh
# or
# cd /home/docker/mplapack/examples/mplapack/03_SymmetricEigenproblems/generic ; bash -x ../../generate.sh
# etc..
FILES=`ls *_generic.cpp`
pushd .. ; _MATFILES=`ls M*.txt` ; popd
MATFILES=`echo $_MATFILES`
MPLIBS="mpfr gmp binary128 binary80 double dd qd"

if [ `uname` = "Darwin" ]; then
    SED=gsed
else
    SED=sed
fi

_FILE=`echo $FILES | awk '{print $1}' | $SED 's/_generic\.cpp//g'`
$SED -e "s|%%ROUTINE%%|$_FILE|g" ../../generic/Makefile.freebsd.in > ../Makefile.freebsd.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" ../../generic/Makefile.macos.in   > ../Makefile.macos.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" ../../generic/Makefile.linux.in   > ../Makefile.linux.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" ../../generic/Makefile.linux.inteloneAPI.in > ../Makefile.linux.inteloneAPI.in
$SED -e "s|%%ROUTINE%%|$_FILE|g" ../../generic/Makefile.mingw.in   > ../Makefile.mingw.in

SOURCEFILES=""
for _file in $FILES; do
    MPLIBS="mpfr gmp binary128 binary80 double dd qd"
    if [ "$_file" = "Cgeev_NPR_generic.cpp" ]; then
        MPLIBS="mpfr binary128 binary80 double dd qd"
    fi
    for _mplib in $MPLIBS; do
        resultfilename=`echo $_file | $SED "s/generic/${_mplib}/g"`
        if echo $_file | grep ^C ; then
            cat ../../generic/header_${_mplib}_complex ${_file} > ../$resultfilename
        else
            cat ../../generic/header_${_mplib} ${_file} > ../$resultfilename
        fi
        SOURCEFILES=`echo $SOURCEFILES ${resultfilename}`
        if [ x"$_mplib" = x"gmp" ]; then
            $SED -i -e "s/REAL/mpf_class/g" -e "s/COMPLEX/mpc_class/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"mpfr" ]; then
            $SED -i -e "s/REAL/mpreal/g" -e "s/COMPLEX/mpcomplex/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"double" ]; then
            $SED -i -e "s/REAL/double/g" -e "s/COMPLEX/std::complex<double>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi 
        if [ x"$_mplib" = x"dd" ]; then
            $SED -i -e "s/REAL/dd_real/g" -e "s/COMPLEX/dd_complex/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi 
        if [ x"$_mplib" = x"qd" ]; then
            $SED -i -e "s/REAL/qd_real/g" -e "s/COMPLEX/qd_complex/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"binary128" ]; then
            $SED -i -e "s/REAL/mplapack_binary128_t/g" -e "s/COMPLEX/std::complex<mplapack_binary128_t>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi
        if [ x"$_mplib" = x"binary80" ]; then
            $SED -i -e "s/REAL/mplapack_binary80_t/g" -e "s/COMPLEX/std::complex<mplapack_binary80_t>/g" -e "s/INTEGER/mplapackint/g" -e "s/InTEGER/INTEGER_${_mplib}/g" -e "s/ReAL/REAL_${_mplib}/g" -e "s/Rlamch/Rlamch_${_mplib}/g" -e "s/%%MPLIB%%/${_mplib}/g" ../$resultfilename
        fi
    done
done

echo "mplapackexamples_PROGRAMS =" > ../Makefile.am

MPLIBS="mpfr gmp binary128 binary80 double dd qd"
for _mplib in $MPLIBS; do
    if [ x"$_mplib" = x"mpfr" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_MPFR" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(top_srcdir)/mpfrc++ -I\$(GMP_INCLUDEDIR) -I\$(MPFR_INCLUDEDIR) -I\$(MPC_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_${_mplib} -L\$(MPC_LIBDIR) -L\$(MPFR_LIBDIR) -L\$(GMP_LIBDIR) -lmpfr -lmpc -lgmp"  >> ../Makefile.am
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
        executefilenames=`echo $FILES | $SED 's/ Cgeev_NPR_generic\.cpp//g' | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(GMP_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_${_mplib} -L\$(GMP_LIBDIR) -lgmp"  >> ../Makefile.am
        echo ""               >> ../Makefile.am
        for _file in $FILES; do
            if [ "$_file" = "Cgeev_NPR_generic.cpp" ]; then
                continue
            fi
            A=`echo $_file | $SED "s/generic\.cpp/${_mplib}/g"`
            echo "${A}_SOURCES = ${A}.cpp" >> ../Makefile.am
            echo "${A}_CXXFLAGS = \$(${_mplib}_cxxflags)" >> ../Makefile.am
            echo "${A}_LDFLAGS = \$(${_mplib}_libdepends)" >> ../Makefile.am
            echo ""               >> ../Makefile.am
        done
        echo "endif"             >> ../Makefile.am
    fi

    if [ x"$_mplib" = x"binary128" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_BINARY128" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "if MPLAPACK_BINARY128_MODE_QUADMATH" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_${_mplib} -lquadmath"  >> ../Makefile.am
        echo "else" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_${_mplib}"  >> ../Makefile.am
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

    if [ x"$_mplib" = x"binary80" ]; then
        echo ""               >> ../Makefile.am
        echo "if ENABLE_BINARY80" >> ../Makefile.am
        executefilenames=`echo $FILES | $SED 's/\.cpp//g' | $SED "s/generic/${_mplib}/g"`
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_${_mplib}"  >> ../Makefile.am
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
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_${_mplib}"  >> ../Makefile.am
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
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(QD_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_dd -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am
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
        echo "mplapackexamples_PROGRAMS += $executefilenames" >> ../Makefile.am
        echo ""               >> ../Makefile.am
        echo "${_mplib}_cxxflags = \$(OPENMP_CXXFLAGS) -I\$(top_srcdir)/include -I\$(QD_INCLUDEDIR)" >> ../Makefile.am
        echo "${_mplib}_libdepends = -Wl,-rpath,\$(libdir) -L\$(top_builddir)/mplapack/reference -lmplapack_qd -L\$(QD_LIBDIR) -lqd"  >> ../Makefile.am
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
array=( `echo $path | tr -s '/' ' '`)
kind_index=`expr ${#array[@]} - 2`
echo "mplapackexamplesdir=\$(prefix)/share/examples/mplapack/${array[${kind_index}]}"   >> ../Makefile.am
echo ""               >> ../Makefile.am
echo "mplapackexamples_DATA = $SOURCEFILES \\" >> ../Makefile.am
echo "$MATFILES Makefile.freebsd Makefile.linux Makefile.linux.inteloneAPI Makefile.macos Makefile.mingw" >> ../Makefile.am
echo ""               >> ../Makefile.am
cat >> ../Makefile.am << 'EOF'
install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif

EXTRA_DIST = \
	Makefile.freebsd.in \
	Makefile.linux.in \
	Makefile.linux.inteloneAPI.in \
	Makefile.macos.in \
	Makefile.mingw.in

GENERATED_MAKEFILES = \
	Makefile.freebsd \
	Makefile.linux \
	Makefile.linux.inteloneAPI \
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
        -e 's|%%LIBQUADMATH%%| -lquadmath|g' \
        -e 's|%%QUADMATH_CPPFLAGS%%|$(QUADMATH_CPPFLAGS)|g' \
        -e 's|%%QUADMATH_LDFLAGS%%|$(QUADMATH_LDFLAGS)|g'
else
MPLAPACK_SED_SUBST = \
        -e 's|%%MPLAPACKDIR%%|$(prefix)|g' \
        -e 's|%%LIBQUADMATH%%||g' \
        -e 's|%%QUADMATH_CPPFLAGS%%||g' \
        -e 's|%%QUADMATH_LDFLAGS%%||g'
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

Makefile.linux.inteloneAPI: Makefile.linux.inteloneAPI.in
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
