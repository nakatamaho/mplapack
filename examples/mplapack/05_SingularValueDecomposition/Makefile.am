mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgesvd_test_mpfr Rgesvd_random_highcond_mpfr Rgesvd_readfromfile_mpfr Rgesvd_test2_mpfr Rgesvd_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgesvd_test_mpfr_SOURCES = Cgesvd_test_mpfr.cpp
Cgesvd_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgesvd_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_random_highcond_mpfr_SOURCES = Rgesvd_random_highcond_mpfr.cpp
Rgesvd_random_highcond_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_random_highcond_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_readfromfile_mpfr_SOURCES = Rgesvd_readfromfile_mpfr.cpp
Rgesvd_readfromfile_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_readfromfile_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_test2_mpfr_SOURCES = Rgesvd_test2_mpfr.cpp
Rgesvd_test2_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_test2_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesvd_test_mpfr_SOURCES = Rgesvd_test_mpfr.cpp
Rgesvd_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesvd_test_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cgesvd_test_gmp Rgesvd_random_highcond_gmp Rgesvd_readfromfile_gmp Rgesvd_test2_gmp Rgesvd_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Cgesvd_test_gmp_SOURCES = Cgesvd_test_gmp.cpp
Cgesvd_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgesvd_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_random_highcond_gmp_SOURCES = Rgesvd_random_highcond_gmp.cpp
Rgesvd_random_highcond_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_random_highcond_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_readfromfile_gmp_SOURCES = Rgesvd_readfromfile_gmp.cpp
Rgesvd_readfromfile_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_readfromfile_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_test2_gmp_SOURCES = Rgesvd_test2_gmp.cpp
Rgesvd_test2_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_test2_gmp_LDFLAGS = $(gmp_libdepends)

Rgesvd_test_gmp_SOURCES = Rgesvd_test_gmp.cpp
Rgesvd_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesvd_test_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cgesvd_test_binary128 Rgesvd_random_highcond_binary128 Rgesvd_readfromfile_binary128 Rgesvd_test2_binary128 Rgesvd_test_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Cgesvd_test_binary128_SOURCES = Cgesvd_test_binary128.cpp
Cgesvd_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgesvd_test_binary128_LDFLAGS = $(binary128_libdepends)

Rgesvd_random_highcond_binary128_SOURCES = Rgesvd_random_highcond_binary128.cpp
Rgesvd_random_highcond_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesvd_random_highcond_binary128_LDFLAGS = $(binary128_libdepends)

Rgesvd_readfromfile_binary128_SOURCES = Rgesvd_readfromfile_binary128.cpp
Rgesvd_readfromfile_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesvd_readfromfile_binary128_LDFLAGS = $(binary128_libdepends)

Rgesvd_test2_binary128_SOURCES = Rgesvd_test2_binary128.cpp
Rgesvd_test2_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesvd_test2_binary128_LDFLAGS = $(binary128_libdepends)

Rgesvd_test_binary128_SOURCES = Rgesvd_test_binary128.cpp
Rgesvd_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesvd_test_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cgesvd_test_binary80 Rgesvd_random_highcond_binary80 Rgesvd_readfromfile_binary80 Rgesvd_test2_binary80 Rgesvd_test_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Cgesvd_test_binary80_SOURCES = Cgesvd_test_binary80.cpp
Cgesvd_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgesvd_test_binary80_LDFLAGS = $(binary80_libdepends)

Rgesvd_random_highcond_binary80_SOURCES = Rgesvd_random_highcond_binary80.cpp
Rgesvd_random_highcond_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesvd_random_highcond_binary80_LDFLAGS = $(binary80_libdepends)

Rgesvd_readfromfile_binary80_SOURCES = Rgesvd_readfromfile_binary80.cpp
Rgesvd_readfromfile_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesvd_readfromfile_binary80_LDFLAGS = $(binary80_libdepends)

Rgesvd_test2_binary80_SOURCES = Rgesvd_test2_binary80.cpp
Rgesvd_test2_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesvd_test2_binary80_LDFLAGS = $(binary80_libdepends)

Rgesvd_test_binary80_SOURCES = Rgesvd_test_binary80.cpp
Rgesvd_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesvd_test_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgesvd_test_double Rgesvd_random_highcond_double Rgesvd_readfromfile_double Rgesvd_test2_double Rgesvd_test_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Cgesvd_test_double_SOURCES = Cgesvd_test_double.cpp
Cgesvd_test_double_CXXFLAGS = $(double_cxxflags)
Cgesvd_test_double_LDFLAGS = $(double_libdepends)

Rgesvd_random_highcond_double_SOURCES = Rgesvd_random_highcond_double.cpp
Rgesvd_random_highcond_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_random_highcond_double_LDFLAGS = $(double_libdepends)

Rgesvd_readfromfile_double_SOURCES = Rgesvd_readfromfile_double.cpp
Rgesvd_readfromfile_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_readfromfile_double_LDFLAGS = $(double_libdepends)

Rgesvd_test2_double_SOURCES = Rgesvd_test2_double.cpp
Rgesvd_test2_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_test2_double_LDFLAGS = $(double_libdepends)

Rgesvd_test_double_SOURCES = Rgesvd_test_double.cpp
Rgesvd_test_double_CXXFLAGS = $(double_cxxflags)
Rgesvd_test_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cgesvd_test_dd Rgesvd_random_highcond_dd Rgesvd_readfromfile_dd Rgesvd_test2_dd Rgesvd_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Cgesvd_test_dd_SOURCES = Cgesvd_test_dd.cpp
Cgesvd_test_dd_CXXFLAGS = $(dd_cxxflags)
Cgesvd_test_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_random_highcond_dd_SOURCES = Rgesvd_random_highcond_dd.cpp
Rgesvd_random_highcond_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_random_highcond_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_readfromfile_dd_SOURCES = Rgesvd_readfromfile_dd.cpp
Rgesvd_readfromfile_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_readfromfile_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_test2_dd_SOURCES = Rgesvd_test2_dd.cpp
Rgesvd_test2_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_test2_dd_LDFLAGS = $(dd_libdepends)

Rgesvd_test_dd_SOURCES = Rgesvd_test_dd.cpp
Rgesvd_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgesvd_test_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cgesvd_test_qd Rgesvd_random_highcond_qd Rgesvd_readfromfile_qd Rgesvd_test2_qd Rgesvd_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Cgesvd_test_qd_SOURCES = Cgesvd_test_qd.cpp
Cgesvd_test_qd_CXXFLAGS = $(qd_cxxflags)
Cgesvd_test_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_random_highcond_qd_SOURCES = Rgesvd_random_highcond_qd.cpp
Rgesvd_random_highcond_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_random_highcond_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_readfromfile_qd_SOURCES = Rgesvd_readfromfile_qd.cpp
Rgesvd_readfromfile_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_readfromfile_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_test2_qd_SOURCES = Rgesvd_test2_qd.cpp
Rgesvd_test2_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_test2_qd_LDFLAGS = $(qd_libdepends)

Rgesvd_test_qd_SOURCES = Rgesvd_test_qd.cpp
Rgesvd_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgesvd_test_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/05_SingularValueDecomposition

mplapackexamples_DATA = Cgesvd_test_mpfr.cpp Cgesvd_test_gmp.cpp Cgesvd_test_binary128.cpp Cgesvd_test_binary80.cpp Cgesvd_test_double.cpp Cgesvd_test_dd.cpp Cgesvd_test_qd.cpp Rgesvd_random_highcond_mpfr.cpp Rgesvd_random_highcond_gmp.cpp Rgesvd_random_highcond_binary128.cpp Rgesvd_random_highcond_binary80.cpp Rgesvd_random_highcond_double.cpp Rgesvd_random_highcond_dd.cpp Rgesvd_random_highcond_qd.cpp Rgesvd_readfromfile_mpfr.cpp Rgesvd_readfromfile_gmp.cpp Rgesvd_readfromfile_binary128.cpp Rgesvd_readfromfile_binary80.cpp Rgesvd_readfromfile_double.cpp Rgesvd_readfromfile_dd.cpp Rgesvd_readfromfile_qd.cpp Rgesvd_test2_mpfr.cpp Rgesvd_test2_gmp.cpp Rgesvd_test2_binary128.cpp Rgesvd_test2_binary80.cpp Rgesvd_test2_double.cpp Rgesvd_test2_dd.cpp Rgesvd_test2_qd.cpp Rgesvd_test_mpfr.cpp Rgesvd_test_gmp.cpp Rgesvd_test_binary128.cpp Rgesvd_test_binary80.cpp Rgesvd_test_double.cpp Rgesvd_test_dd.cpp Rgesvd_test_qd.cpp \
Matrix_01.txt Matrix_02.txt Matrix_03.txt Makefile.freebsd Makefile.linux Makefile.linux.inteloneAPI Makefile.macos Makefile.mingw

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
