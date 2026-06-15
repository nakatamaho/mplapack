mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cheev_test_mpfr Rsyev_Frank_mpfr Rsyev_test_mpfr Rsyevd_DingDong_mpfr Rsyevd_Frank_mpfr Rsyevr_Frank_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cheev_test_mpfr_SOURCES = Cheev_test_mpfr.cpp
Cheev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cheev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyev_Frank_mpfr_SOURCES = Rsyev_Frank_mpfr.cpp
Rsyev_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyev_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyev_test_mpfr_SOURCES = Rsyev_test_mpfr.cpp
Rsyev_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyev_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyevd_DingDong_mpfr_SOURCES = Rsyevd_DingDong_mpfr.cpp
Rsyevd_DingDong_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyevd_DingDong_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyevd_Frank_mpfr_SOURCES = Rsyevd_Frank_mpfr.cpp
Rsyevd_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyevd_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsyevr_Frank_mpfr_SOURCES = Rsyevr_Frank_mpfr.cpp
Rsyevr_Frank_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsyevr_Frank_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cheev_test_gmp Rsyev_Frank_gmp Rsyev_test_gmp Rsyevd_DingDong_gmp Rsyevd_Frank_gmp Rsyevr_Frank_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cheev_test_gmp_SOURCES = Cheev_test_gmp.cpp
Cheev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cheev_test_gmp_LDFLAGS = $(gmp_libdepends)

Rsyev_Frank_gmp_SOURCES = Rsyev_Frank_gmp.cpp
Rsyev_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyev_Frank_gmp_LDFLAGS = $(gmp_libdepends)

Rsyev_test_gmp_SOURCES = Rsyev_test_gmp.cpp
Rsyev_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyev_test_gmp_LDFLAGS = $(gmp_libdepends)

Rsyevd_DingDong_gmp_SOURCES = Rsyevd_DingDong_gmp.cpp
Rsyevd_DingDong_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyevd_DingDong_gmp_LDFLAGS = $(gmp_libdepends)

Rsyevd_Frank_gmp_SOURCES = Rsyevd_Frank_gmp.cpp
Rsyevd_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyevd_Frank_gmp_LDFLAGS = $(gmp_libdepends)

Rsyevr_Frank_gmp_SOURCES = Rsyevr_Frank_gmp.cpp
Rsyevr_Frank_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsyevr_Frank_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cheev_test_binary128 Rsyev_Frank_binary128 Rsyev_test_binary128 Rsyevd_DingDong_binary128 Rsyevd_Frank_binary128 Rsyevr_Frank_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128
endif

Cheev_test_binary128_SOURCES = Cheev_test_binary128.cpp
Cheev_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Cheev_test_binary128_LDFLAGS = $(binary128_libdepends)

Rsyev_Frank_binary128_SOURCES = Rsyev_Frank_binary128.cpp
Rsyev_Frank_binary128_CXXFLAGS = $(binary128_cxxflags)
Rsyev_Frank_binary128_LDFLAGS = $(binary128_libdepends)

Rsyev_test_binary128_SOURCES = Rsyev_test_binary128.cpp
Rsyev_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Rsyev_test_binary128_LDFLAGS = $(binary128_libdepends)

Rsyevd_DingDong_binary128_SOURCES = Rsyevd_DingDong_binary128.cpp
Rsyevd_DingDong_binary128_CXXFLAGS = $(binary128_cxxflags)
Rsyevd_DingDong_binary128_LDFLAGS = $(binary128_libdepends)

Rsyevd_Frank_binary128_SOURCES = Rsyevd_Frank_binary128.cpp
Rsyevd_Frank_binary128_CXXFLAGS = $(binary128_cxxflags)
Rsyevd_Frank_binary128_LDFLAGS = $(binary128_libdepends)

Rsyevr_Frank_binary128_SOURCES = Rsyevr_Frank_binary128.cpp
Rsyevr_Frank_binary128_CXXFLAGS = $(binary128_cxxflags)
Rsyevr_Frank_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cheev_test_binary80 Rsyev_Frank_binary80 Rsyev_test_binary80 Rsyevd_DingDong_binary80 Rsyevd_Frank_binary80 Rsyevr_Frank_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -L$(top_builddir)/mpblas/reference -lmpblas_binary80

Cheev_test_binary80_SOURCES = Cheev_test_binary80.cpp
Cheev_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Cheev_test_binary80_LDFLAGS = $(binary80_libdepends)

Rsyev_Frank_binary80_SOURCES = Rsyev_Frank_binary80.cpp
Rsyev_Frank_binary80_CXXFLAGS = $(binary80_cxxflags)
Rsyev_Frank_binary80_LDFLAGS = $(binary80_libdepends)

Rsyev_test_binary80_SOURCES = Rsyev_test_binary80.cpp
Rsyev_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Rsyev_test_binary80_LDFLAGS = $(binary80_libdepends)

Rsyevd_DingDong_binary80_SOURCES = Rsyevd_DingDong_binary80.cpp
Rsyevd_DingDong_binary80_CXXFLAGS = $(binary80_cxxflags)
Rsyevd_DingDong_binary80_LDFLAGS = $(binary80_libdepends)

Rsyevd_Frank_binary80_SOURCES = Rsyevd_Frank_binary80.cpp
Rsyevd_Frank_binary80_CXXFLAGS = $(binary80_cxxflags)
Rsyevd_Frank_binary80_LDFLAGS = $(binary80_libdepends)

Rsyevr_Frank_binary80_SOURCES = Rsyevr_Frank_binary80.cpp
Rsyevr_Frank_binary80_CXXFLAGS = $(binary80_cxxflags)
Rsyevr_Frank_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cheev_test_double Rsyev_Frank_double Rsyev_test_double Rsyevd_DingDong_double Rsyevd_Frank_double Rsyevr_Frank_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Cheev_test_double_SOURCES = Cheev_test_double.cpp
Cheev_test_double_CXXFLAGS = $(double_cxxflags)
Cheev_test_double_LDFLAGS = $(double_libdepends)

Rsyev_Frank_double_SOURCES = Rsyev_Frank_double.cpp
Rsyev_Frank_double_CXXFLAGS = $(double_cxxflags)
Rsyev_Frank_double_LDFLAGS = $(double_libdepends)

Rsyev_test_double_SOURCES = Rsyev_test_double.cpp
Rsyev_test_double_CXXFLAGS = $(double_cxxflags)
Rsyev_test_double_LDFLAGS = $(double_libdepends)

Rsyevd_DingDong_double_SOURCES = Rsyevd_DingDong_double.cpp
Rsyevd_DingDong_double_CXXFLAGS = $(double_cxxflags)
Rsyevd_DingDong_double_LDFLAGS = $(double_libdepends)

Rsyevd_Frank_double_SOURCES = Rsyevd_Frank_double.cpp
Rsyevd_Frank_double_CXXFLAGS = $(double_cxxflags)
Rsyevd_Frank_double_LDFLAGS = $(double_libdepends)

Rsyevr_Frank_double_SOURCES = Rsyevr_Frank_double.cpp
Rsyevr_Frank_double_CXXFLAGS = $(double_cxxflags)
Rsyevr_Frank_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cheev_test_dd Rsyev_Frank_dd Rsyev_test_dd Rsyevd_DingDong_dd Rsyevd_Frank_dd Rsyevr_Frank_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cheev_test_dd_SOURCES = Cheev_test_dd.cpp
Cheev_test_dd_CXXFLAGS = $(dd_cxxflags)
Cheev_test_dd_LDFLAGS = $(dd_libdepends)

Rsyev_Frank_dd_SOURCES = Rsyev_Frank_dd.cpp
Rsyev_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rsyev_Frank_dd_LDFLAGS = $(dd_libdepends)

Rsyev_test_dd_SOURCES = Rsyev_test_dd.cpp
Rsyev_test_dd_CXXFLAGS = $(dd_cxxflags)
Rsyev_test_dd_LDFLAGS = $(dd_libdepends)

Rsyevd_DingDong_dd_SOURCES = Rsyevd_DingDong_dd.cpp
Rsyevd_DingDong_dd_CXXFLAGS = $(dd_cxxflags)
Rsyevd_DingDong_dd_LDFLAGS = $(dd_libdepends)

Rsyevd_Frank_dd_SOURCES = Rsyevd_Frank_dd.cpp
Rsyevd_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rsyevd_Frank_dd_LDFLAGS = $(dd_libdepends)

Rsyevr_Frank_dd_SOURCES = Rsyevr_Frank_dd.cpp
Rsyevr_Frank_dd_CXXFLAGS = $(dd_cxxflags)
Rsyevr_Frank_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cheev_test_qd Rsyev_Frank_qd Rsyev_test_qd Rsyevd_DingDong_qd Rsyevd_Frank_qd Rsyevr_Frank_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cheev_test_qd_SOURCES = Cheev_test_qd.cpp
Cheev_test_qd_CXXFLAGS = $(qd_cxxflags)
Cheev_test_qd_LDFLAGS = $(qd_libdepends)

Rsyev_Frank_qd_SOURCES = Rsyev_Frank_qd.cpp
Rsyev_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rsyev_Frank_qd_LDFLAGS = $(qd_libdepends)

Rsyev_test_qd_SOURCES = Rsyev_test_qd.cpp
Rsyev_test_qd_CXXFLAGS = $(qd_cxxflags)
Rsyev_test_qd_LDFLAGS = $(qd_libdepends)

Rsyevd_DingDong_qd_SOURCES = Rsyevd_DingDong_qd.cpp
Rsyevd_DingDong_qd_CXXFLAGS = $(qd_cxxflags)
Rsyevd_DingDong_qd_LDFLAGS = $(qd_libdepends)

Rsyevd_Frank_qd_SOURCES = Rsyevd_Frank_qd.cpp
Rsyevd_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rsyevd_Frank_qd_LDFLAGS = $(qd_libdepends)

Rsyevr_Frank_qd_SOURCES = Rsyevr_Frank_qd.cpp
Rsyevr_Frank_qd_CXXFLAGS = $(qd_cxxflags)
Rsyevr_Frank_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/03_SymmetricEigenproblems

mplapackexamples_DATA = Cheev_test_mpfr.cpp Cheev_test_gmp.cpp Cheev_test_binary128.cpp Cheev_test_binary80.cpp Cheev_test_double.cpp Cheev_test_dd.cpp Cheev_test_qd.cpp Rsyev_Frank_mpfr.cpp Rsyev_Frank_gmp.cpp Rsyev_Frank_binary128.cpp Rsyev_Frank_binary80.cpp Rsyev_Frank_double.cpp Rsyev_Frank_dd.cpp Rsyev_Frank_qd.cpp Rsyev_test_mpfr.cpp Rsyev_test_gmp.cpp Rsyev_test_binary128.cpp Rsyev_test_binary80.cpp Rsyev_test_double.cpp Rsyev_test_dd.cpp Rsyev_test_qd.cpp Rsyevd_DingDong_mpfr.cpp Rsyevd_DingDong_gmp.cpp Rsyevd_DingDong_binary128.cpp Rsyevd_DingDong_binary80.cpp Rsyevd_DingDong_double.cpp Rsyevd_DingDong_dd.cpp Rsyevd_DingDong_qd.cpp Rsyevd_Frank_mpfr.cpp Rsyevd_Frank_gmp.cpp Rsyevd_Frank_binary128.cpp Rsyevd_Frank_binary80.cpp Rsyevd_Frank_double.cpp Rsyevd_Frank_dd.cpp Rsyevd_Frank_qd.cpp Rsyevr_Frank_mpfr.cpp Rsyevr_Frank_gmp.cpp Rsyevr_Frank_binary128.cpp Rsyevr_Frank_binary80.cpp Rsyevr_Frank_double.cpp Rsyevr_Frank_dd.cpp Rsyevr_Frank_qd.cpp \
 Makefile.freebsd Makefile.linux Makefile.macos Makefile.mingw

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif

EXTRA_DIST = \
	Makefile.freebsd.in \
	Makefile.linux.in \
	Makefile.macos.in \
	Makefile.mingw.in

GENERATED_MAKEFILES = \
	Makefile.freebsd \
	Makefile.linux \
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

Makefile.macos: Makefile.macos.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

Makefile.mingw: Makefile.mingw.in
	sed $(MPLAPACK_SED_SUBST) $< > $@
	bash $(top_srcdir)/misc/filter-example-programs.sh $@ $(DISABLED_EXAMPLE_SUFFIXES)

all-local: $(GENERATED_MAKEFILES)

clean-local:
	rm -f $(GENERATED_MAKEFILES)
