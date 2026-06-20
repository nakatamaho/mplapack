mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgetri_test_mpfr Rgesv_test_mpfr Rgetri_Hilbert_mpfr Rgetri_test2_mpfr Rgetri_test3_mpfr Rgetri_test_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgetri_test_mpfr_SOURCES = Cgetri_test_mpfr.cpp
Cgetri_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgetri_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgesv_test_mpfr_SOURCES = Rgesv_test_mpfr.cpp
Rgesv_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgesv_test_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_Hilbert_mpfr_SOURCES = Rgetri_Hilbert_mpfr.cpp
Rgetri_Hilbert_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_Hilbert_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_test2_mpfr_SOURCES = Rgetri_test2_mpfr.cpp
Rgetri_test2_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_test2_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_test3_mpfr_SOURCES = Rgetri_test3_mpfr.cpp
Rgetri_test3_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_test3_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgetri_test_mpfr_SOURCES = Rgetri_test_mpfr.cpp
Rgetri_test_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgetri_test_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cgetri_test_gmp Rgesv_test_gmp Rgetri_Hilbert_gmp Rgetri_test2_gmp Rgetri_test3_gmp Rgetri_test_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

Cgetri_test_gmp_SOURCES = Cgetri_test_gmp.cpp
Cgetri_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgetri_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgesv_test_gmp_SOURCES = Rgesv_test_gmp.cpp
Rgesv_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgesv_test_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_Hilbert_gmp_SOURCES = Rgetri_Hilbert_gmp.cpp
Rgetri_Hilbert_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_Hilbert_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_test2_gmp_SOURCES = Rgetri_test2_gmp.cpp
Rgetri_test2_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_test2_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_test3_gmp_SOURCES = Rgetri_test3_gmp.cpp
Rgetri_test3_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_test3_gmp_LDFLAGS = $(gmp_libdepends)

Rgetri_test_gmp_SOURCES = Rgetri_test_gmp.cpp
Rgetri_test_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgetri_test_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cgetri_test_binary128 Rgesv_test_binary128 Rgetri_Hilbert_binary128 Rgetri_test2_binary128 Rgetri_test3_binary128 Rgetri_test_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128
endif

Cgetri_test_binary128_SOURCES = Cgetri_test_binary128.cpp
Cgetri_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgetri_test_binary128_LDFLAGS = $(binary128_libdepends)

Rgesv_test_binary128_SOURCES = Rgesv_test_binary128.cpp
Rgesv_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgesv_test_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_Hilbert_binary128_SOURCES = Rgetri_Hilbert_binary128.cpp
Rgetri_Hilbert_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_Hilbert_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_test2_binary128_SOURCES = Rgetri_test2_binary128.cpp
Rgetri_test2_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_test2_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_test3_binary128_SOURCES = Rgetri_test3_binary128.cpp
Rgetri_test3_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_test3_binary128_LDFLAGS = $(binary128_libdepends)

Rgetri_test_binary128_SOURCES = Rgetri_test_binary128.cpp
Rgetri_test_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgetri_test_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cgetri_test_binary80 Rgesv_test_binary80 Rgetri_Hilbert_binary80 Rgetri_test2_binary80 Rgetri_test3_binary80 Rgetri_test_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -L$(top_builddir)/mpblas/reference -lmpblas_binary80

Cgetri_test_binary80_SOURCES = Cgetri_test_binary80.cpp
Cgetri_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgetri_test_binary80_LDFLAGS = $(binary80_libdepends)

Rgesv_test_binary80_SOURCES = Rgesv_test_binary80.cpp
Rgesv_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgesv_test_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_Hilbert_binary80_SOURCES = Rgetri_Hilbert_binary80.cpp
Rgetri_Hilbert_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_Hilbert_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_test2_binary80_SOURCES = Rgetri_test2_binary80.cpp
Rgetri_test2_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_test2_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_test3_binary80_SOURCES = Rgetri_test3_binary80.cpp
Rgetri_test3_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_test3_binary80_LDFLAGS = $(binary80_libdepends)

Rgetri_test_binary80_SOURCES = Rgetri_test_binary80.cpp
Rgetri_test_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgetri_test_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgetri_test_double Rgesv_test_double Rgetri_Hilbert_double Rgetri_test2_double Rgetri_test3_double Rgetri_test_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

Cgetri_test_double_SOURCES = Cgetri_test_double.cpp
Cgetri_test_double_CXXFLAGS = $(double_cxxflags)
Cgetri_test_double_LDFLAGS = $(double_libdepends)

Rgesv_test_double_SOURCES = Rgesv_test_double.cpp
Rgesv_test_double_CXXFLAGS = $(double_cxxflags)
Rgesv_test_double_LDFLAGS = $(double_libdepends)

Rgetri_Hilbert_double_SOURCES = Rgetri_Hilbert_double.cpp
Rgetri_Hilbert_double_CXXFLAGS = $(double_cxxflags)
Rgetri_Hilbert_double_LDFLAGS = $(double_libdepends)

Rgetri_test2_double_SOURCES = Rgetri_test2_double.cpp
Rgetri_test2_double_CXXFLAGS = $(double_cxxflags)
Rgetri_test2_double_LDFLAGS = $(double_libdepends)

Rgetri_test3_double_SOURCES = Rgetri_test3_double.cpp
Rgetri_test3_double_CXXFLAGS = $(double_cxxflags)
Rgetri_test3_double_LDFLAGS = $(double_libdepends)

Rgetri_test_double_SOURCES = Rgetri_test_double.cpp
Rgetri_test_double_CXXFLAGS = $(double_cxxflags)
Rgetri_test_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cgetri_test_dd Rgesv_test_dd Rgetri_Hilbert_dd Rgetri_test2_dd Rgetri_test3_dd Rgetri_test_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

Cgetri_test_dd_SOURCES = Cgetri_test_dd.cpp
Cgetri_test_dd_CXXFLAGS = $(dd_cxxflags)
Cgetri_test_dd_LDFLAGS = $(dd_libdepends)

Rgesv_test_dd_SOURCES = Rgesv_test_dd.cpp
Rgesv_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgesv_test_dd_LDFLAGS = $(dd_libdepends)

Rgetri_Hilbert_dd_SOURCES = Rgetri_Hilbert_dd.cpp
Rgetri_Hilbert_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_Hilbert_dd_LDFLAGS = $(dd_libdepends)

Rgetri_test2_dd_SOURCES = Rgetri_test2_dd.cpp
Rgetri_test2_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_test2_dd_LDFLAGS = $(dd_libdepends)

Rgetri_test3_dd_SOURCES = Rgetri_test3_dd.cpp
Rgetri_test3_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_test3_dd_LDFLAGS = $(dd_libdepends)

Rgetri_test_dd_SOURCES = Rgetri_test_dd.cpp
Rgetri_test_dd_CXXFLAGS = $(dd_cxxflags)
Rgetri_test_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cgetri_test_qd Rgesv_test_qd Rgetri_Hilbert_qd Rgetri_test2_qd Rgetri_test3_qd Rgetri_test_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

Cgetri_test_qd_SOURCES = Cgetri_test_qd.cpp
Cgetri_test_qd_CXXFLAGS = $(qd_cxxflags)
Cgetri_test_qd_LDFLAGS = $(qd_libdepends)

Rgesv_test_qd_SOURCES = Rgesv_test_qd.cpp
Rgesv_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgesv_test_qd_LDFLAGS = $(qd_libdepends)

Rgetri_Hilbert_qd_SOURCES = Rgetri_Hilbert_qd.cpp
Rgetri_Hilbert_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_Hilbert_qd_LDFLAGS = $(qd_libdepends)

Rgetri_test2_qd_SOURCES = Rgetri_test2_qd.cpp
Rgetri_test2_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_test2_qd_LDFLAGS = $(qd_libdepends)

Rgetri_test3_qd_SOURCES = Rgetri_test3_qd.cpp
Rgetri_test3_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_test3_qd_LDFLAGS = $(qd_libdepends)

Rgetri_test_qd_SOURCES = Rgetri_test_qd.cpp
Rgetri_test_qd_CXXFLAGS = $(qd_cxxflags)
Rgetri_test_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/00_LinearEquations

mplapackexamples_DATA = Cgetri_test_mpfr.cpp Cgetri_test_gmp.cpp Cgetri_test_binary128.cpp Cgetri_test_binary80.cpp Cgetri_test_double.cpp Cgetri_test_dd.cpp Cgetri_test_qd.cpp Rgesv_test_mpfr.cpp Rgesv_test_gmp.cpp Rgesv_test_binary128.cpp Rgesv_test_binary80.cpp Rgesv_test_double.cpp Rgesv_test_dd.cpp Rgesv_test_qd.cpp Rgetri_Hilbert_mpfr.cpp Rgetri_Hilbert_gmp.cpp Rgetri_Hilbert_binary128.cpp Rgetri_Hilbert_binary80.cpp Rgetri_Hilbert_double.cpp Rgetri_Hilbert_dd.cpp Rgetri_Hilbert_qd.cpp Rgetri_test2_mpfr.cpp Rgetri_test2_gmp.cpp Rgetri_test2_binary128.cpp Rgetri_test2_binary80.cpp Rgetri_test2_double.cpp Rgetri_test2_dd.cpp Rgetri_test2_qd.cpp Rgetri_test3_mpfr.cpp Rgetri_test3_gmp.cpp Rgetri_test3_binary128.cpp Rgetri_test3_binary80.cpp Rgetri_test3_double.cpp Rgetri_test3_dd.cpp Rgetri_test3_qd.cpp Rgetri_test_mpfr.cpp Rgetri_test_gmp.cpp Rgetri_test_binary128.cpp Rgetri_test_binary80.cpp Rgetri_test_double.cpp Rgetri_test_dd.cpp Rgetri_test_qd.cpp \
 Makefile.freebsd Makefile.linux Makefile.linux.inteloneAPI Makefile.macos Makefile.mingw

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
