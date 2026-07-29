mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cposv_example_mpfr Cpotrf_potrs_example_mpfr Rpocon_example_mpfr Rposv_example_mpfr Rpotrf_potrs_example_mpfr Rptsv_example_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cposv_example_mpfr_SOURCES = Cposv_example_mpfr.cpp
Cposv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cposv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Cpotrf_potrs_example_mpfr_SOURCES = Cpotrf_potrs_example_mpfr.cpp
Cpotrf_potrs_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cpotrf_potrs_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rpocon_example_mpfr_SOURCES = Rpocon_example_mpfr.cpp
Rpocon_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rpocon_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rposv_example_mpfr_SOURCES = Rposv_example_mpfr.cpp
Rposv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rposv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rpotrf_potrs_example_mpfr_SOURCES = Rpotrf_potrs_example_mpfr.cpp
Rpotrf_potrs_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rpotrf_potrs_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rptsv_example_mpfr_SOURCES = Rptsv_example_mpfr.cpp
Rptsv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rptsv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cposv_example_gmp Cpotrf_potrs_example_gmp Rpocon_example_gmp Rposv_example_gmp Rpotrf_potrs_example_gmp Rptsv_example_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Cposv_example_gmp_SOURCES = Cposv_example_gmp.cpp
Cposv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cposv_example_gmp_LDFLAGS = $(gmp_libdepends)

Cpotrf_potrs_example_gmp_SOURCES = Cpotrf_potrs_example_gmp.cpp
Cpotrf_potrs_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cpotrf_potrs_example_gmp_LDFLAGS = $(gmp_libdepends)

Rpocon_example_gmp_SOURCES = Rpocon_example_gmp.cpp
Rpocon_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rpocon_example_gmp_LDFLAGS = $(gmp_libdepends)

Rposv_example_gmp_SOURCES = Rposv_example_gmp.cpp
Rposv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rposv_example_gmp_LDFLAGS = $(gmp_libdepends)

Rpotrf_potrs_example_gmp_SOURCES = Rpotrf_potrs_example_gmp.cpp
Rpotrf_potrs_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rpotrf_potrs_example_gmp_LDFLAGS = $(gmp_libdepends)

Rptsv_example_gmp_SOURCES = Rptsv_example_gmp.cpp
Rptsv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rptsv_example_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cposv_example_binary128 Cpotrf_potrs_example_binary128 Rpocon_example_binary128 Rposv_example_binary128 Rpotrf_potrs_example_binary128 Rptsv_example_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Cposv_example_binary128_SOURCES = Cposv_example_binary128.cpp
Cposv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cposv_example_binary128_LDFLAGS = $(binary128_libdepends)

Cpotrf_potrs_example_binary128_SOURCES = Cpotrf_potrs_example_binary128.cpp
Cpotrf_potrs_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cpotrf_potrs_example_binary128_LDFLAGS = $(binary128_libdepends)

Rpocon_example_binary128_SOURCES = Rpocon_example_binary128.cpp
Rpocon_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rpocon_example_binary128_LDFLAGS = $(binary128_libdepends)

Rposv_example_binary128_SOURCES = Rposv_example_binary128.cpp
Rposv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rposv_example_binary128_LDFLAGS = $(binary128_libdepends)

Rpotrf_potrs_example_binary128_SOURCES = Rpotrf_potrs_example_binary128.cpp
Rpotrf_potrs_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rpotrf_potrs_example_binary128_LDFLAGS = $(binary128_libdepends)

Rptsv_example_binary128_SOURCES = Rptsv_example_binary128.cpp
Rptsv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rptsv_example_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cposv_example_binary80 Cpotrf_potrs_example_binary80 Rpocon_example_binary80 Rposv_example_binary80 Rpotrf_potrs_example_binary80 Rptsv_example_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Cposv_example_binary80_SOURCES = Cposv_example_binary80.cpp
Cposv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cposv_example_binary80_LDFLAGS = $(binary80_libdepends)

Cpotrf_potrs_example_binary80_SOURCES = Cpotrf_potrs_example_binary80.cpp
Cpotrf_potrs_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cpotrf_potrs_example_binary80_LDFLAGS = $(binary80_libdepends)

Rpocon_example_binary80_SOURCES = Rpocon_example_binary80.cpp
Rpocon_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rpocon_example_binary80_LDFLAGS = $(binary80_libdepends)

Rposv_example_binary80_SOURCES = Rposv_example_binary80.cpp
Rposv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rposv_example_binary80_LDFLAGS = $(binary80_libdepends)

Rpotrf_potrs_example_binary80_SOURCES = Rpotrf_potrs_example_binary80.cpp
Rpotrf_potrs_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rpotrf_potrs_example_binary80_LDFLAGS = $(binary80_libdepends)

Rptsv_example_binary80_SOURCES = Rptsv_example_binary80.cpp
Rptsv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rptsv_example_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cposv_example_double Cpotrf_potrs_example_double Rpocon_example_double Rposv_example_double Rpotrf_potrs_example_double Rptsv_example_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Cposv_example_double_SOURCES = Cposv_example_double.cpp
Cposv_example_double_CXXFLAGS = $(double_cxxflags)
Cposv_example_double_LDFLAGS = $(double_libdepends)

Cpotrf_potrs_example_double_SOURCES = Cpotrf_potrs_example_double.cpp
Cpotrf_potrs_example_double_CXXFLAGS = $(double_cxxflags)
Cpotrf_potrs_example_double_LDFLAGS = $(double_libdepends)

Rpocon_example_double_SOURCES = Rpocon_example_double.cpp
Rpocon_example_double_CXXFLAGS = $(double_cxxflags)
Rpocon_example_double_LDFLAGS = $(double_libdepends)

Rposv_example_double_SOURCES = Rposv_example_double.cpp
Rposv_example_double_CXXFLAGS = $(double_cxxflags)
Rposv_example_double_LDFLAGS = $(double_libdepends)

Rpotrf_potrs_example_double_SOURCES = Rpotrf_potrs_example_double.cpp
Rpotrf_potrs_example_double_CXXFLAGS = $(double_cxxflags)
Rpotrf_potrs_example_double_LDFLAGS = $(double_libdepends)

Rptsv_example_double_SOURCES = Rptsv_example_double.cpp
Rptsv_example_double_CXXFLAGS = $(double_cxxflags)
Rptsv_example_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cposv_example_dd Cpotrf_potrs_example_dd Rpocon_example_dd Rposv_example_dd Rpotrf_potrs_example_dd Rptsv_example_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Cposv_example_dd_SOURCES = Cposv_example_dd.cpp
Cposv_example_dd_CXXFLAGS = $(dd_cxxflags)
Cposv_example_dd_LDFLAGS = $(dd_libdepends)

Cpotrf_potrs_example_dd_SOURCES = Cpotrf_potrs_example_dd.cpp
Cpotrf_potrs_example_dd_CXXFLAGS = $(dd_cxxflags)
Cpotrf_potrs_example_dd_LDFLAGS = $(dd_libdepends)

Rpocon_example_dd_SOURCES = Rpocon_example_dd.cpp
Rpocon_example_dd_CXXFLAGS = $(dd_cxxflags)
Rpocon_example_dd_LDFLAGS = $(dd_libdepends)

Rposv_example_dd_SOURCES = Rposv_example_dd.cpp
Rposv_example_dd_CXXFLAGS = $(dd_cxxflags)
Rposv_example_dd_LDFLAGS = $(dd_libdepends)

Rpotrf_potrs_example_dd_SOURCES = Rpotrf_potrs_example_dd.cpp
Rpotrf_potrs_example_dd_CXXFLAGS = $(dd_cxxflags)
Rpotrf_potrs_example_dd_LDFLAGS = $(dd_libdepends)

Rptsv_example_dd_SOURCES = Rptsv_example_dd.cpp
Rptsv_example_dd_CXXFLAGS = $(dd_cxxflags)
Rptsv_example_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cposv_example_qd Cpotrf_potrs_example_qd Rpocon_example_qd Rposv_example_qd Rpotrf_potrs_example_qd Rptsv_example_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Cposv_example_qd_SOURCES = Cposv_example_qd.cpp
Cposv_example_qd_CXXFLAGS = $(qd_cxxflags)
Cposv_example_qd_LDFLAGS = $(qd_libdepends)

Cpotrf_potrs_example_qd_SOURCES = Cpotrf_potrs_example_qd.cpp
Cpotrf_potrs_example_qd_CXXFLAGS = $(qd_cxxflags)
Cpotrf_potrs_example_qd_LDFLAGS = $(qd_libdepends)

Rpocon_example_qd_SOURCES = Rpocon_example_qd.cpp
Rpocon_example_qd_CXXFLAGS = $(qd_cxxflags)
Rpocon_example_qd_LDFLAGS = $(qd_libdepends)

Rposv_example_qd_SOURCES = Rposv_example_qd.cpp
Rposv_example_qd_CXXFLAGS = $(qd_cxxflags)
Rposv_example_qd_LDFLAGS = $(qd_libdepends)

Rpotrf_potrs_example_qd_SOURCES = Rpotrf_potrs_example_qd.cpp
Rpotrf_potrs_example_qd_CXXFLAGS = $(qd_cxxflags)
Rpotrf_potrs_example_qd_LDFLAGS = $(qd_libdepends)

Rptsv_example_qd_SOURCES = Rptsv_example_qd.cpp
Rptsv_example_qd_CXXFLAGS = $(qd_cxxflags)
Rptsv_example_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/01_PositiveDefiniteLinearEquations

mplapackexamples_DATA = Cposv_example_mpfr.cpp Cposv_example_gmp.cpp Cposv_example_binary128.cpp Cposv_example_binary80.cpp Cposv_example_double.cpp Cposv_example_dd.cpp Cposv_example_qd.cpp Cpotrf_potrs_example_mpfr.cpp Cpotrf_potrs_example_gmp.cpp Cpotrf_potrs_example_binary128.cpp Cpotrf_potrs_example_binary80.cpp Cpotrf_potrs_example_double.cpp Cpotrf_potrs_example_dd.cpp Cpotrf_potrs_example_qd.cpp Rpocon_example_mpfr.cpp Rpocon_example_gmp.cpp Rpocon_example_binary128.cpp Rpocon_example_binary80.cpp Rpocon_example_double.cpp Rpocon_example_dd.cpp Rpocon_example_qd.cpp Rposv_example_mpfr.cpp Rposv_example_gmp.cpp Rposv_example_binary128.cpp Rposv_example_binary80.cpp Rposv_example_double.cpp Rposv_example_dd.cpp Rposv_example_qd.cpp Rpotrf_potrs_example_mpfr.cpp Rpotrf_potrs_example_gmp.cpp Rpotrf_potrs_example_binary128.cpp Rpotrf_potrs_example_binary80.cpp Rpotrf_potrs_example_double.cpp Rpotrf_potrs_example_dd.cpp Rpotrf_potrs_example_qd.cpp Rptsv_example_mpfr.cpp Rptsv_example_gmp.cpp Rptsv_example_binary128.cpp Rptsv_example_binary80.cpp Rptsv_example_double.cpp Rptsv_example_dd.cpp Rptsv_example_qd.cpp \
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
