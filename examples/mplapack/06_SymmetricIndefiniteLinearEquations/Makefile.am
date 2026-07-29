mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Chesv_example_mpfr Csysv_example_mpfr Rsysv_example_mpfr

mpfr_cxxflags = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Chesv_example_mpfr_SOURCES = Chesv_example_mpfr.cpp
Chesv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Chesv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Csysv_example_mpfr_SOURCES = Csysv_example_mpfr.cpp
Csysv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Csysv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rsysv_example_mpfr_SOURCES = Rsysv_example_mpfr.cpp
Rsysv_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rsysv_example_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Chesv_example_gmp Csysv_example_gmp Rsysv_example_gmp

gmp_cxxflags = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Chesv_example_gmp_SOURCES = Chesv_example_gmp.cpp
Chesv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Chesv_example_gmp_LDFLAGS = $(gmp_libdepends)

Csysv_example_gmp_SOURCES = Csysv_example_gmp.cpp
Csysv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Csysv_example_gmp_LDFLAGS = $(gmp_libdepends)

Rsysv_example_gmp_SOURCES = Rsysv_example_gmp.cpp
Rsysv_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rsysv_example_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Chesv_example_binary128 Csysv_example_binary128 Rsysv_example_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Chesv_example_binary128_SOURCES = Chesv_example_binary128.cpp
Chesv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Chesv_example_binary128_LDFLAGS = $(binary128_libdepends)

Csysv_example_binary128_SOURCES = Csysv_example_binary128.cpp
Csysv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Csysv_example_binary128_LDFLAGS = $(binary128_libdepends)

Rsysv_example_binary128_SOURCES = Rsysv_example_binary128.cpp
Rsysv_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rsysv_example_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Chesv_example_binary80 Csysv_example_binary80 Rsysv_example_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Chesv_example_binary80_SOURCES = Chesv_example_binary80.cpp
Chesv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Chesv_example_binary80_LDFLAGS = $(binary80_libdepends)

Csysv_example_binary80_SOURCES = Csysv_example_binary80.cpp
Csysv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Csysv_example_binary80_LDFLAGS = $(binary80_libdepends)

Rsysv_example_binary80_SOURCES = Rsysv_example_binary80.cpp
Rsysv_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rsysv_example_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Chesv_example_double Csysv_example_double Rsysv_example_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Chesv_example_double_SOURCES = Chesv_example_double.cpp
Chesv_example_double_CXXFLAGS = $(double_cxxflags)
Chesv_example_double_LDFLAGS = $(double_libdepends)

Csysv_example_double_SOURCES = Csysv_example_double.cpp
Csysv_example_double_CXXFLAGS = $(double_cxxflags)
Csysv_example_double_LDFLAGS = $(double_libdepends)

Rsysv_example_double_SOURCES = Rsysv_example_double.cpp
Rsysv_example_double_CXXFLAGS = $(double_cxxflags)
Rsysv_example_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Chesv_example_dd Csysv_example_dd Rsysv_example_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Chesv_example_dd_SOURCES = Chesv_example_dd.cpp
Chesv_example_dd_CXXFLAGS = $(dd_cxxflags)
Chesv_example_dd_LDFLAGS = $(dd_libdepends)

Csysv_example_dd_SOURCES = Csysv_example_dd.cpp
Csysv_example_dd_CXXFLAGS = $(dd_cxxflags)
Csysv_example_dd_LDFLAGS = $(dd_libdepends)

Rsysv_example_dd_SOURCES = Rsysv_example_dd.cpp
Rsysv_example_dd_CXXFLAGS = $(dd_cxxflags)
Rsysv_example_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Chesv_example_qd Csysv_example_qd Rsysv_example_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Chesv_example_qd_SOURCES = Chesv_example_qd.cpp
Chesv_example_qd_CXXFLAGS = $(qd_cxxflags)
Chesv_example_qd_LDFLAGS = $(qd_libdepends)

Csysv_example_qd_SOURCES = Csysv_example_qd.cpp
Csysv_example_qd_CXXFLAGS = $(qd_cxxflags)
Csysv_example_qd_LDFLAGS = $(qd_libdepends)

Rsysv_example_qd_SOURCES = Rsysv_example_qd.cpp
Rsysv_example_qd_CXXFLAGS = $(qd_cxxflags)
Rsysv_example_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/06_SymmetricIndefiniteLinearEquations

mplapackexamples_DATA = Chesv_example_mpfr.cpp Chesv_example_gmp.cpp Chesv_example_binary128.cpp Chesv_example_binary80.cpp Chesv_example_double.cpp Chesv_example_dd.cpp Chesv_example_qd.cpp Csysv_example_mpfr.cpp Csysv_example_gmp.cpp Csysv_example_binary128.cpp Csysv_example_binary80.cpp Csysv_example_double.cpp Csysv_example_dd.cpp Csysv_example_qd.cpp Rsysv_example_mpfr.cpp Rsysv_example_gmp.cpp Rsysv_example_binary128.cpp Rsysv_example_binary80.cpp Rsysv_example_double.cpp Rsysv_example_dd.cpp Rsysv_example_qd.cpp \
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
