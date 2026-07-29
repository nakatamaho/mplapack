mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cggev_example_mpfr Rgges_example_mpfr Rggev_example_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cggev_example_mpfr_SOURCES = Cggev_example_mpfr.cpp
Cggev_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cggev_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgges_example_mpfr_SOURCES = Rgges_example_mpfr.cpp
Rgges_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgges_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rggev_example_mpfr_SOURCES = Rggev_example_mpfr.cpp
Rggev_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rggev_example_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cggev_example_gmp Rgges_example_gmp Rggev_example_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Cggev_example_gmp_SOURCES = Cggev_example_gmp.cpp
Cggev_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cggev_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgges_example_gmp_SOURCES = Rgges_example_gmp.cpp
Rgges_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgges_example_gmp_LDFLAGS = $(gmp_libdepends)

Rggev_example_gmp_SOURCES = Rggev_example_gmp.cpp
Rggev_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rggev_example_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cggev_example_binary128 Rgges_example_binary128 Rggev_example_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Cggev_example_binary128_SOURCES = Cggev_example_binary128.cpp
Cggev_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cggev_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgges_example_binary128_SOURCES = Rgges_example_binary128.cpp
Rgges_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgges_example_binary128_LDFLAGS = $(binary128_libdepends)

Rggev_example_binary128_SOURCES = Rggev_example_binary128.cpp
Rggev_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rggev_example_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cggev_example_binary80 Rgges_example_binary80 Rggev_example_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Cggev_example_binary80_SOURCES = Cggev_example_binary80.cpp
Cggev_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cggev_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgges_example_binary80_SOURCES = Rgges_example_binary80.cpp
Rgges_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgges_example_binary80_LDFLAGS = $(binary80_libdepends)

Rggev_example_binary80_SOURCES = Rggev_example_binary80.cpp
Rggev_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rggev_example_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cggev_example_double Rgges_example_double Rggev_example_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Cggev_example_double_SOURCES = Cggev_example_double.cpp
Cggev_example_double_CXXFLAGS = $(double_cxxflags)
Cggev_example_double_LDFLAGS = $(double_libdepends)

Rgges_example_double_SOURCES = Rgges_example_double.cpp
Rgges_example_double_CXXFLAGS = $(double_cxxflags)
Rgges_example_double_LDFLAGS = $(double_libdepends)

Rggev_example_double_SOURCES = Rggev_example_double.cpp
Rggev_example_double_CXXFLAGS = $(double_cxxflags)
Rggev_example_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cggev_example_dd Rgges_example_dd Rggev_example_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Cggev_example_dd_SOURCES = Cggev_example_dd.cpp
Cggev_example_dd_CXXFLAGS = $(dd_cxxflags)
Cggev_example_dd_LDFLAGS = $(dd_libdepends)

Rgges_example_dd_SOURCES = Rgges_example_dd.cpp
Rgges_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgges_example_dd_LDFLAGS = $(dd_libdepends)

Rggev_example_dd_SOURCES = Rggev_example_dd.cpp
Rggev_example_dd_CXXFLAGS = $(dd_cxxflags)
Rggev_example_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cggev_example_qd Rgges_example_qd Rggev_example_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Cggev_example_qd_SOURCES = Cggev_example_qd.cpp
Cggev_example_qd_CXXFLAGS = $(qd_cxxflags)
Cggev_example_qd_LDFLAGS = $(qd_libdepends)

Rgges_example_qd_SOURCES = Rgges_example_qd.cpp
Rgges_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgges_example_qd_LDFLAGS = $(qd_libdepends)

Rggev_example_qd_SOURCES = Rggev_example_qd.cpp
Rggev_example_qd_CXXFLAGS = $(qd_cxxflags)
Rggev_example_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/08_GeneralizedNonsymmetricEigenproblems

mplapackexamples_DATA = Cggev_example_mpfr.cpp Cggev_example_gmp.cpp Cggev_example_binary128.cpp Cggev_example_binary80.cpp Cggev_example_double.cpp Cggev_example_dd.cpp Cggev_example_qd.cpp Rgges_example_mpfr.cpp Rgges_example_gmp.cpp Rgges_example_binary128.cpp Rgges_example_binary80.cpp Rgges_example_double.cpp Rgges_example_dd.cpp Rgges_example_qd.cpp Rggev_example_mpfr.cpp Rggev_example_gmp.cpp Rggev_example_binary128.cpp Rggev_example_binary80.cpp Rggev_example_double.cpp Rggev_example_dd.cpp Rggev_example_qd.cpp \
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
