mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Hilbert_inverse_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Hilbert_inverse_mpfr_SOURCES = Hilbert_inverse_mpfr.cpp
Hilbert_inverse_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Hilbert_inverse_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Hilbert_inverse_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Hilbert_inverse_gmp_SOURCES = Hilbert_inverse_gmp.cpp
Hilbert_inverse_gmp_CXXFLAGS = $(gmp_cxxflags)
Hilbert_inverse_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Hilbert_inverse_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Hilbert_inverse_binary128_SOURCES = Hilbert_inverse_binary128.cpp
Hilbert_inverse_binary128_CXXFLAGS = $(binary128_cxxflags)
Hilbert_inverse_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Hilbert_inverse_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Hilbert_inverse_binary80_SOURCES = Hilbert_inverse_binary80.cpp
Hilbert_inverse_binary80_CXXFLAGS = $(binary80_cxxflags)
Hilbert_inverse_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Hilbert_inverse_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Hilbert_inverse_double_SOURCES = Hilbert_inverse_double.cpp
Hilbert_inverse_double_CXXFLAGS = $(double_cxxflags)
Hilbert_inverse_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Hilbert_inverse_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Hilbert_inverse_dd_SOURCES = Hilbert_inverse_dd.cpp
Hilbert_inverse_dd_CXXFLAGS = $(dd_cxxflags)
Hilbert_inverse_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Hilbert_inverse_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Hilbert_inverse_qd_SOURCES = Hilbert_inverse_qd.cpp
Hilbert_inverse_qd_CXXFLAGS = $(qd_cxxflags)
Hilbert_inverse_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/90_PrecisionComparison

mplapackexamples_DATA = Hilbert_inverse_mpfr.cpp Hilbert_inverse_gmp.cpp Hilbert_inverse_binary128.cpp Hilbert_inverse_binary80.cpp Hilbert_inverse_double.cpp Hilbert_inverse_dd.cpp Hilbert_inverse_qd.cpp \
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
