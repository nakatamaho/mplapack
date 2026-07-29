mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Cgels_example_mpfr Cgelss_example_mpfr Rgels_example_mpfr Rgelss_example_mpfr Rgelsy_example_mpfr Rggglm_example_mpfr Rgglse_example_mpfr

mpfr_cxxflags = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgels_example_mpfr_SOURCES = Cgels_example_mpfr.cpp
Cgels_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgels_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgelss_example_mpfr_SOURCES = Cgelss_example_mpfr.cpp
Cgelss_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgelss_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgels_example_mpfr_SOURCES = Rgels_example_mpfr.cpp
Rgels_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgels_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgelss_example_mpfr_SOURCES = Rgelss_example_mpfr.cpp
Rgelss_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgelss_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgelsy_example_mpfr_SOURCES = Rgelsy_example_mpfr.cpp
Rgelsy_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgelsy_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rggglm_example_mpfr_SOURCES = Rggglm_example_mpfr.cpp
Rggglm_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rggglm_example_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgglse_example_mpfr_SOURCES = Rgglse_example_mpfr.cpp
Rgglse_example_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgglse_example_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Cgels_example_gmp Cgelss_example_gmp Rgels_example_gmp Rgelss_example_gmp Rgelsy_example_gmp Rggglm_example_gmp Rgglse_example_gmp

gmp_cxxflags = $(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Cgels_example_gmp_SOURCES = Cgels_example_gmp.cpp
Cgels_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgels_example_gmp_LDFLAGS = $(gmp_libdepends)

Cgelss_example_gmp_SOURCES = Cgelss_example_gmp.cpp
Cgelss_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgelss_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgels_example_gmp_SOURCES = Rgels_example_gmp.cpp
Rgels_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgels_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgelss_example_gmp_SOURCES = Rgelss_example_gmp.cpp
Rgelss_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgelss_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgelsy_example_gmp_SOURCES = Rgelsy_example_gmp.cpp
Rgelsy_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgelsy_example_gmp_LDFLAGS = $(gmp_libdepends)

Rggglm_example_gmp_SOURCES = Rggglm_example_gmp.cpp
Rggglm_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rggglm_example_gmp_LDFLAGS = $(gmp_libdepends)

Rgglse_example_gmp_SOURCES = Rgglse_example_gmp.cpp
Rgglse_example_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgglse_example_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Cgels_example_binary128 Cgelss_example_binary128 Rgels_example_binary128 Rgelss_example_binary128 Rgelsy_example_binary128 Rggglm_example_binary128 Rgglse_example_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Cgels_example_binary128_SOURCES = Cgels_example_binary128.cpp
Cgels_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgels_example_binary128_LDFLAGS = $(binary128_libdepends)

Cgelss_example_binary128_SOURCES = Cgelss_example_binary128.cpp
Cgelss_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgelss_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgels_example_binary128_SOURCES = Rgels_example_binary128.cpp
Rgels_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgels_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgelss_example_binary128_SOURCES = Rgelss_example_binary128.cpp
Rgelss_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgelss_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgelsy_example_binary128_SOURCES = Rgelsy_example_binary128.cpp
Rgelsy_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgelsy_example_binary128_LDFLAGS = $(binary128_libdepends)

Rggglm_example_binary128_SOURCES = Rggglm_example_binary128.cpp
Rggglm_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rggglm_example_binary128_LDFLAGS = $(binary128_libdepends)

Rgglse_example_binary128_SOURCES = Rgglse_example_binary128.cpp
Rgglse_example_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgglse_example_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Cgels_example_binary80 Cgelss_example_binary80 Rgels_example_binary80 Rgelss_example_binary80 Rgelsy_example_binary80 Rggglm_example_binary80 Rgglse_example_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Cgels_example_binary80_SOURCES = Cgels_example_binary80.cpp
Cgels_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgels_example_binary80_LDFLAGS = $(binary80_libdepends)

Cgelss_example_binary80_SOURCES = Cgelss_example_binary80.cpp
Cgelss_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgelss_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgels_example_binary80_SOURCES = Rgels_example_binary80.cpp
Rgels_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgels_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgelss_example_binary80_SOURCES = Rgelss_example_binary80.cpp
Rgelss_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgelss_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgelsy_example_binary80_SOURCES = Rgelsy_example_binary80.cpp
Rgelsy_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgelsy_example_binary80_LDFLAGS = $(binary80_libdepends)

Rggglm_example_binary80_SOURCES = Rggglm_example_binary80.cpp
Rggglm_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rggglm_example_binary80_LDFLAGS = $(binary80_libdepends)

Rgglse_example_binary80_SOURCES = Rgglse_example_binary80.cpp
Rgglse_example_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgglse_example_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Cgels_example_double Cgelss_example_double Rgels_example_double Rgelss_example_double Rgelsy_example_double Rggglm_example_double Rgglse_example_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Cgels_example_double_SOURCES = Cgels_example_double.cpp
Cgels_example_double_CXXFLAGS = $(double_cxxflags)
Cgels_example_double_LDFLAGS = $(double_libdepends)

Cgelss_example_double_SOURCES = Cgelss_example_double.cpp
Cgelss_example_double_CXXFLAGS = $(double_cxxflags)
Cgelss_example_double_LDFLAGS = $(double_libdepends)

Rgels_example_double_SOURCES = Rgels_example_double.cpp
Rgels_example_double_CXXFLAGS = $(double_cxxflags)
Rgels_example_double_LDFLAGS = $(double_libdepends)

Rgelss_example_double_SOURCES = Rgelss_example_double.cpp
Rgelss_example_double_CXXFLAGS = $(double_cxxflags)
Rgelss_example_double_LDFLAGS = $(double_libdepends)

Rgelsy_example_double_SOURCES = Rgelsy_example_double.cpp
Rgelsy_example_double_CXXFLAGS = $(double_cxxflags)
Rgelsy_example_double_LDFLAGS = $(double_libdepends)

Rggglm_example_double_SOURCES = Rggglm_example_double.cpp
Rggglm_example_double_CXXFLAGS = $(double_cxxflags)
Rggglm_example_double_LDFLAGS = $(double_libdepends)

Rgglse_example_double_SOURCES = Rgglse_example_double.cpp
Rgglse_example_double_CXXFLAGS = $(double_cxxflags)
Rgglse_example_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Cgels_example_dd Cgelss_example_dd Rgels_example_dd Rgelss_example_dd Rgelsy_example_dd Rggglm_example_dd Rgglse_example_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Cgels_example_dd_SOURCES = Cgels_example_dd.cpp
Cgels_example_dd_CXXFLAGS = $(dd_cxxflags)
Cgels_example_dd_LDFLAGS = $(dd_libdepends)

Cgelss_example_dd_SOURCES = Cgelss_example_dd.cpp
Cgelss_example_dd_CXXFLAGS = $(dd_cxxflags)
Cgelss_example_dd_LDFLAGS = $(dd_libdepends)

Rgels_example_dd_SOURCES = Rgels_example_dd.cpp
Rgels_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgels_example_dd_LDFLAGS = $(dd_libdepends)

Rgelss_example_dd_SOURCES = Rgelss_example_dd.cpp
Rgelss_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgelss_example_dd_LDFLAGS = $(dd_libdepends)

Rgelsy_example_dd_SOURCES = Rgelsy_example_dd.cpp
Rgelsy_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgelsy_example_dd_LDFLAGS = $(dd_libdepends)

Rggglm_example_dd_SOURCES = Rggglm_example_dd.cpp
Rggglm_example_dd_CXXFLAGS = $(dd_cxxflags)
Rggglm_example_dd_LDFLAGS = $(dd_libdepends)

Rgglse_example_dd_SOURCES = Rgglse_example_dd.cpp
Rgglse_example_dd_CXXFLAGS = $(dd_cxxflags)
Rgglse_example_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Cgels_example_qd Cgelss_example_qd Rgels_example_qd Rgelss_example_qd Rgelsy_example_qd Rggglm_example_qd Rgglse_example_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Cgels_example_qd_SOURCES = Cgels_example_qd.cpp
Cgels_example_qd_CXXFLAGS = $(qd_cxxflags)
Cgels_example_qd_LDFLAGS = $(qd_libdepends)

Cgelss_example_qd_SOURCES = Cgelss_example_qd.cpp
Cgelss_example_qd_CXXFLAGS = $(qd_cxxflags)
Cgelss_example_qd_LDFLAGS = $(qd_libdepends)

Rgels_example_qd_SOURCES = Rgels_example_qd.cpp
Rgels_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgels_example_qd_LDFLAGS = $(qd_libdepends)

Rgelss_example_qd_SOURCES = Rgelss_example_qd.cpp
Rgelss_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgelss_example_qd_LDFLAGS = $(qd_libdepends)

Rgelsy_example_qd_SOURCES = Rgelsy_example_qd.cpp
Rgelsy_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgelsy_example_qd_LDFLAGS = $(qd_libdepends)

Rggglm_example_qd_SOURCES = Rggglm_example_qd.cpp
Rggglm_example_qd_CXXFLAGS = $(qd_cxxflags)
Rggglm_example_qd_LDFLAGS = $(qd_libdepends)

Rgglse_example_qd_SOURCES = Rgglse_example_qd.cpp
Rgglse_example_qd_CXXFLAGS = $(qd_cxxflags)
Rgglse_example_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/02_LeastSquares

mplapackexamples_DATA = Cgels_example_mpfr.cpp Cgels_example_gmp.cpp Cgels_example_binary128.cpp Cgels_example_binary80.cpp Cgels_example_double.cpp Cgels_example_dd.cpp Cgels_example_qd.cpp Cgelss_example_mpfr.cpp Cgelss_example_gmp.cpp Cgelss_example_binary128.cpp Cgelss_example_binary80.cpp Cgelss_example_double.cpp Cgelss_example_dd.cpp Cgelss_example_qd.cpp Rgels_example_mpfr.cpp Rgels_example_gmp.cpp Rgels_example_binary128.cpp Rgels_example_binary80.cpp Rgels_example_double.cpp Rgels_example_dd.cpp Rgels_example_qd.cpp Rgelss_example_mpfr.cpp Rgelss_example_gmp.cpp Rgelss_example_binary128.cpp Rgelss_example_binary80.cpp Rgelss_example_double.cpp Rgelss_example_dd.cpp Rgelss_example_qd.cpp Rgelsy_example_mpfr.cpp Rgelsy_example_gmp.cpp Rgelsy_example_binary128.cpp Rgelsy_example_binary80.cpp Rgelsy_example_double.cpp Rgelsy_example_dd.cpp Rgelsy_example_qd.cpp Rggglm_example_mpfr.cpp Rggglm_example_gmp.cpp Rggglm_example_binary128.cpp Rggglm_example_binary80.cpp Rggglm_example_double.cpp Rggglm_example_dd.cpp Rggglm_example_qd.cpp Rgglse_example_mpfr.cpp Rgglse_example_gmp.cpp Rgglse_example_binary128.cpp Rgglse_example_binary80.cpp Rgglse_example_double.cpp Rgglse_example_dd.cpp Rgglse_example_qd.cpp \
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
