mpblasexamples_PROGRAMS =

if ENABLE_MPFR
mpblasexamples_PROGRAMS += Cgemm_mpfr Rgemm_mpfr Cgemm_mpfr_opt Rgemm_mpfr_opt

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include $(GMPFRXX_MKII_CPPFLAGS) -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp
mpfr_opt_libdepends = -L$(top_builddir)/mpblas/optimized/mpfr -lmplapack_mpfr_opt -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Cgemm_mpfr_SOURCES = Cgemm_mpfr.cpp
Cgemm_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Cgemm_mpfr_LDFLAGS = $(mpfr_libdepends)

Rgemm_mpfr_SOURCES = Rgemm_mpfr.cpp
Rgemm_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rgemm_mpfr_LDFLAGS = $(mpfr_libdepends)

Cgemm_mpfr_opt_SOURCES = Cgemm_mpfr.cpp
Cgemm_mpfr_opt_CXXFLAGS = $(mpfr_cxxflags)
Cgemm_mpfr_opt_LDFLAGS = $(mpfr_opt_libdepends)

Rgemm_mpfr_opt_SOURCES = Rgemm_mpfr.cpp
Rgemm_mpfr_opt_CXXFLAGS = $(mpfr_cxxflags)
Rgemm_mpfr_opt_LDFLAGS = $(mpfr_opt_libdepends)

endif

if ENABLE_GMP
mpblasexamples_PROGRAMS += Cgemm_gmp Rgemm_gmp Cgemm_gmp_opt Rgemm_gmp_opt

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp
gmp_opt_libdepends = -L$(top_builddir)/mpblas/optimized/gmp -lmplapack_gmp_opt -L$(GMP_LIBDIR) -lgmp

Cgemm_gmp_SOURCES = Cgemm_gmp.cpp
Cgemm_gmp_CXXFLAGS = $(gmp_cxxflags)
Cgemm_gmp_LDFLAGS = $(gmp_libdepends)

Rgemm_gmp_SOURCES = Rgemm_gmp.cpp
Rgemm_gmp_CXXFLAGS = $(gmp_cxxflags)
Rgemm_gmp_LDFLAGS = $(gmp_libdepends)

Cgemm_gmp_opt_SOURCES = Cgemm_gmp.cpp
Cgemm_gmp_opt_CXXFLAGS = $(gmp_cxxflags)
Cgemm_gmp_opt_LDFLAGS = $(gmp_opt_libdepends)

Rgemm_gmp_opt_SOURCES = Rgemm_gmp.cpp
Rgemm_gmp_opt_CXXFLAGS = $(gmp_cxxflags)
Rgemm_gmp_opt_LDFLAGS = $(gmp_opt_libdepends)

endif

if ENABLE_BINARY128
mpblasexamples_PROGRAMS += Cgemm_binary128 Rgemm_binary128 Cgemm_binary128_opt Rgemm_binary128_opt

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_binary128
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends += -lquadmath
endif
binary128_opt_libdepends = -L$(top_builddir)/mpblas/optimized/binary128 -lmplapack_binary128_opt
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_opt_libdepends += -lquadmath
endif

Cgemm_binary128_SOURCES = Cgemm_binary128.cpp
Cgemm_binary128_CXXFLAGS = $(binary128_cxxflags)
Cgemm_binary128_LDFLAGS = $(binary128_libdepends)

Rgemm_binary128_SOURCES = Rgemm_binary128.cpp
Rgemm_binary128_CXXFLAGS = $(binary128_cxxflags)
Rgemm_binary128_LDFLAGS = $(binary128_libdepends)

Cgemm_binary128_opt_SOURCES = Cgemm_binary128.cpp
Cgemm_binary128_opt_CXXFLAGS = $(binary128_cxxflags)
Cgemm_binary128_opt_LDFLAGS = $(binary128_opt_libdepends)

Rgemm_binary128_opt_SOURCES = Rgemm_binary128.cpp
Rgemm_binary128_opt_CXXFLAGS = $(binary128_cxxflags)
Rgemm_binary128_opt_LDFLAGS = $(binary128_opt_libdepends)

endif

if ENABLE_BINARY80
mpblasexamples_PROGRAMS += Cgemm_binary80 Rgemm_binary80 Cgemm_binary80_opt Rgemm_binary80_opt

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_binary80
binary80_opt_libdepends = -L$(top_builddir)/mpblas/optimized/binary80 -lmplapack_binary80_opt

Cgemm_binary80_SOURCES = Cgemm_binary80.cpp
Cgemm_binary80_CXXFLAGS = $(binary80_cxxflags)
Cgemm_binary80_LDFLAGS = $(binary80_libdepends)

Rgemm_binary80_SOURCES = Rgemm_binary80.cpp
Rgemm_binary80_CXXFLAGS = $(binary80_cxxflags)
Rgemm_binary80_LDFLAGS = $(binary80_libdepends)

Cgemm_binary80_opt_SOURCES = Cgemm_binary80.cpp
Cgemm_binary80_opt_CXXFLAGS = $(binary80_cxxflags)
Cgemm_binary80_opt_LDFLAGS = $(binary80_opt_libdepends)

Rgemm_binary80_opt_SOURCES = Rgemm_binary80.cpp
Rgemm_binary80_opt_CXXFLAGS = $(binary80_cxxflags)
Rgemm_binary80_opt_LDFLAGS = $(binary80_opt_libdepends)

endif

if ENABLE_DOUBLE
mpblasexamples_PROGRAMS += Cgemm_double Rgemm_double Cgemm_double_opt Rgemm_double_opt

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_double
double_opt_libdepends = -L$(top_builddir)/mpblas/optimized/double -lmplapack_double_opt

Cgemm_double_SOURCES = Cgemm_double.cpp
Cgemm_double_CXXFLAGS = $(double_cxxflags)
Cgemm_double_LDFLAGS = $(double_libdepends)

Rgemm_double_SOURCES = Rgemm_double.cpp
Rgemm_double_CXXFLAGS = $(double_cxxflags)
Rgemm_double_LDFLAGS = $(double_libdepends)

Cgemm_double_opt_SOURCES = Cgemm_double.cpp
Cgemm_double_opt_CXXFLAGS = $(double_cxxflags)
Cgemm_double_opt_LDFLAGS = $(double_opt_libdepends)

Rgemm_double_opt_SOURCES = Rgemm_double.cpp
Rgemm_double_opt_CXXFLAGS = $(double_cxxflags)
Rgemm_double_opt_LDFLAGS = $(double_opt_libdepends)

endif

if ENABLE_DD
mpblasexamples_PROGRAMS += Cgemm_dd Rgemm_dd Cgemm_dd_opt Rgemm_dd_opt

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd
dd_opt_libdepends = -L$(top_builddir)/mpblas/optimized/dd -lmplapack_dd_opt -L$(QD_LIBDIR) -lqd

Cgemm_dd_SOURCES = Cgemm_dd.cpp
Cgemm_dd_CXXFLAGS = $(dd_cxxflags)
Cgemm_dd_LDFLAGS = $(dd_libdepends)

Rgemm_dd_SOURCES = Rgemm_dd.cpp
Rgemm_dd_CXXFLAGS = $(dd_cxxflags)
Rgemm_dd_LDFLAGS = $(dd_libdepends)

Cgemm_dd_opt_SOURCES = Cgemm_dd.cpp
Cgemm_dd_opt_CXXFLAGS = $(dd_cxxflags)
Cgemm_dd_opt_LDFLAGS = $(dd_opt_libdepends)

Rgemm_dd_opt_SOURCES = Rgemm_dd.cpp
Rgemm_dd_opt_CXXFLAGS = $(dd_cxxflags)
Rgemm_dd_opt_LDFLAGS = $(dd_opt_libdepends)

endif

if ENABLE_QD
mpblasexamples_PROGRAMS += Cgemm_qd Rgemm_qd Cgemm_qd_opt Rgemm_qd_opt

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd
qd_opt_libdepends = -L$(top_builddir)/mpblas/optimized/qd -lmplapack_qd_opt -L$(QD_LIBDIR) -lqd

Cgemm_qd_SOURCES = Cgemm_qd.cpp
Cgemm_qd_CXXFLAGS = $(qd_cxxflags)
Cgemm_qd_LDFLAGS = $(qd_libdepends)

Rgemm_qd_SOURCES = Rgemm_qd.cpp
Rgemm_qd_CXXFLAGS = $(qd_cxxflags)
Rgemm_qd_LDFLAGS = $(qd_libdepends)

Cgemm_qd_opt_SOURCES = Cgemm_qd.cpp
Cgemm_qd_opt_CXXFLAGS = $(qd_cxxflags)
Cgemm_qd_opt_LDFLAGS = $(qd_opt_libdepends)

Rgemm_qd_opt_SOURCES = Rgemm_qd.cpp
Rgemm_qd_opt_CXXFLAGS = $(qd_cxxflags)
Rgemm_qd_opt_LDFLAGS = $(qd_opt_libdepends)

endif

mpblasexamplesdir=$(prefix)/share/examples/mpblas/

mpblasexamples_DATA = Cgemm_mpfr.cpp Cgemm_gmp.cpp Cgemm_binary128.cpp Cgemm_binary80.cpp Cgemm_double.cpp Cgemm_dd.cpp Cgemm_qd.cpp Rgemm_mpfr.cpp Rgemm_gmp.cpp Rgemm_binary128.cpp Rgemm_binary80.cpp Rgemm_double.cpp Rgemm_dd.cpp Rgemm_qd.cpp \
Makefile.freebsd Makefile.linux Makefile.linux.inteloneAPI Makefile.linux_cuda Makefile.macos Makefile.mingw

install-data-hook:
if IS_MACOS
	bash $(top_srcdir)/misc/fix_dylib_macOS.sh $(mpblasexamplesdir) $(prefix)
endif

EXTRA_DIST = \
	Makefile.freebsd.in \
	Makefile.linux.in \
	Makefile.linux.inteloneAPI.in \
	Makefile.linux_cuda.in \
	Makefile.macos.in \
	Makefile.mingw.in

GENERATED_MAKEFILES = \
	Makefile.freebsd \
	Makefile.linux \
	Makefile.linux.inteloneAPI \
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
