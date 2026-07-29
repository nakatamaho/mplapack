mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += Rlamch_dump_mpfr Frank_eigenvalues_mpfr Gejsv_vs_gesvd_mpfr Hilbert_inverse_mpfr Hilbert_linear_solve_mpfr Kahan_rank_reveal_mpfr Pascal_inverse_mpfr Vandermonde_solve_mpfr Wilkinson_roots_mpfr

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

Rlamch_dump_mpfr_SOURCES = Rlamch_dump_mpfr.cpp
Rlamch_dump_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Rlamch_dump_mpfr_LDFLAGS = $(mpfr_libdepends)

Frank_eigenvalues_mpfr_SOURCES = Frank_eigenvalues_mpfr.cpp
Frank_eigenvalues_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Frank_eigenvalues_mpfr_LDFLAGS = $(mpfr_libdepends)

Gejsv_vs_gesvd_mpfr_SOURCES = Gejsv_vs_gesvd_mpfr.cpp
Gejsv_vs_gesvd_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Gejsv_vs_gesvd_mpfr_LDFLAGS = $(mpfr_libdepends)

Hilbert_inverse_mpfr_SOURCES = Hilbert_inverse_mpfr.cpp
Hilbert_inverse_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Hilbert_inverse_mpfr_LDFLAGS = $(mpfr_libdepends)

Hilbert_linear_solve_mpfr_SOURCES = Hilbert_linear_solve_mpfr.cpp
Hilbert_linear_solve_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Hilbert_linear_solve_mpfr_LDFLAGS = $(mpfr_libdepends)

Kahan_rank_reveal_mpfr_SOURCES = Kahan_rank_reveal_mpfr.cpp
Kahan_rank_reveal_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Kahan_rank_reveal_mpfr_LDFLAGS = $(mpfr_libdepends)

Pascal_inverse_mpfr_SOURCES = Pascal_inverse_mpfr.cpp
Pascal_inverse_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Pascal_inverse_mpfr_LDFLAGS = $(mpfr_libdepends)

Vandermonde_solve_mpfr_SOURCES = Vandermonde_solve_mpfr.cpp
Vandermonde_solve_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Vandermonde_solve_mpfr_LDFLAGS = $(mpfr_libdepends)

Wilkinson_roots_mpfr_SOURCES = Wilkinson_roots_mpfr.cpp
Wilkinson_roots_mpfr_CXXFLAGS = $(mpfr_cxxflags)
Wilkinson_roots_mpfr_LDFLAGS = $(mpfr_libdepends)

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += Rlamch_dump_gmp Frank_eigenvalues_gmp Gejsv_vs_gesvd_gmp Hilbert_inverse_gmp Hilbert_linear_solve_gmp Kahan_rank_reveal_gmp Pascal_inverse_gmp Vandermonde_solve_gmp Wilkinson_roots_gmp

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(GMP_LIBDIR) -lgmp

Rlamch_dump_gmp_SOURCES = Rlamch_dump_gmp.cpp
Rlamch_dump_gmp_CXXFLAGS = $(gmp_cxxflags)
Rlamch_dump_gmp_LDFLAGS = $(gmp_libdepends)

Frank_eigenvalues_gmp_SOURCES = Frank_eigenvalues_gmp.cpp
Frank_eigenvalues_gmp_CXXFLAGS = $(gmp_cxxflags)
Frank_eigenvalues_gmp_LDFLAGS = $(gmp_libdepends)

Gejsv_vs_gesvd_gmp_SOURCES = Gejsv_vs_gesvd_gmp.cpp
Gejsv_vs_gesvd_gmp_CXXFLAGS = $(gmp_cxxflags)
Gejsv_vs_gesvd_gmp_LDFLAGS = $(gmp_libdepends)

Hilbert_inverse_gmp_SOURCES = Hilbert_inverse_gmp.cpp
Hilbert_inverse_gmp_CXXFLAGS = $(gmp_cxxflags)
Hilbert_inverse_gmp_LDFLAGS = $(gmp_libdepends)

Hilbert_linear_solve_gmp_SOURCES = Hilbert_linear_solve_gmp.cpp
Hilbert_linear_solve_gmp_CXXFLAGS = $(gmp_cxxflags)
Hilbert_linear_solve_gmp_LDFLAGS = $(gmp_libdepends)

Kahan_rank_reveal_gmp_SOURCES = Kahan_rank_reveal_gmp.cpp
Kahan_rank_reveal_gmp_CXXFLAGS = $(gmp_cxxflags)
Kahan_rank_reveal_gmp_LDFLAGS = $(gmp_libdepends)

Pascal_inverse_gmp_SOURCES = Pascal_inverse_gmp.cpp
Pascal_inverse_gmp_CXXFLAGS = $(gmp_cxxflags)
Pascal_inverse_gmp_LDFLAGS = $(gmp_libdepends)

Vandermonde_solve_gmp_SOURCES = Vandermonde_solve_gmp.cpp
Vandermonde_solve_gmp_CXXFLAGS = $(gmp_cxxflags)
Vandermonde_solve_gmp_LDFLAGS = $(gmp_libdepends)

Wilkinson_roots_gmp_SOURCES = Wilkinson_roots_gmp.cpp
Wilkinson_roots_gmp_CXXFLAGS = $(gmp_cxxflags)
Wilkinson_roots_gmp_LDFLAGS = $(gmp_libdepends)

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += Rlamch_dump_binary128 Frank_eigenvalues_binary128 Gejsv_vs_gesvd_binary128 Hilbert_inverse_binary128 Hilbert_linear_solve_binary128 Kahan_rank_reveal_binary128 Pascal_inverse_binary128 Vandermonde_solve_binary128 Wilkinson_roots_binary128

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128
endif

Rlamch_dump_binary128_SOURCES = Rlamch_dump_binary128.cpp
Rlamch_dump_binary128_CXXFLAGS = $(binary128_cxxflags)
Rlamch_dump_binary128_LDFLAGS = $(binary128_libdepends)

Frank_eigenvalues_binary128_SOURCES = Frank_eigenvalues_binary128.cpp
Frank_eigenvalues_binary128_CXXFLAGS = $(binary128_cxxflags)
Frank_eigenvalues_binary128_LDFLAGS = $(binary128_libdepends)

Gejsv_vs_gesvd_binary128_SOURCES = Gejsv_vs_gesvd_binary128.cpp
Gejsv_vs_gesvd_binary128_CXXFLAGS = $(binary128_cxxflags)
Gejsv_vs_gesvd_binary128_LDFLAGS = $(binary128_libdepends)

Hilbert_inverse_binary128_SOURCES = Hilbert_inverse_binary128.cpp
Hilbert_inverse_binary128_CXXFLAGS = $(binary128_cxxflags)
Hilbert_inverse_binary128_LDFLAGS = $(binary128_libdepends)

Hilbert_linear_solve_binary128_SOURCES = Hilbert_linear_solve_binary128.cpp
Hilbert_linear_solve_binary128_CXXFLAGS = $(binary128_cxxflags)
Hilbert_linear_solve_binary128_LDFLAGS = $(binary128_libdepends)

Kahan_rank_reveal_binary128_SOURCES = Kahan_rank_reveal_binary128.cpp
Kahan_rank_reveal_binary128_CXXFLAGS = $(binary128_cxxflags)
Kahan_rank_reveal_binary128_LDFLAGS = $(binary128_libdepends)

Pascal_inverse_binary128_SOURCES = Pascal_inverse_binary128.cpp
Pascal_inverse_binary128_CXXFLAGS = $(binary128_cxxflags)
Pascal_inverse_binary128_LDFLAGS = $(binary128_libdepends)

Vandermonde_solve_binary128_SOURCES = Vandermonde_solve_binary128.cpp
Vandermonde_solve_binary128_CXXFLAGS = $(binary128_cxxflags)
Vandermonde_solve_binary128_LDFLAGS = $(binary128_libdepends)

Wilkinson_roots_binary128_SOURCES = Wilkinson_roots_binary128.cpp
Wilkinson_roots_binary128_CXXFLAGS = $(binary128_cxxflags)
Wilkinson_roots_binary128_LDFLAGS = $(binary128_libdepends)

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += Rlamch_dump_binary80 Frank_eigenvalues_binary80 Gejsv_vs_gesvd_binary80 Hilbert_inverse_binary80 Hilbert_linear_solve_binary80 Kahan_rank_reveal_binary80 Pascal_inverse_binary80 Vandermonde_solve_binary80 Wilkinson_roots_binary80

binary80_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80

Rlamch_dump_binary80_SOURCES = Rlamch_dump_binary80.cpp
Rlamch_dump_binary80_CXXFLAGS = $(binary80_cxxflags)
Rlamch_dump_binary80_LDFLAGS = $(binary80_libdepends)

Frank_eigenvalues_binary80_SOURCES = Frank_eigenvalues_binary80.cpp
Frank_eigenvalues_binary80_CXXFLAGS = $(binary80_cxxflags)
Frank_eigenvalues_binary80_LDFLAGS = $(binary80_libdepends)

Gejsv_vs_gesvd_binary80_SOURCES = Gejsv_vs_gesvd_binary80.cpp
Gejsv_vs_gesvd_binary80_CXXFLAGS = $(binary80_cxxflags)
Gejsv_vs_gesvd_binary80_LDFLAGS = $(binary80_libdepends)

Hilbert_inverse_binary80_SOURCES = Hilbert_inverse_binary80.cpp
Hilbert_inverse_binary80_CXXFLAGS = $(binary80_cxxflags)
Hilbert_inverse_binary80_LDFLAGS = $(binary80_libdepends)

Hilbert_linear_solve_binary80_SOURCES = Hilbert_linear_solve_binary80.cpp
Hilbert_linear_solve_binary80_CXXFLAGS = $(binary80_cxxflags)
Hilbert_linear_solve_binary80_LDFLAGS = $(binary80_libdepends)

Kahan_rank_reveal_binary80_SOURCES = Kahan_rank_reveal_binary80.cpp
Kahan_rank_reveal_binary80_CXXFLAGS = $(binary80_cxxflags)
Kahan_rank_reveal_binary80_LDFLAGS = $(binary80_libdepends)

Pascal_inverse_binary80_SOURCES = Pascal_inverse_binary80.cpp
Pascal_inverse_binary80_CXXFLAGS = $(binary80_cxxflags)
Pascal_inverse_binary80_LDFLAGS = $(binary80_libdepends)

Vandermonde_solve_binary80_SOURCES = Vandermonde_solve_binary80.cpp
Vandermonde_solve_binary80_CXXFLAGS = $(binary80_cxxflags)
Vandermonde_solve_binary80_LDFLAGS = $(binary80_libdepends)

Wilkinson_roots_binary80_SOURCES = Wilkinson_roots_binary80.cpp
Wilkinson_roots_binary80_CXXFLAGS = $(binary80_cxxflags)
Wilkinson_roots_binary80_LDFLAGS = $(binary80_libdepends)

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += Rlamch_dump_double Frank_eigenvalues_double Gejsv_vs_gesvd_double Hilbert_inverse_double Hilbert_linear_solve_double Kahan_rank_reveal_double Pascal_inverse_double Vandermonde_solve_double Wilkinson_roots_double

double_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double

Rlamch_dump_double_SOURCES = Rlamch_dump_double.cpp
Rlamch_dump_double_CXXFLAGS = $(double_cxxflags)
Rlamch_dump_double_LDFLAGS = $(double_libdepends)

Frank_eigenvalues_double_SOURCES = Frank_eigenvalues_double.cpp
Frank_eigenvalues_double_CXXFLAGS = $(double_cxxflags)
Frank_eigenvalues_double_LDFLAGS = $(double_libdepends)

Gejsv_vs_gesvd_double_SOURCES = Gejsv_vs_gesvd_double.cpp
Gejsv_vs_gesvd_double_CXXFLAGS = $(double_cxxflags)
Gejsv_vs_gesvd_double_LDFLAGS = $(double_libdepends)

Hilbert_inverse_double_SOURCES = Hilbert_inverse_double.cpp
Hilbert_inverse_double_CXXFLAGS = $(double_cxxflags)
Hilbert_inverse_double_LDFLAGS = $(double_libdepends)

Hilbert_linear_solve_double_SOURCES = Hilbert_linear_solve_double.cpp
Hilbert_linear_solve_double_CXXFLAGS = $(double_cxxflags)
Hilbert_linear_solve_double_LDFLAGS = $(double_libdepends)

Kahan_rank_reveal_double_SOURCES = Kahan_rank_reveal_double.cpp
Kahan_rank_reveal_double_CXXFLAGS = $(double_cxxflags)
Kahan_rank_reveal_double_LDFLAGS = $(double_libdepends)

Pascal_inverse_double_SOURCES = Pascal_inverse_double.cpp
Pascal_inverse_double_CXXFLAGS = $(double_cxxflags)
Pascal_inverse_double_LDFLAGS = $(double_libdepends)

Vandermonde_solve_double_SOURCES = Vandermonde_solve_double.cpp
Vandermonde_solve_double_CXXFLAGS = $(double_cxxflags)
Vandermonde_solve_double_LDFLAGS = $(double_libdepends)

Wilkinson_roots_double_SOURCES = Wilkinson_roots_double.cpp
Wilkinson_roots_double_CXXFLAGS = $(double_cxxflags)
Wilkinson_roots_double_LDFLAGS = $(double_libdepends)

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += Rlamch_dump_dd Frank_eigenvalues_dd Gejsv_vs_gesvd_dd Hilbert_inverse_dd Hilbert_linear_solve_dd Kahan_rank_reveal_dd Pascal_inverse_dd Vandermonde_solve_dd Wilkinson_roots_dd

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(QD_LIBDIR) -lqd

Rlamch_dump_dd_SOURCES = Rlamch_dump_dd.cpp
Rlamch_dump_dd_CXXFLAGS = $(dd_cxxflags)
Rlamch_dump_dd_LDFLAGS = $(dd_libdepends)

Frank_eigenvalues_dd_SOURCES = Frank_eigenvalues_dd.cpp
Frank_eigenvalues_dd_CXXFLAGS = $(dd_cxxflags)
Frank_eigenvalues_dd_LDFLAGS = $(dd_libdepends)

Gejsv_vs_gesvd_dd_SOURCES = Gejsv_vs_gesvd_dd.cpp
Gejsv_vs_gesvd_dd_CXXFLAGS = $(dd_cxxflags)
Gejsv_vs_gesvd_dd_LDFLAGS = $(dd_libdepends)

Hilbert_inverse_dd_SOURCES = Hilbert_inverse_dd.cpp
Hilbert_inverse_dd_CXXFLAGS = $(dd_cxxflags)
Hilbert_inverse_dd_LDFLAGS = $(dd_libdepends)

Hilbert_linear_solve_dd_SOURCES = Hilbert_linear_solve_dd.cpp
Hilbert_linear_solve_dd_CXXFLAGS = $(dd_cxxflags)
Hilbert_linear_solve_dd_LDFLAGS = $(dd_libdepends)

Kahan_rank_reveal_dd_SOURCES = Kahan_rank_reveal_dd.cpp
Kahan_rank_reveal_dd_CXXFLAGS = $(dd_cxxflags)
Kahan_rank_reveal_dd_LDFLAGS = $(dd_libdepends)

Pascal_inverse_dd_SOURCES = Pascal_inverse_dd.cpp
Pascal_inverse_dd_CXXFLAGS = $(dd_cxxflags)
Pascal_inverse_dd_LDFLAGS = $(dd_libdepends)

Vandermonde_solve_dd_SOURCES = Vandermonde_solve_dd.cpp
Vandermonde_solve_dd_CXXFLAGS = $(dd_cxxflags)
Vandermonde_solve_dd_LDFLAGS = $(dd_libdepends)

Wilkinson_roots_dd_SOURCES = Wilkinson_roots_dd.cpp
Wilkinson_roots_dd_CXXFLAGS = $(dd_cxxflags)
Wilkinson_roots_dd_LDFLAGS = $(dd_libdepends)

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += Rlamch_dump_qd Frank_eigenvalues_qd Gejsv_vs_gesvd_qd Hilbert_inverse_qd Hilbert_linear_solve_qd Kahan_rank_reveal_qd Pascal_inverse_qd Vandermonde_solve_qd Wilkinson_roots_qd

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(QD_LIBDIR) -lqd

Rlamch_dump_qd_SOURCES = Rlamch_dump_qd.cpp
Rlamch_dump_qd_CXXFLAGS = $(qd_cxxflags)
Rlamch_dump_qd_LDFLAGS = $(qd_libdepends)

Frank_eigenvalues_qd_SOURCES = Frank_eigenvalues_qd.cpp
Frank_eigenvalues_qd_CXXFLAGS = $(qd_cxxflags)
Frank_eigenvalues_qd_LDFLAGS = $(qd_libdepends)

Gejsv_vs_gesvd_qd_SOURCES = Gejsv_vs_gesvd_qd.cpp
Gejsv_vs_gesvd_qd_CXXFLAGS = $(qd_cxxflags)
Gejsv_vs_gesvd_qd_LDFLAGS = $(qd_libdepends)

Hilbert_inverse_qd_SOURCES = Hilbert_inverse_qd.cpp
Hilbert_inverse_qd_CXXFLAGS = $(qd_cxxflags)
Hilbert_inverse_qd_LDFLAGS = $(qd_libdepends)

Hilbert_linear_solve_qd_SOURCES = Hilbert_linear_solve_qd.cpp
Hilbert_linear_solve_qd_CXXFLAGS = $(qd_cxxflags)
Hilbert_linear_solve_qd_LDFLAGS = $(qd_libdepends)

Kahan_rank_reveal_qd_SOURCES = Kahan_rank_reveal_qd.cpp
Kahan_rank_reveal_qd_CXXFLAGS = $(qd_cxxflags)
Kahan_rank_reveal_qd_LDFLAGS = $(qd_libdepends)

Pascal_inverse_qd_SOURCES = Pascal_inverse_qd.cpp
Pascal_inverse_qd_CXXFLAGS = $(qd_cxxflags)
Pascal_inverse_qd_LDFLAGS = $(qd_libdepends)

Vandermonde_solve_qd_SOURCES = Vandermonde_solve_qd.cpp
Vandermonde_solve_qd_CXXFLAGS = $(qd_cxxflags)
Vandermonde_solve_qd_LDFLAGS = $(qd_libdepends)

Wilkinson_roots_qd_SOURCES = Wilkinson_roots_qd.cpp
Wilkinson_roots_qd_CXXFLAGS = $(qd_cxxflags)
Wilkinson_roots_qd_LDFLAGS = $(qd_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/90_PrecisionComparison

mplapackexamples_DATA = Rlamch_dump_mpfr.cpp Rlamch_dump_gmp.cpp Rlamch_dump_binary128.cpp Rlamch_dump_binary80.cpp Rlamch_dump_double.cpp Rlamch_dump_dd.cpp Rlamch_dump_qd.cpp Frank_eigenvalues_mpfr.cpp Frank_eigenvalues_gmp.cpp Frank_eigenvalues_binary128.cpp Frank_eigenvalues_binary80.cpp Frank_eigenvalues_double.cpp Frank_eigenvalues_dd.cpp Frank_eigenvalues_qd.cpp Gejsv_vs_gesvd_mpfr.cpp Gejsv_vs_gesvd_gmp.cpp Gejsv_vs_gesvd_binary128.cpp Gejsv_vs_gesvd_binary80.cpp Gejsv_vs_gesvd_double.cpp Gejsv_vs_gesvd_dd.cpp Gejsv_vs_gesvd_qd.cpp Hilbert_inverse_mpfr.cpp Hilbert_inverse_gmp.cpp Hilbert_inverse_binary128.cpp Hilbert_inverse_binary80.cpp Hilbert_inverse_double.cpp Hilbert_inverse_dd.cpp Hilbert_inverse_qd.cpp Hilbert_linear_solve_mpfr.cpp Hilbert_linear_solve_gmp.cpp Hilbert_linear_solve_binary128.cpp Hilbert_linear_solve_binary80.cpp Hilbert_linear_solve_double.cpp Hilbert_linear_solve_dd.cpp Hilbert_linear_solve_qd.cpp Kahan_rank_reveal_mpfr.cpp Kahan_rank_reveal_gmp.cpp Kahan_rank_reveal_binary128.cpp Kahan_rank_reveal_binary80.cpp Kahan_rank_reveal_double.cpp Kahan_rank_reveal_dd.cpp Kahan_rank_reveal_qd.cpp Pascal_inverse_mpfr.cpp Pascal_inverse_gmp.cpp Pascal_inverse_binary128.cpp Pascal_inverse_binary80.cpp Pascal_inverse_double.cpp Pascal_inverse_dd.cpp Pascal_inverse_qd.cpp Vandermonde_solve_mpfr.cpp Vandermonde_solve_gmp.cpp Vandermonde_solve_binary128.cpp Vandermonde_solve_binary80.cpp Vandermonde_solve_double.cpp Vandermonde_solve_dd.cpp Vandermonde_solve_qd.cpp Wilkinson_roots_mpfr.cpp Wilkinson_roots_gmp.cpp Wilkinson_roots_binary128.cpp Wilkinson_roots_binary80.cpp Wilkinson_roots_double.cpp Wilkinson_roots_dd.cpp Wilkinson_roots_qd.cpp \
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
