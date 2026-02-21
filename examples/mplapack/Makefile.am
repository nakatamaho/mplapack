mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += 

mpfr_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += 

gmp_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

endif

if ENABLE_BINARY128
mplapackexamples_PROGRAMS += 

binary128_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if MPLAPACK_BINARY128_MODE_QUADMATH
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128 -lquadmath
else
binary128_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary128 -L$(top_builddir)/mpblas/reference -lmpblas_binary128
endif

endif

if ENABLE_BINARY80
mplapackexamples_PROGRAMS += 

binary80_cxxflags = $(OPENMP_CXXFLAGS)
binary80_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_binary80 -L$(top_builddir)/mpblas/reference -lmpblas_binary80

endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += 

double_cxxflags = $(OPENMP_CXXFLAGS)
double_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

endif

if ENABLE_DD
mplapackexamples_PROGRAMS += 

dd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

endif

if ENABLE_QD
mplapackexamples_PROGRAMS += 

qd_cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -Wl,-rpath,$(libdir) -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

endif

mplapackexamplesdir=$(prefix)/share/examples/mplapack/mplapack

mplapackexamples_DATA =  \
 Makefile.freebsd Makefile.linux Makefile.macos Makefile.mingw

install-data-hook:
if IS_MACOS
	bash $(top_builddir)/misc/fix_dylib_macOS.sh $(mplapackexamplesdir) $(prefix)
endif
