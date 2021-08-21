mplapackexamples_PROGRAMS =

if ENABLE_MPFR
mplapackexamples_PROGRAMS += eigenvalue_mpfr inv_mat_mpfr \
                        hilbert_mpfr condition_num_mpfr
mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/mpfrc++ -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)
mpfr_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp

eigenvalue_mpfr_SOURCES= eigenvalue_mpfr.cpp
eigenvalue_mpfr_CXXFLAGS= $(mpfrcxxflags)
eigenvalue_mpfr_LDFLAGS= $(mpfr_libdepends)

inv_mat_mpfr_SOURCES= inv_mat_mpfr.cpp
inv_mat_mpfr_CXXFLAGS= $(mpfrcxxflags)
inv_mat_mpfr_LDFLAGS= $(mpfr_libdepends)

hilbert_mpfr_SOURCES= hilbert_mpfr.cpp
hilbert_mpfr_CXXFLAGS= $(mpfrcxxflags)
hilbert_mpfr_LDFLAGS= $(mpfr_libdepends)

condition_num_mpfr_SOURCES= condition_num_mpfr.cpp
condition_num_mpfr_CXXFLAGS= $(mpfrcxxflags)
condition_num_mpfr_LDFLAGS= $(mpfr_libdepends)
endif

if ENABLE_GMP
mplapackexamples_PROGRAMS += eigenvalue_gmp inv_mat_gmp \
                        hilbert_gmp condition_num_gmp
mpfrcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

eigenvalue_gmp_SOURCES= eigenvalue_gmp.cpp
eigenvalue_gmp_CXXFLAGS= $(mpfrcxxflags)
eigenvalue_gmp_LDFLAGS= $(gmp_libdepends)

inv_mat_gmp_SOURCES= inv_mat_gmp.cpp
inv_mat_gmp_CXXFLAGS= $(mpfrcxxflags)
inv_mat_gmp_LDFLAGS= $(gmp_libdepends)

hilbert_gmp_SOURCES= hilbert_gmp.cpp
hilbert_gmp_CXXFLAGS= $(mpfrcxxflags)
hilbert_gmp_LDFLAGS= $(gmp_libdepends)

condition_num_gmp_SOURCES= condition_num_gmp.cpp
condition_num_gmp_CXXFLAGS= $(mpfrcxxflags)
condition_num_gmp_LDFLAGS= $(gmp_libdepends)
endif

if ENABLE__FLOAT128
mplapackexamples_PROGRAMS += eigenvalue__Float128 inv_mat__Float128 \
                        hilbert__Float128 condition_num__Float128
_Float128cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
if WANT_QUADMATH
_Float128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 -lquadmath
else
_Float128_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128
endif

eigenvalue__Float128_SOURCES= eigenvalue__Float128.cpp
eigenvalue__Float128_CXXFLAGS= $(_Float128cxxflags)
eigenvalue__Float128_LDFLAGS= $(_Float128_libdepends)

inv_mat__Float128_SOURCES= inv_mat__Float128.cpp
inv_mat__Float128_CXXFLAGS= $(_Float128cxxflags)
inv_mat__Float128_LDFLAGS= $(_Float128_libdepends)

hilbert__Float128_SOURCES= hilbert__Float128.cpp
hilbert__Float128_CXXFLAGS= $(_Float128cxxflags)
hilbert__Float128_LDFLAGS= $(_Float128_libdepends)

condition_num__Float128_SOURCES= condition_num__Float128.cpp
condition_num__Float128_CXXFLAGS= $(_Float128cxxflags)
condition_num__Float128_LDFLAGS= $(_Float128_libdepends)

endif

if ENABLE__FLOAT64X
mplapackexamples_PROGRAMS += eigenvalue__Float64x inv_mat__Float64x \
                        hilbert__Float64x condition_num__Float64x
_Float64xcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
_Float64x_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack__Float64x -L$(top_builddir)/mpblas/reference -lmpblas__Float64x

eigenvalue__Float64x_SOURCES= eigenvalue__Float64x.cpp
eigenvalue__Float64x_CXXFLAGS= $(_Float64xcxxflags)
eigenvalue__Float64x_LDFLAGS= $(_Float64x_libdepends)

inv_mat__Float64x_SOURCES= inv_mat__Float64x.cpp
inv_mat__Float64x_CXXFLAGS= $(_Float64xcxxflags)
inv_mat__Float64x_LDFLAGS= $(_Float64x_libdepends)

hilbert__Float64x_SOURCES= hilbert__Float64x.cpp
hilbert__Float64x_CXXFLAGS= $(_Float64xcxxflags)
hilbert__Float64x_LDFLAGS= $(_Float64x_libdepends)

condition_num__Float64x_SOURCES= condition_num__Float64x.cpp
condition_num__Float64x_CXXFLAGS= $(_Float64xcxxflags)
condition_num__Float64x_LDFLAGS= $(_Float64x_libdepends)

endif

mplapackexamplesdir=$(prefix)/share/mplapack/examples/99_Miscellaneous
