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
gmpcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR)
gmp_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp -L$(GMP_LIBDIR) -lgmp

eigenvalue_gmp_SOURCES= eigenvalue_gmp.cpp
eigenvalue_gmp_CXXFLAGS= $(gmpcxxflags)
eigenvalue_gmp_LDFLAGS= $(gmp_libdepends)

inv_mat_gmp_SOURCES= inv_mat_gmp.cpp
inv_mat_gmp_CXXFLAGS= $(gmpcxxflags)
inv_mat_gmp_LDFLAGS= $(gmp_libdepends)

hilbert_gmp_SOURCES= hilbert_gmp.cpp
hilbert_gmp_CXXFLAGS= $(gmpcxxflags)
hilbert_gmp_LDFLAGS= $(gmp_libdepends)

condition_num_gmp_SOURCES= condition_num_gmp.cpp
condition_num_gmp_CXXFLAGS= $(gmpcxxflags)
condition_num_gmp_LDFLAGS= $(gmp_libdepends)
endif

if ENABLE_DD
mplapackexamples_PROGRAMS += eigenvalue_dd inv_mat_dd \
                        hilbert_dd condition_num_dd
ddcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
dd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd -L$(QD_LIBDIR) -lqd

eigenvalue_dd_SOURCES= eigenvalue_dd.cpp
eigenvalue_dd_CXXFLAGS= $(ddcxxflags)
eigenvalue_dd_LDFLAGS= $(dd_libdepends)

inv_mat_dd_SOURCES= inv_mat_dd.cpp
inv_mat_dd_CXXFLAGS= $(ddcxxflags)
inv_mat_dd_LDFLAGS= $(dd_libdepends)

hilbert_dd_SOURCES= hilbert_dd.cpp
hilbert_dd_CXXFLAGS= $(ddcxxflags)
hilbert_dd_LDFLAGS= $(dd_libdepends)

condition_num_dd_SOURCES= condition_num_dd.cpp
condition_num_dd_CXXFLAGS= $(ddcxxflags)
condition_num_dd_LDFLAGS= $(dd_libdepends)
endif

if ENABLE_QD
mplapackexamples_PROGRAMS += eigenvalue_qd inv_mat_qd \
                        hilbert_qd condition_num_qd
qdcxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR)
qd_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd -L$(QD_LIBDIR) -lqd

eigenvalue_qd_SOURCES= eigenvalue_qd.cpp
eigenvalue_qd_CXXFLAGS= $(qdcxxflags)
eigenvalue_qd_LDFLAGS= $(qd_libdepends)

inv_mat_qd_SOURCES= inv_mat_qd.cpp
inv_mat_qd_CXXFLAGS= $(qdcxxflags)
inv_mat_qd_LDFLAGS= $(qd_libdepends)

hilbert_qd_SOURCES= hilbert_qd.cpp
hilbert_qd_CXXFLAGS= $(qdcxxflags)
hilbert_qd_LDFLAGS= $(qd_libdepends)

condition_num_qd_SOURCES= condition_num_qd.cpp
condition_num_qd_CXXFLAGS= $(qdcxxflags)
condition_num_qd_LDFLAGS= $(qd_libdepends)
endif

if ENABLE_DOUBLE
mplapackexamples_PROGRAMS += eigenvalue_double inv_mat_double \
                        hilbert_double condition_num_double
doublecxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include
double_libdepends = -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double

eigenvalue_double_SOURCES= eigenvalue_double.cpp
eigenvalue_double_CXXFLAGS= $(doublecxxflags)
eigenvalue_double_LDFLAGS= $(double_libdepends)

inv_mat_double_SOURCES= inv_mat_double.cpp
inv_mat_double_CXXFLAGS= $(doublecxxflags)
inv_mat_double_LDFLAGS= $(double_libdepends)

hilbert_double_SOURCES= hilbert_double.cpp
hilbert_double_CXXFLAGS= $(doublecxxflags)
hilbert_double_LDFLAGS= $(double_libdepends)

condition_num_double_SOURCES= condition_num_double.cpp
condition_num_double_CXXFLAGS= $(doublecxxflags)
condition_num_double_LDFLAGS= $(double_libdepends)

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
