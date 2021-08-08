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

mplapackexamplesdir=$(prefix)/share/mplapack/examples/99_Miscellaneous
