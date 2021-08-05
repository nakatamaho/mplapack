bin_PROGRAMS = xeigtstR_mpfr \
xeigtstC_mpfr

mplapack_mpfr_eigtestdir = $(prefix)/share/mplapack/test/mplapack/eig/mpfr

$(bin_PROGRAMS): ../libeig_mpfr.la ../libeig_override_mpfr.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig_mpfr -leig_override_mpfr -L$(top_builddir)/mplapack/test/matgen -lmatgen_mpfr -L$(top_builddir)/mplapack/reference -lmplapack_mpfr -L$(top_builddir)/mpblas/reference -lmpblas_mpfr 
mplibs = -L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR_mpfr_SOURCES  = ../common/Rchkee.cpp
xeigtstR_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_mpfr_CXXFLAGS = $(cxxflags)

xeigtstC_mpfr_SOURCES  = ../common/Cchkee.cpp
xeigtstC_mpfr_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_mpfr_CXXFLAGS = $(cxxflags)

test: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: nep.in xeigtstR_mpfr
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR_mpfr < nep.in > $@ 2>&1

Rsep.out: sep.in xeigtstR_mpfr
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_mpfr < sep.in > $@ 2>&1

Rse2.out: se2.in xeigtstR_mpfr
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_mpfr < se2.in > $@ 2>&1

Rsvd.out: svd.in xeigtstR_mpfr
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR_mpfr < svd.in > $@ 2>&1

Rec.out: Rec.in xeigtstR_mpfr
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR_mpfr < Rec.in > $@ 2>&1

Red.out: Red.in xeigtstR_mpfr
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR_mpfr < Red.in > $@ 2>&1

Rgg.out: Rgg.in xeigtstR_mpfr
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR_mpfr < Rgg.in > $@ 2>&1

Rgd.out: Rgd.in xeigtstR_mpfr
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR_mpfr < Rgd.in > $@ 2>&1

Rsb.out: Rsb.in xeigtstR_mpfr
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_mpfr < Rsb.in > $@ 2>&1

Rsg.out: Rsg.in xeigtstR_mpfr
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR_mpfr < Rsg.in > $@ 2>&1

Rbal.out: Rbal.in xeigtstR_mpfr
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR_mpfr < Rbal.in > $@ 2>&1

Rbak.out: Rbak.in xeigtstR_mpfr
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR_mpfr < Rbak.in > $@ 2>&1

Rgbal.out: Rgbal.in xeigtstR_mpfr
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR_mpfr < Rgbal.in > $@ 2>&1

Rgbak.out: Rgbak.in xeigtstR_mpfr
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR_mpfr < Rgbak.in > $@ 2>&1

Rbb.out: Rbb.in xeigtstR_mpfr
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR_mpfr < Rbb.in > $@ 2>&1

Rglm.out: glm.in xeigtstR_mpfr
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR_mpfr < glm.in > $@ 2>&1

Rgqr.out: gqr.in xeigtstR_mpfr
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR_mpfr < gqr.in > $@ 2>&1

Rgsv.out: gsv.in xeigtstR_mpfr
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR_mpfr < gsv.in > $@ 2>&1

Rcsd.out: csd.in xeigtstR_mpfr
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR_mpfr < csd.in > $@ 2>&1

Rlse.out: lse.in xeigtstR_mpfr
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR_mpfr < lse.in > $@ 2>&1
