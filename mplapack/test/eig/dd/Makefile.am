bin_PROGRAMS = xeigtstR_dd \
xeigtstC_dd

mplapack_dd_eigtestdir = $(prefix)/share/mplapack/test/mplapack/eig/dd

$(bin_PROGRAMS): ../libeig_dd.la ../libeig_override.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig_dd -leig_override -L$(top_builddir)/mplapack/test/matgen -lmatgen_dd -L$(top_builddir)/mplapack/reference -lmplapack_dd -L$(top_builddir)/mpblas/reference -lmpblas_dd 
mplibs = -L$(MPC_LIBDIR) -L$(DD_LIBDIR) -L$(GMP_LIBDIR) -lmpc -ldd -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(DD_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/ddc++ -D___MPLAPACK_BUILD_WITH_DD___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR_dd_SOURCES  = ../common/Rchkee.cpp  ../../../../mpblas/test/common/mplapack.test.cpp
xeigtstR_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_dd_CXXFLAGS = $(cxxflags)

xeigtstC_dd_SOURCES  = ../common/Cchkee.cpp ../../../../mpblas/test/common/mplapack.test.cpp
xeigtstC_dd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_dd_CXXFLAGS = $(cxxflags)

test: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: nep.in xeigtstR_dd
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR_dd < nep.in > $@ 2>&1

Rsep.out: sep.in xeigtstR_dd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_dd < sep.in > $@ 2>&1

Rse2.out: se2.in xeigtstR_dd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_dd < se2.in > $@ 2>&1

Rsvd.out: svd.in xeigtstR_dd
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR_dd < svd.in > $@ 2>&1

Rec.out: Rec.in xeigtstR_dd
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR_dd < Rec.in > $@ 2>&1

Red.out: Red.in xeigtstR_dd
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR_dd < Red.in > $@ 2>&1

Rgg.out: Rgg.in xeigtstR_dd
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR_dd < Rgg.in > $@ 2>&1

Rgd.out: Rgd.in xeigtstR_dd
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR_dd < Rgd.in > $@ 2>&1

Rsb.out: Rsb.in xeigtstR_dd
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_dd < Rsb.in > $@ 2>&1

Rsg.out: Rsg.in xeigtstR_dd
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR_dd < Rsg.in > $@ 2>&1

Rbal.out: Rbal.in xeigtstR_dd
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR_dd < Rbal.in > $@ 2>&1

Rbak.out: Rbak.in xeigtstR_dd
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR_dd < Rbak.in > $@ 2>&1

Rgbal.out: Rgbal.in xeigtstR_dd
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR_dd < Rgbal.in > $@ 2>&1

Rgbak.out: Rgbak.in xeigtstR_dd
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR_dd < Rgbak.in > $@ 2>&1

Rbb.out: Rbb.in xeigtstR_dd
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR_dd < Rbb.in > $@ 2>&1

Rglm.out: glm.in xeigtstR_dd
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR_dd < glm.in > $@ 2>&1

Rgqr.out: gqr.in xeigtstR_dd
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR_dd < gqr.in > $@ 2>&1

Rgsv.out: gsv.in xeigtstR_dd
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR_dd < gsv.in > $@ 2>&1

Rcsd.out: csd.in xeigtstR_dd
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR_dd < csd.in > $@ 2>&1

Rlse.out: lse.in xeigtstR_dd
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR_dd < lse.in > $@ 2>&1
