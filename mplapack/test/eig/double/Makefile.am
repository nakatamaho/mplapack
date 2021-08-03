bin_PROGRAMS = xeigtstR_double \
xeigtstC_double

mplapack_double_eigtestdir = $(prefix)/share/mplapack/test/mplapack/eig/double

$(bin_PROGRAMS): ../libeig_double.la ../libeig_override.la

libdepends = -L$(top_buildoubleir)/mplapack/test/eig -leig_double -leig_override -L$(top_buildoubleir)/mplapack/test/matgen -lmatgen_double -L$(top_buildoubleir)/mplapack/reference -lmplapack_double -L$(top_buildoubleir)/mpblas/reference -lmpblas_double 
mplibs = -L$(MPC_LIBDIR) -L$(DOUBLE_LIBDIR) -L$(GMP_LIBDIR) -lmpc -ldouble -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -I$(DOUBLE_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/doublec++ -D___MPLAPACK_BUILD_WITH_DOUBLE___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR_double_SOURCES  = ../common/Rchkee.cpp
xeigtstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_double_CXXFLAGS = $(cxxflags)

xeigtstC_double_SOURCES  = ../common/Cchkee.cpp
xeigtstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_double_CXXFLAGS = $(cxxflags)

test: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: nep.in xeigtstR_double
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR_double < nep.in > $@ 2>&1

Rsep.out: sep.in xeigtstR_double
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_double < sep.in > $@ 2>&1

Rse2.out: se2.in xeigtstR_double
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_double < se2.in > $@ 2>&1

Rsvd.out: svd.in xeigtstR_double
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR_double < svd.in > $@ 2>&1

Rec.out: Rec.in xeigtstR_double
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR_double < Rec.in > $@ 2>&1

Red.out: Red.in xeigtstR_double
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR_double < Red.in > $@ 2>&1

Rgg.out: Rgg.in xeigtstR_double
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR_double < Rgg.in > $@ 2>&1

Rgd.out: Rgd.in xeigtstR_double
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR_double < Rgd.in > $@ 2>&1

Rsb.out: Rsb.in xeigtstR_double
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_double < Rsb.in > $@ 2>&1

Rsg.out: Rsg.in xeigtstR_double
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR_double < Rsg.in > $@ 2>&1

Rbal.out: Rbal.in xeigtstR_double
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR_double < Rbal.in > $@ 2>&1

Rbak.out: Rbak.in xeigtstR_double
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR_double < Rbak.in > $@ 2>&1

Rgbal.out: Rgbal.in xeigtstR_double
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR_double < Rgbal.in > $@ 2>&1

Rgbak.out: Rgbak.in xeigtstR_double
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR_double < Rgbak.in > $@ 2>&1

Rbb.out: Rbb.in xeigtstR_double
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR_double < Rbb.in > $@ 2>&1

Rglm.out: glm.in xeigtstR_double
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR_double < glm.in > $@ 2>&1

Rgqr.out: gqr.in xeigtstR_double
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR_double < gqr.in > $@ 2>&1

Rgsv.out: gsv.in xeigtstR_double
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR_double < gsv.in > $@ 2>&1

Rcsd.out: csd.in xeigtstR_double
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR_double < csd.in > $@ 2>&1

Rlse.out: lse.in xeigtstR_double
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR_double < lse.in > $@ 2>&1
