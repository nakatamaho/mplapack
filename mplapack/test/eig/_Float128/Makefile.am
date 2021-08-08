bin_PROGRAMS = \
xeigtstR__Float128 \
xeigtstC__Float128

mplapack__Float128_eigtestdir = $(prefix)/share/mplapack/test/mplapack/eig/_Float128

$(bin_PROGRAMS): ../libeig__Float128.la ../libeig_override__Float128.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig__Float128 -leig_override__Float128 -L$(top_builddir)/mplapack/test/matgen -lmatgen__Float128 -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 
mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR__Float128_SOURCES  = ../common/Rchkee.cpp
xeigtstR__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR__Float128_CXXFLAGS = $(cxxflags)

xeigtstC__Float128_SOURCES  = ../common/Cchkee.cpp
xeigtstC__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC__Float128_CXXFLAGS = $(cxxflags)

test: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: nep.in xeigtstR__Float128
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < nep.in > $@ 2>&1

Rsep.out: sep.in xeigtstR__Float128
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < sep.in > $@ 2>&1

Rse2.out: se2.in xeigtstR__Float128
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < se2.in > $@ 2>&1

Rsvd.out: svd.in xeigtstR__Float128
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR__Float128 < svd.in > $@ 2>&1

Rec.out: Rec.in xeigtstR__Float128
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR__Float128 < Rec.in > $@ 2>&1

Red.out: Red.in xeigtstR__Float128
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR__Float128 < Red.in > $@ 2>&1

Rgg.out: Rgg.in xeigtstR__Float128
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR__Float128 < Rgg.in > $@ 2>&1

Rgd.out: Rgd.in xeigtstR__Float128
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR__Float128 < Rgd.in > $@ 2>&1

Rsb.out: Rsb.in xeigtstR__Float128
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < Rsb.in > $@ 2>&1

Rsg.out: Rsg.in xeigtstR__Float128
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR__Float128 < Rsg.in > $@ 2>&1

Rbal.out: Rbal.in xeigtstR__Float128
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR__Float128 < Rbal.in > $@ 2>&1

Rbak.out: Rbak.in xeigtstR__Float128
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR__Float128 < Rbak.in > $@ 2>&1

Rgbal.out: Rgbal.in xeigtstR__Float128
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR__Float128 < Rgbal.in > $@ 2>&1

Rgbak.out: Rgbak.in xeigtstR__Float128
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR__Float128 < Rgbak.in > $@ 2>&1

Rbb.out: Rbb.in xeigtstR__Float128
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR__Float128 < Rbb.in > $@ 2>&1

Rglm.out: glm.in xeigtstR__Float128
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR__Float128 < glm.in > $@ 2>&1

Rgqr.out: gqr.in xeigtstR__Float128
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR__Float128 < gqr.in > $@ 2>&1

Rgsv.out: gsv.in xeigtstR__Float128
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR__Float128 < gsv.in > $@ 2>&1

Rcsd.out: csd.in xeigtstR__Float128
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR__Float128 < csd.in > $@ 2>&1

Rlse.out: lse.in xeigtstR__Float128
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR__Float128 < lse.in > $@ 2>&1
