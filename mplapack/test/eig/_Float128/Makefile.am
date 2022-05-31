check_PROGRAMS = $(mplapack__Float128_test_PROGRAMS)

mplapack__Float128_test_PROGRAMS = xeigtstR__Float128 xeigtstC__Float128
mplapack__Float128_testdir = $(prefix)/lib/$(target)/mplapack/eig

$(check_PROGRAMS): ../libeig__Float128.la ../libeig_override__Float128.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig__Float128 -leig_override__Float128 -L$(top_builddir)/mplapack/test/matgen -lmatgen__Float128 -L$(top_builddir)/mplapack/reference -lmplapack__Float128 -L$(top_builddir)/mpblas/reference -lmpblas__Float128 
if WANT_QUADMATH
mplibs = -lm -lquadmath
else
mplibs = -lm
endif
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH__FLOAT128___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR__Float128_SOURCES  = ../common/Rchkee.cpp
xeigtstR__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR__Float128_CXXFLAGS = $(cxxflags)

xeigtstC__Float128_SOURCES  = ../common/Cchkee.cpp
xeigtstC__Float128_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC__Float128_CXXFLAGS = $(cxxflags)

check: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out \
       Cnep.out Csep.out Cse2.out Csvd.out Cec.out Ced.out Cgg.out Cgd.out Csb.out Csg.out Cbal.out Cbak.out Cgbal.out Cgbak.out Cbb.out Cglm.out Cgqr.out Cgsv.out Ccsd.out Clse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: ../nep.in xeigtstR__Float128
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < ../nep.in > $@ 2>&1

Rsep.out: ../sep.in xeigtstR__Float128
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < ../sep.in > $@ 2>&1

Rse2.out: ../se2.in xeigtstR__Float128
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < ../se2.in > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR__Float128
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR__Float128 < ../svd.in > $@ 2>&1

Rec.out: ../Rec.in xeigtstR__Float128
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR__Float128 < ../Rec.in > $@ 2>&1

Red.out: ../Red.in xeigtstR__Float128
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR__Float128 < ../Red.in > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR__Float128
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR__Float128 < ../Rgg.in > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR__Float128
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR__Float128 < ../Rgd.in > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR__Float128
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR__Float128 < ../Rsb.in > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR__Float128
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR__Float128 < ../Rsg.in > $@ 2>&1

Rbal.out: ../Rbal.in xeigtstR__Float128
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR__Float128 < ../Rbal.in > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR__Float128
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR__Float128 < ../Rbak.in > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR__Float128
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR__Float128 < ../Rgbal.in > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR__Float128
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR__Float128 < ../Rgbak.in > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR__Float128
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR__Float128 < ../Rbb.in > $@ 2>&1

Rglm.out: ../glm.in xeigtstR__Float128
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR__Float128 < ../glm.in > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR__Float128
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR__Float128 < ../gqr.in > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR__Float128
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR__Float128 < ../gsv.in > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR__Float128
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR__Float128 < ../csd.in > $@ 2>&1

Rlse.out: ../lse.in xeigtstR__Float128
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR__Float128 < ../lse.in > $@ 2>&1

Cnep.out: ../nep.in ./xeigtstC__Float128
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstC__Float128 < ../nep.in > $@ 2>&1

Csep.out: ../sep.in ./xeigtstC__Float128
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC__Float128 < ../sep.in > $@ 2>&1

Cse2.out: ../se2.in ./xeigtstC__Float128
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC__Float128 < ../se2.in > $@ 2>&1

Csvd.out: ../svd.in ./xeigtstC__Float128
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstC__Float128 < ../svd.in > $@ 2>&1

Cec.out: ../Cec.in ./xeigtstC__Float128
	@echo CEC: Testing COMPLEX Eigen Condition Routines
	./xeigtstC__Float128 < ../Cec.in > $@ 2>&1

Ced.out: ../Ced.in ./xeigtstC__Float128
	@echo CES: Testing COMPLEX Nonsymmetric Schur Form Driver
	./xeigtstC__Float128 < ../Ced.in > $@ 2>&1

Cgg.out: ../Cgg.in ./xeigtstC__Float128
	@echo CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstC__Float128 < ../Cgg.in > $@ 2>&1

Cgd.out: ../Cgd.in ./xeigtstC__Float128
	@echo CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstC__Float128 < ../Cgd.in > $@ 2>&1

Csb.out: ../Csb.in ./xeigtstC__Float128
	@echo CHB: Testing Hermitian Eigenvalue Problem routines
	./xeigtstC__Float128 < ../Csb.in > $@ 2>&1

Csg.out: ../Csg.in ./xeigtstC__Float128
	@echo CSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstC__Float128 < ../Csg.in > $@ 2>&1

Cbal.out: ../Cbal.in ./xeigtstC__Float128
	@echo CGEBAL: Testing the balancing of a COMPLEX general matrix
	./xeigtstC__Float128 < ../Cbal.in > $@ 2>&1

Cbak.out: ../Cbak.in ./xeigtstC__Float128
	@echo CGEBAK: Testing the back transformation of a COMPLEX balanced matrix
	./xeigtstC__Float128 < ../Cbak.in > $@ 2>&1

Cgbal.out: ../Cgbal.in ./xeigtstC__Float128
	@echo CGGBAL: Testing the balancing of a pair of COMPLEX general matrices
	./xeigtstC__Float128 < ../Cgbal.in > $@ 2>&1

Cgbak.out: ../Cgbak.in ./xeigtstC__Float128
	@echo CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices
	./xeigtstC__Float128 < ../Cgbak.in > $@ 2>&1

Cbb.out: ../Cbb.in ./xeigtstC__Float128
	@echo CBB: Testing banded Singular Value Decomposition routines
	./xeigtstC__Float128 < ../Cbb.in > $@ 2>&1

Cglm.out: ../glm.in ./xeigtstC__Float128
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstC__Float128 < ../glm.in > $@ 2>&1

Cgqr.out: ../gqr.in ./xeigtstC__Float128
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstC__Float128 < ../gqr.in > $@ 2>&1

Cgsv.out: ../gsv.in ./xeigtstC__Float128
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstC__Float128 < ../gsv.in > $@ 2>&1

Ccsd.out: ../csd.in ./xeigtstC__Float128
	@echo CSD: Testing CS Decomposition routines
	./xeigtstC__Float128 < ../csd.in > $@ 2>&1

Clse.out: ../lse.in ./xeigtstC__Float128
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstC__Float128 < ../lse.in > $@ 2>&1
