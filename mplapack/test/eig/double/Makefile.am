check_PROGRAMS = $(mplapack_double_test_PROGRAMS)

mplapack_double_test_PROGRAMS = xeigtstR_double xeigtstC_double
mplapack_double_testdir = $(prefix)/lib/$(target)/mplapack/eig

$(check_PROGRAMS): ../libeig_double.la ../libeig_override_double.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig_double -leig_override_double -L$(top_builddir)/mplapack/test/matgen -lmatgen_double -L$(top_builddir)/mplapack/reference -lmplapack_double -L$(top_builddir)/mpblas/reference -lmpblas_double 
mplibs = -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DOUBLE___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR_double_SOURCES  = ../common/Rchkee.cpp
xeigtstR_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_double_CXXFLAGS = $(cxxflags)

xeigtstC_double_SOURCES  = ../common/Cchkee.cpp
xeigtstC_double_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_double_CXXFLAGS = $(cxxflags)

check: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out \
       Cnep.out Csep.out Cse2.out Csvd.out Cec.out Ced.out Cgg.out Cgd.out Csb.out Csg.out Cbal.out Cbak.out Cgbal.out Cgbak.out Cbb.out Cglm.out Cgqr.out Cgsv.out Ccsd.out Clse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: ../nep.in xeigtstR_double
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR_double < ../nep.in > $@ 2>&1

Rsep.out: ../sep.in xeigtstR_double
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_double < ../sep.in > $@ 2>&1

Rse2.out: ../se2.in xeigtstR_double
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_double < ../se2.in > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR_double
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR_double < ../svd.in > $@ 2>&1

Rec.out: ../Rec.in xeigtstR_double
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR_double < ../Rec.in > $@ 2>&1

Red.out: ../Red.in xeigtstR_double
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR_double < ../Red.in > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR_double
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR_double < ../Rgg.in > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR_double
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR_double < ../Rgd.in > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR_double
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_double < ../Rsb.in > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR_double
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR_double < ../Rsg.in > $@ 2>&1

Rbal.out: ../Rbal_double.in xeigtstR_double
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR_double < ../Rbal_double.in > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR_double
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR_double < ../Rbak.in > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR_double
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR_double < ../Rgbal.in > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR_double
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR_double < ../Rgbak.in > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR_double
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR_double < ../Rbb.in > $@ 2>&1

Rglm.out: ../glm.in xeigtstR_double
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR_double < ../glm.in > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR_double
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR_double < ../gqr.in > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR_double
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR_double < ../gsv.in > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR_double
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR_double < ../csd.in > $@ 2>&1

Rlse.out: ../lse.in xeigtstR_double
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR_double < ../lse.in > $@ 2>&1

Cnep.out: ../nep.in ./xeigtstC_double
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstC_double < ../nep.in > $@ 2>&1

Csep.out: ../sep.in ./xeigtstC_double
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC_double < ../sep.in > $@ 2>&1

Cse2.out: ../se2.in ./xeigtstC_double
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC_double < ../se2.in > $@ 2>&1

Csvd.out: ../svd.in ./xeigtstC_double
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstC_double < ../svd.in > $@ 2>&1

Cec.out: ../Cec.in ./xeigtstC_double
	@echo CEC: Testing COMPLEX Eigen Condition Routines
	./xeigtstC_double < ../Cec.in > $@ 2>&1

Ced.out: ../Ced.in ./xeigtstC_double
	@echo CES: Testing COMPLEX Nonsymmetric Schur Form Driver
	./xeigtstC_double < ../Ced.in > $@ 2>&1

Cgg.out: ../Cgg.in ./xeigtstC_double
	@echo CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstC_double < ../Cgg.in > $@ 2>&1

Cgd.out: ../Cgd.in ./xeigtstC_double
	@echo CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstC_double < ../Cgd.in > $@ 2>&1

Csb.out: ../Csb.in ./xeigtstC_double
	@echo CHB: Testing Hermitian Eigenvalue Problem routines
	./xeigtstC_double < ../Csb.in > $@ 2>&1

Csg.out: ../Csg.in ./xeigtstC_double
	@echo CSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstC_double < ../Csg.in > $@ 2>&1

Cbal.out: ../Cbal_double.in ./xeigtstC_double
	@echo CGEBAL: Testing the balancing of a COMPLEX general matrix
	./xeigtstC_double < ../Cbal_double.in > $@ 2>&1

Cbak.out: ../Cbak.in ./xeigtstC_double
	@echo CGEBAK: Testing the back transformation of a COMPLEX balanced matrix
	./xeigtstC_double < ../Cbak.in > $@ 2>&1

Cgbal.out: ../Cgbal.in ./xeigtstC_double
	@echo CGGBAL: Testing the balancing of a pair of COMPLEX general matrices
	./xeigtstC_double < ../Cgbal.in > $@ 2>&1

Cgbak.out: ../Cgbak.in ./xeigtstC_double
	@echo CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices
	./xeigtstC_double < ../Cgbak.in > $@ 2>&1

Cbb.out: ../Cbb.in ./xeigtstC_double
	@echo CBB: Testing banded Singular Value Decomposition routines
	./xeigtstC_double < ../Cbb.in > $@ 2>&1

Cglm.out: ../glm.in ./xeigtstC_double
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstC_double < ../glm.in > $@ 2>&1

Cgqr.out: ../gqr.in ./xeigtstC_double
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstC_double < ../gqr.in > $@ 2>&1

Cgsv.out: ../gsv.in ./xeigtstC_double
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstC_double < ../gsv.in > $@ 2>&1

Ccsd.out: ../csd.in ./xeigtstC_double
	@echo CSD: Testing CS Decomposition routines
	./xeigtstC_double < ../csd.in > $@ 2>&1

Clse.out: ../lse.in ./xeigtstC_double
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstC_double < ../lse.in > $@ 2>&1
