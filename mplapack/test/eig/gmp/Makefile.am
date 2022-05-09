bin_PROGRAMS = xeigtstR_gmp \
xeigtstC_gmp

mplapack_gmp_eigtestdir = $(prefix)/share/mplapack/test/mplapack/eig/gmp

$(bin_PROGRAMS): ../libeig_gmp.la ../libeig_override_gmp.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig_gmp -leig_override_gmp -L$(top_builddir)/mplapack/test/matgen -lmatgen_gmp -L$(top_builddir)/mplapack/reference -lmplapack_gmp -L$(top_builddir)/mpblas/reference -lmpblas_gmp 
mplibs = -L$(GMP_LIBDIR) -lgmpxx -lgmp -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR_gmp_SOURCES  = ../common/Rchkee.cpp
xeigtstR_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_gmp_CXXFLAGS = $(cxxflags)

xeigtstC_gmp_SOURCES  = ../common/Cchkee.cpp
xeigtstC_gmp_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_gmp_CXXFLAGS = $(cxxflags)

check: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out \
       Cnep.out Csep.out Cse2.out Csvd.out Cec.out Ced.out Cgg.out Cgd.out Csb.out Csg.out Cbal.out Cbak.out Cgbal.out Cgbak.out Cbb.out Cglm.out Cgqr.out Cgsv.out Ccsd.out Clse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: ../nep.in xeigtstR_gmp
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR_gmp < ../nep.in > $@ 2>&1

Rsep.out: ../sep.in xeigtstR_gmp
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_gmp < ../sep.in > $@ 2>&1

Rse2.out: ../se2.in xeigtstR_gmp
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_gmp < ../se2.in > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR_gmp
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR_gmp < ../svd.in > $@ 2>&1

Rec.out: ../Rec.in xeigtstR_gmp
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR_gmp < ../Rec.in > $@ 2>&1

Red.out: ../Red.in xeigtstR_gmp
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR_gmp < ../Red.in > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR_gmp
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR_gmp < ../Rgg.in > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR_gmp
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR_gmp < ../Rgd.in > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR_gmp
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_gmp < ../Rsb.in > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR_gmp
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR_gmp < ../Rsg.in > $@ 2>&1

Rbal.out: ../Rbal.in xeigtstR_gmp
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR_gmp < ../Rbal.in > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR_gmp
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR_gmp < ../Rbak.in > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR_gmp
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR_gmp < ../Rgbal.in > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR_gmp
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR_gmp < ../Rgbak.in > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR_gmp
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR_gmp < ../Rbb.in > $@ 2>&1

Rglm.out: ../glm.in xeigtstR_gmp
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR_gmp < ../glm.in > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR_gmp
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR_gmp < ../gqr.in > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR_gmp
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR_gmp < ../gsv.in > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR_gmp
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR_gmp < ../csd.in > $@ 2>&1

Rlse.out: ../lse.in xeigtstR_gmp
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR_gmp < ../lse.in > $@ 2>&1

Cnep.out: ../nep.in ./xeigtstC_gmp
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstC_gmp < ../nep.in > $@ 2>&1

Csep.out: ../sep.in ./xeigtstC_gmp
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC_gmp < ../sep.in > $@ 2>&1

Cse2.out: ../se2.in ./xeigtstC_gmp
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC_gmp < ../se2.in > $@ 2>&1

Csvd.out: ../svd.in ./xeigtstC_gmp
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstC_gmp < ../svd.in > $@ 2>&1

Cec.out: ../Cec.in ./xeigtstC_gmp
	@echo CEC: Testing COMPLEX Eigen Condition Routines
	./xeigtstC_gmp < ../Cec.in > $@ 2>&1

Ced.out: ../Ced.in ./xeigtstC_gmp
	@echo CES: Testing COMPLEX Nonsymmetric Schur Form Driver
	./xeigtstC_gmp < ../Ced.in > $@ 2>&1

Cgg.out: ../Cgg.in ./xeigtstC_gmp
	@echo CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstC_gmp < ../Cgg.in > $@ 2>&1

Cgd.out: ../Cgd.in ./xeigtstC_gmp
	@echo CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstC_gmp < ../Cgd.in > $@ 2>&1

Csb.out: ../Csb.in ./xeigtstC_gmp
	@echo CHB: Testing Hermitian Eigenvalue Problem routines
	./xeigtstC_gmp < ../Csb.in > $@ 2>&1

Csg.out: ../Csg.in ./xeigtstC_gmp
	@echo CSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstC_gmp < ../Csg.in > $@ 2>&1

Cbal.out: ../Cbal.in ./xeigtstC_gmp
	@echo CGEBAL: Testing the balancing of a COMPLEX general matrix
	./xeigtstC_gmp < ../Cbal.in > $@ 2>&1

Cbak.out: ../Cbak.in ./xeigtstC_gmp
	@echo CGEBAK: Testing the back transformation of a COMPLEX balanced matrix
	./xeigtstC_gmp < ../Cbak.in > $@ 2>&1

Cgbal.out: ../Cgbal.in ./xeigtstC_gmp
	@echo CGGBAL: Testing the balancing of a pair of COMPLEX general matrices
	./xeigtstC_gmp < ../Cgbal.in > $@ 2>&1

Cgbak.out: ../Cgbak.in ./xeigtstC_gmp
	@echo CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices
	./xeigtstC_gmp < ../Cgbak.in > $@ 2>&1

Cbb.out: ../Cbb.in ./xeigtstC_gmp
	@echo CBB: Testing banded Singular Value Decomposition routines
	./xeigtstC_gmp < ../Cbb.in > $@ 2>&1

Cglm.out: ../glm.in ./xeigtstC_gmp
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstC_gmp < ../glm.in > $@ 2>&1

Cgqr.out: ../gqr.in ./xeigtstC_gmp
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstC_gmp < ../gqr.in > $@ 2>&1

Cgsv.out: ../gsv.in ./xeigtstC_gmp
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstC_gmp < ../gsv.in > $@ 2>&1

Ccsd.out: ../csd.in ./xeigtstC_gmp
	@echo CSD: Testing CS Decomposition routines
	./xeigtstC_gmp < ../csd.in > $@ 2>&1

Clse.out: ../lse.in ./xeigtstC_gmp
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstC_gmp < ../lse.in > $@ 2>&1
