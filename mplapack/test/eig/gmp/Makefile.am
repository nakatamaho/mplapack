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

check: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out

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
