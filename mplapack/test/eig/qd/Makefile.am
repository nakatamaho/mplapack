check_PROGRAMS = $(mplapack_qd_test_PROGRAMS)

mplapack_qd_test_PROGRAMS = xeigtstR_qd xeigtstC_qd
mplapack_qd_testdir = $(prefix)/lib/$(target)/mplapack/test/eig

$(check_PROGRAMS): ../libeig_qd.la ../libeig_override_qd.la

libdepends = -L$(top_builddir)/mplapack/test/eig -leig_qd -leig_override_qd -L$(top_builddir)/mplapack/test/matgen -lmatgen_qd -L$(top_builddir)/mplapack/reference -lmplapack_qd -L$(top_builddir)/mpblas/reference -lmpblas_qd 
mplibs = -L$(QD_LIBDIR) -lqd -lm
cxxflags = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(LAPACK_INCLUDEDIR) -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___ -D___MPLAPACK_INTERNAL___ -I$(top_srcdir)/mplapack/test/fable

xeigtstR_qd_SOURCES  = ../common/Rchkee.cpp
xeigtstR_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstR_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstR_qd_CXXFLAGS = $(cxxflags)

xeigtstC_qd_SOURCES  = ../common/Cchkee.cpp
xeigtstC_qd_LDFLAGS  = $(libdepends) $(mplibs) -L$(LAPACK_LIBDIR) -llapacke -llapack -lblas -lgfortran
xeigtstC_qd_LDFLAGS += -Wl,-rpath,$(libdir)
xeigtstC_qd_CXXFLAGS = $(cxxflags)

check: Rnep.out Rsep.out Rse2.out Rsvd.out Rec.out Red.out Rgg.out Rgd.out Rsb.out Rsg.out Rbal.out Rbak.out Rgbal.out Rgbak.out Rbb.out Rglm.out Rgqr.out Rgsv.out Rcsd.out Rlse.out \
       Cnep.out Csep.out Cse2.out Csvd.out Cec.out Ced.out Cgg.out Cgd.out Csb.out Csg.out Cbal.out Cbak.out Cgbal.out Cgbak.out Cbb.out Cglm.out Cgqr.out Cgsv.out Ccsd.out Clse.out

clean-local:
	rm -f *.out *~ \#*

Rnep.out: ../nep.in xeigtstR_qd
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstR_qd < ../nep.in > $@ 2>&1

Rsep.out: ../sep.in xeigtstR_qd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_qd < ../sep.in > $@ 2>&1

Rse2.out: ../se2.in xeigtstR_qd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_qd < ../se2.in > $@ 2>&1

Rsvd.out: ../svd.in xeigtstR_qd
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstR_qd < ../svd.in > $@ 2>&1

Rec.out: ../Rec.in xeigtstR_qd
	@echo DEC: Testing Eigen Condition Routines
	./xeigtstR_qd < ../Rec.in > $@ 2>&1

Red.out: ../Red.in xeigtstR_qd
	@echo DEV: Testing Nonsymmetric Eigenvalue Driver
	./xeigtstR_qd < ../Red.in > $@ 2>&1

Rgg.out: ../Rgg.in xeigtstR_qd
	@echo DGG: Testing Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstR_qd < ../Rgg.in > $@ 2>&1

Rgd.out: ../Rgd.in xeigtstR_qd
	@echo DGD: Testing Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstR_qd < ../Rgd.in > $@ 2>&1

Rsb.out: ../Rsb.in xeigtstR_qd
	@echo DSB: Testing Symmetric Eigenvalue Problem routines
	./xeigtstR_qd < ../Rsb.in > $@ 2>&1

Rsg.out: ../Rsg.in xeigtstR_qd
	@echo DSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstR_qd < ../Rsg.in > $@ 2>&1

Rbal.out: ../Rbal.in xeigtstR_qd
	@echo Rgebal: Testing the balancing of a general matrix
	./xeigtstR_qd < ../Rbal.in > $@ 2>&1

Rbak.out: ../Rbak.in xeigtstR_qd
	@echo Rgebak: Testing the back transformation of balanced matrix
	./xeigtstR_qd < ../Rbak.in > $@ 2>&1

Rgbal.out: ../Rgbal.in xeigtstR_qd
	@echo Rggbal: Testing the balancing of a pair of general matrices
	./xeigtstR_qd < ../Rgbal.in > $@ 2>&1

Rgbak.out: ../Rgbak.in xeigtstR_qd
	@echo Rggbak: Testing the back transformation of a pair of balanced matrices
	./xeigtstR_qd < ../Rgbak.in > $@ 2>&1

Rbb.out: ../Rbb.in xeigtstR_qd
	@echo DBB: Testing banded Singular Value Decomposition routines
	./xeigtstR_qd < ../Rbb.in > $@ 2>&1

Rglm.out: ../glm.in xeigtstR_qd
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstR_qd < ../glm.in > $@ 2>&1

Rgqr.out: ../gqr.in xeigtstR_qd
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstR_qd < ../gqr.in > $@ 2>&1

Rgsv.out: ../gsv.in xeigtstR_qd
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstR_qd < ../gsv.in > $@ 2>&1

Rcsd.out: ../csd.in xeigtstR_qd
	@echo CSD: Testing CS Decomposition routines
	./xeigtstR_qd < ../csd.in > $@ 2>&1

Rlse.out: ../lse.in xeigtstR_qd
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstR_qd < ../lse.in > $@ 2>&1

Cnep.out: ../nep.in ./xeigtstC_qd
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstC_qd < ../nep.in > $@ 2>&1

Csep.out: ../sep.in ./xeigtstC_qd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC_qd < ../sep.in > $@ 2>&1

Cse2.out: ../se2.in ./xeigtstC_qd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstC_qd < ../se2.in > $@ 2>&1

Csvd.out: ../svd.in ./xeigtstC_qd
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstC_qd < ../svd.in > $@ 2>&1

Cec.out: ../Cec.in ./xeigtstC_qd
	@echo CEC: Testing COMPLEX Eigen Condition Routines
	./xeigtstC_qd < ../Cec.in > $@ 2>&1

Ced.out: ../Ced.in ./xeigtstC_qd
	@echo CES: Testing COMPLEX Nonsymmetric Schur Form Driver
	./xeigtstC_qd < ../Ced.in > $@ 2>&1

Cgg.out: ../Cgg.in ./xeigtstC_qd
	@echo CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstC_qd < ../Cgg.in > $@ 2>&1

Cgd.out: ../Cgd.in ./xeigtstC_qd
	@echo CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstC_qd < ../Cgd.in > $@ 2>&1

Csb.out: ../Csb.in ./xeigtstC_qd
	@echo CHB: Testing Hermitian Eigenvalue Problem routines
	./xeigtstC_qd < ../Csb.in > $@ 2>&1

Csg.out: ../Csg.in ./xeigtstC_qd
	@echo CSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstC_qd < ../Csg.in > $@ 2>&1

Cbal.out: ../Cbal.in ./xeigtstC_qd
	@echo CGEBAL: Testing the balancing of a COMPLEX general matrix
	./xeigtstC_qd < ../Cbal.in > $@ 2>&1

Cbak.out: ../Cbak.in ./xeigtstC_qd
	@echo CGEBAK: Testing the back transformation of a COMPLEX balanced matrix
	./xeigtstC_qd < ../Cbak.in > $@ 2>&1

Cgbal.out: ../Cgbal.in ./xeigtstC_qd
	@echo CGGBAL: Testing the balancing of a pair of COMPLEX general matrices
	./xeigtstC_qd < ../Cgbal.in > $@ 2>&1

Cgbak.out: ../Cgbak.in ./xeigtstC_qd
	@echo CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices
	./xeigtstC_qd < ../Cgbak.in > $@ 2>&1

Cbb.out: ../Cbb.in ./xeigtstC_qd
	@echo CBB: Testing banded Singular Value Decomposition routines
	./xeigtstC_qd < ../Cbb.in > $@ 2>&1

Cglm.out: ../glm.in ./xeigtstC_qd
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstC_qd < ../glm.in > $@ 2>&1

Cgqr.out: ../gqr.in ./xeigtstC_qd
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstC_qd < ../gqr.in > $@ 2>&1

Cgsv.out: ../gsv.in ./xeigtstC_qd
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstC_qd < ../gsv.in > $@ 2>&1

Ccsd.out: ../csd.in ./xeigtstC_qd
	@echo CSD: Testing CS Decomposition routines
	./xeigtstC_qd < ../csd.in > $@ 2>&1

Clse.out: ../lse.in ./xeigtstC_qd
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstC_qd < ../lse.in > $@ 2>&1
