lib_LTLIBRARIES = libmpblas__Float128.la

libmpblas__Float128_la_CPPFLAGS = -g -I. -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) $(OPENMP_CXXFLAGS) -D___MPLAPACK_BUILD_WITH__FLOAT128___ -finline-functions 
libmpblas__Float128_la_SOURCES = \
Rcopy.cpp Rcopy_ref.cpp openmp/Rcopy_omp.cpp \
Rdot.cpp  Rdot_ref.cpp  openmp/Rdot_omp.cpp \
Rgemm.cpp Rgemm_ref.cpp openmp/Rgemm_NN_omp.cpp openmp/Rgemm_NT_omp.cpp openmp/Rgemm_TN_omp.cpp openmp/Rgemm_TT_omp.cpp \
Raxpy.cpp Raxpy_ref.cpp openmp/Raxpy_omp.cpp \
../../reference/Rgemv.cpp \
../../reference/Caxpy.cpp \
../../reference/Ccopy.cpp \
../../reference/Cdotc.cpp \
../../reference/Cdotu.cpp \
../../reference/Cgbmv.cpp \
../../reference/Cgemm.cpp \
../../reference/Cgemv.cpp \
../../reference/Cgerc.cpp \
../../reference/Cgeru.cpp \
../../reference/Chbmv.cpp \
../../reference/Chemm.cpp \
../../reference/Chemv.cpp \
../../reference/Cher.cpp \
../../reference/Cher2.cpp \
../../reference/Cher2k.cpp \
../../reference/Cherk.cpp \
../../reference/Chpmv.cpp \
../../reference/Chpr.cpp \
../../reference/Chpr2.cpp \
../../reference/Crotg.cpp \
../../reference/Cscal.cpp \
../../reference/CRrot.cpp \
../../reference/CRscal.cpp \
../../reference/Cswap.cpp \
../../reference/Csymm.cpp \
../../reference/Csyr2k.cpp \
../../reference/Csyrk.cpp \
../../reference/Ctbmv.cpp \
../../reference/Ctbsv.cpp \
../../reference/Ctpmv.cpp \
../../reference/Ctpsv.cpp \
../../reference/Ctrmm.cpp \
../../reference/Ctrmv.cpp \
../../reference/Ctrsm.cpp \
../../reference/Ctrsv.cpp \
../../reference/RCabs1.cpp \
../../reference/RCasum.cpp \
../../reference/RCnrm2.cpp \
../../reference/Rasum.cpp \
../../reference/Rgbmv.cpp \
../../reference/Rger.cpp \
../../reference/Rnrm2.cpp \
../../reference/Rrot.cpp \
../../reference/Rrotg.cpp \
../../reference/Rrotm.cpp \
../../reference/Rrotmg.cpp \
../../reference/Rsbmv.cpp \
../../reference/Rscal.cpp \
../../reference/Rspmv.cpp \
../../reference/Rspr.cpp \
../../reference/Rspr2.cpp \
../../reference/Rswap.cpp \
../../reference/Rsymm.cpp \
../../reference/Rsymv.cpp \
../../reference/Rsyr.cpp \
../../reference/Rsyr2.cpp \
../../reference/Rsyr2k.cpp \
../../reference/Rsyrk.cpp \
../../reference/Rtbmv.cpp \
../../reference/Rtbsv.cpp \
../../reference/Rtpmv.cpp \
../../reference/Rtpsv.cpp \
../../reference/Rtrmm.cpp \
../../reference/Rtrmv.cpp \
../../reference/Rtrsm.cpp \
../../reference/Rtrsv.cpp \
../../reference/Mlsame.cpp \
../../reference/Mxerbla.cpp \
../../reference/iCamax.cpp  \
../../reference/iRamax.cpp \
../../reference/mplapackinit.cpp

if IS_MACOSX
libmpblas__Float128_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
