lib_LTLIBRARIES = 

if ENABLE_GMP
lib_LTLIBRARIES += libblaswrapper_gmp.la
endif

if ENABLE_MPFR
lib_LTLIBRARIES += libblaswrapper_mpfr.la
endif

if ENABLE_QD
lib_LTLIBRARIES += libblaswrapper_qd.la
endif

if ENABLE_DD
lib_LTLIBRARIES +=  libblaswrapper_dd.la
endif

if ENABLE_DOUBLE
lib_LTLIBRARIES += libblaswrapper_double.la
endif

if ENABLE___FLOAT128
lib_LTLIBRARIES += libblaswrapper___float128.la
endif

BLASWRAPPER_SOURCES = zaxpy.cpp ../reference/mplapackinit.cpp
#zcopy.cpp zdotc.cpp zdotu.cpp zgbmv.cpp \
#zgemm.cpp zgemv.cpp zgerc.cpp zgeru.cpp zhbmv.cpp zhemm.cpp zhemv.cpp \
#zher.cpp zher2.cpp zher2k.cpp zherk.cpp zhpmv.cpp zhpr.cpp zhpr2.cpp \
#zrotg.cpp zscal.cpp zdrot.cpp zdscal.cpp zswap.cpp zsymm.cpp zsyr2k.cpp \
#zsyrk.cpp ztbmv.cpp ztbsv.cpp ztpmv.cpp ztpsv.cpp ztrmm.cpp ztrmv.cpp \
#ztrsm.cpp ztrsv.cpp dzabs1.cpp dzasum.cpp dznrm2.cpp dasum.cpp daxpy.cpp\
#dcopy.cpp ddot.cpp dgbmv.cpp dgemm.cpp dgemv.cpp dger.cpp dnrm2.cpp drot.cpp\
#drotg.cpp drotm.cpp drotmg.cpp dsbmv.cpp dscal.cpp dspmv.cpp \
#dspr.cpp dspr2.cpp dswap.cpp dsymm.cpp dsymv.cpp dsyr.cpp dsyr2.cpp \
#dsyr2k.cpp dsyrk.cpp dtbmv.cpp dtbsv.cpp dtpmv.cpp dtpsv.cpp dtrmm.cpp \
#dtrmv.cpp dtrsm.cpp dtrsv.cpp izamax.cpp idamax.cpp \
#dset.cpp lsame.cpp xerbla.cpp

if ENABLE_GMP
libblaswrapper_gmp_la_CPPFLAGS = -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
libblaswrapper_gmp_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif

if ENABLE_MPFR
libblaswrapper_mpfr_la_CPPFLAGS = -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_MPFR___ -I$(top_srcdir)/mpfrc++
libblaswrapper_mpfr_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif

if ENABLE_QD
libblaswrapper_qd_la_CPPFLAGS = -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___
libblaswrapper_qd_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif

if ENABLE_DD
libblaswrapper_dd_la_CPPFLAGS = -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___
libblaswrapper_dd_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif

if ENABLE_DOUBLE
libblaswrapper_double_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___
libblaswrapper_double_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif

if ENABLE___FLOAT128
libblaswrapper___float128_la_CPPFLAGS = -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH___FLOAT128___
libblaswrapper___float128_la_SOURCES = $(BLASWRAPPER_SOURCES)
endif
