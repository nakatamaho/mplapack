lib_LTLIBRARIES = 
if ENABLE_MPFR
lib_LTLIBRARIES += libmplapack_mpfr_ref.la
endif

#if ENABLE_GMP
#lib_LTLIBRARIES += libmplapack_gmp_ref.la
#endif

#if ENABLE_QD
#lib_LTLIBRARIES += libmplapack_qd_ref.la
#endif

#if ENABLE_DD
#lib_LTLIBRARIES += libmplapack_dd_ref.la
#endif

#if ENABLE_DOUBLE
#lib_LTLIBRARIES += libmplapack_double_ref.la
#endif

#if ENABLE__FLOAT64X
#lib_LTLIBRARIES += libmplapack__Float64x_ref.la
#endif

#if ENABLE__FLOAT128
#lib_LTLIBRARIES += libmplapack__Float128_ref.la
#endif

if ENABLE_GMP
libmplapack_gmp_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_GMP___
libmplapack_gmp_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_gmp_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_MPFR
libmplapack_mpfr_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR)  -I$(MPC_INCLUDEDIR) -I$(top_srcdir)/mpfrc++ -D___MPLAPACK_BUILD_WITH_MPFR___
libmplapack_mpfr_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_mpfr_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_QD
libmplapack_qd_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_QD___
libmplapack_qd_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_qd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DD
libmplapack_dd_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) -D___MPLAPACK_BUILD_WITH_DD___
libmplapack_dd_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_dd_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE_DOUBLE
libmplapack_double_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH_DOUBLE___
libmplapack_double_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack_double_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT64X
libmplapack__Float64x_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT64X___
libmplapack__Float64x_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack__Float64x_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

if ENABLE__FLOAT128
libmplapack__Float128_ref_la_CPPFLAGS = $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -D___MPLAPACK_BUILD_WITH__FLOAT128___
libmplapack__Float128_ref_la_SOURCES = $(MPLAPACK_SOURCES)
if IS_MACOSX
libmplapack__Float128_ref_la_LDFLAGS = -Wl,-flat_namespace,-undefined,dynamic_lookup
endif
endif

MPLAPACK_SOURCES = \
Rlals0.cpp \
Rlalsa.cpp \
Rlalsd.cpp \
Rgelsd.cpp \
Rlaln2.cpp \
Ctrevc3.cpp \
Rggevx.cpp \
Rsytrf_aa_2stage.cpp \
Rsytri_3x.cpp \
Csytri_3x.cpp \
Rtgex2.cpp \
Rgesvdx.cpp \
Rgsvj1.cpp \
Csytf2_rook.cpp \
Rgghd3.cpp \
Clasyf_rook.cpp \
Clasyf_rk.cpp \
Csytf2_rk.cpp \
Rlaqr5.cpp \
Rbbcsd.cpp \
Rgsvj0.cpp \
Clahef_rk.cpp \
Rgesvdq.cpp \
Clahef_rook.cpp \
Rhgeqz.cpp \
Rtrevc3.cpp \
Rgesvj.cpp \
Rgejsv.cpp \
Rgesdd.cpp
