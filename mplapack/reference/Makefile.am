lib_LTLIBRARIES = 

#if ENABLE_MPFR
#lib_LTLIBRARIES += libmplapack_mpfr_ref.la
#endif

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

if ENABLE__FLOAT128
lib_LTLIBRARIES += libmplapack__Float128_ref.la
endif

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
Rsteqr.cpp \
Rcombssq.cpp \
Rgesv.cpp \
Rgetf2.cpp \
Rgetrf.cpp \
Rgetrf2.cpp \
Rgetri.cpp \
Rgetrs.cpp \
Risnan.cpp \
Rladiv.cpp \
Rlae2.cpp \
Rlaev2.cpp \
Rlaisnan.cpp \
Rlanst.cpp \
Rlansy.cpp \
Rlaorhr_col_getrfnp.cpp \
Rlaorhr_col_getrfnp2.cpp \
Rlapy2.cpp \
Rlapy3.cpp \
Rlarf.cpp \
Rlarfb.cpp \
Rlarfb_gett.cpp \
Rlarfg.cpp \
Rlarfgp.cpp \
Rlarft.cpp \
Rlarfy.cpp \
Rlartg.cpp \
Rlartgp.cpp \
Rlartgs.cpp \
Rlartv.cpp \
Rlascl.cpp \
Rlascl2.cpp \
Rlaset.cpp \
Rlasr.cpp \
Rlasrt.cpp \
Rlassq.cpp \
Rlaswp.cpp \
Rlasyf.cpp \
Rlasyf_aa.cpp \
Rlasyf_rk.cpp \
Rlasyf_rook.cpp \
Rlatrd.cpp \
Rorg2l.cpp \
Rorg2r.cpp \
Rorgql.cpp \
Rorgqr.cpp \
Rorgtr.cpp \
Rpotf2.cpp \
Rpotrf.cpp \
Rpotrf2.cpp \
Rsprfs.cpp \
Rsterf.cpp \
Rsyev.cpp \
Rsyev_2stage.cpp \
Rsyevd.cpp \
Rsyevd_2stage.cpp \
Rsyevr.cpp \
Rsyevr_2stage.cpp \
Rsyevx_2stage.cpp \
Rsyrfs.cpp \
Rsytd2.cpp \
Rsytrd.cpp \
Rsytrd_2stage.cpp \
Rsytrd_sy2sb.cpp \
Rtrti2.cpp \
Rtrtri.cpp \
Rtrtrs.cpp \
iladlc.cpp \
iladlr.cpp \
iCmax1.cpp \
Cgesv.cpp \
Cgetf2.cpp \
Cgetrf.cpp \
Cgetrf2.cpp \
Cgetri.cpp \
Cgetrs.cpp \
Cheev.cpp \
Cheev_2stage.cpp \
Cheevd.cpp \
Cheevd_2stage.cpp \
Cheevr.cpp \
Cheevr_2stage.cpp \
Cheevx_2stage.cpp \
Chetd2.cpp \
Chetrd.cpp \
Chetrd_2stage.cpp \
Chetrd_he2hb.cpp \
Clacgv.cpp \
Clacrm.cpp \
Clacrt.cpp \
Cladiv.cpp \
Claesy.cpp \
Claev2.cpp \
Clahef.cpp \
Clahef_aa.cpp \
Clahef_rk.cpp \
Clahef_rook.cpp \
Clanhe.cpp \
Clansy.cpp \
Clarf.cpp \
Clarfb.cpp \
Clarfb_gett.cpp \
Clarfg.cpp \
Clarfgp.cpp \
Clarft.cpp \
Clarfy.cpp \
Clartg.cpp \
Clartv.cpp \
Clascl.cpp \
Clascl2.cpp \
Claset.cpp \
Clasr.cpp \
Classq.cpp \
Claswp.cpp \
Clasyf.cpp \
Clasyf_aa.cpp \
Clasyf_rk.cpp \
Clasyf_rook.cpp \
Clatrd.cpp \
Claunhr_col_getrfnp.cpp \
Claunhr_col_getrfnp2.cpp \
Cpotf2.cpp \
Cpotrf.cpp \
Cpotrf2.cpp \
Crot.cpp \
Cspmv.cpp \
Cspr.cpp \
Csprfs.cpp \
Csteqr.cpp \
Csymv.cpp \
Csyr.cpp \
Csyrfs.cpp \
Ctrti2.cpp \
Ctrtri.cpp \
Ctrtrs.cpp \
Cung2l.cpp \
Cung2r.cpp \
Cungql.cpp \
Cungqr.cpp \
Cungtr.cpp \
Rlamch.cpp 
