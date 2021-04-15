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
Rlatrz.cpp \
Clatrz.cpp \
Rgeqr2p.cpp \
Rgerq2.cpp \
Rgeql2.cpp \
Rgehd2.cpp \
Rgesc2.cpp \
Rorgr2.cpp \
Rormr2.cpp \
Rsycon.cpp \
Clauu2.cpp \
Rormhr.cpp \
Cunmhr.cpp \
Rgeqrt2.cpp \
Rpptrf.cpp \
Cpptrf.cpp \
Rlaed1.cpp \
Ctrexc.cpp \
Rlange.cpp \
Ropgtr.cpp \
Clange.cpp \
Cupgtr.cpp \
Rpbtf2.cpp \
Rormtr.cpp \
Rpocon.cpp \
Rgecon.cpp \
Cpbtf2.cpp \
Rgbtrs.cpp \
Claunhr_col_getrfnp2.cpp \
Rgbtf2.cpp \
Rtpqrt2.cpp \
Ctpqrt2.cpp \
Rgelqt3.cpp \
Clacn2.cpp \
Rtplqt2.cpp \
Rlacn2.cpp \
Rlarzb.cpp \
Rlaed7.cpp \
Cunmqr.cpp \
Rsptrd.cpp \
Rormlq.cpp \
Rlahr2.cpp \
Cunmlq.cpp \
Cunmbr.cpp \
Clarzb.cpp \
Clahr2.cpp \
Rlansb.cpp \
Ctplqt2.cpp \
Ctrsen.cpp \
Rlaed3.cpp \
Cungbr.cpp \
Rlansp.cpp \
Rlatdf.cpp \
Clarfgp.cpp \
Clanhp.cpp \
Clanhb.cpp \
Clargv.cpp \
Clamtsqr.cpp \
Rgebal.cpp \
Rsyequb.cpp \
Ctgex2.cpp \
Rhseqr.cpp \
Rlag2.cpp \
Rlagv2.cpp \
Rtrsen.cpp \
Rptrfs.cpp \
Rlaed6.cpp \
Rporfs.cpp \
Rsyrfs.cpp \
Rlaed0.cpp \
Clabrd.cpp \
Rstedc.cpp \
Rlaed2.cpp \
Ctrsyl.cpp \
Rlasy2.cpp \
Rlaexc.cpp \
Rsptrs.cpp \
Csptrs.cpp \
Rsytrf_aa.cpp \
Rpbtrf.cpp \
Cpbtrf.cpp \
Rlaqr2.cpp \
Rlaqr3.cpp \
Rlaln2.cpp \
Ctrevc3.cpp \
Rgbtrf.cpp \
Clahqr.cpp \
Csytf2.cpp \
Rsptrf.cpp \
Csptrf.cpp \
Rlahqr.cpp \
Rsytri2x.cpp \
Csytri2x.cpp \
Rsytrf_aa_2stage.cpp \
Rsytri_3x.cpp \
Csytri_3x.cpp \
Rtgex2.cpp \
Rsbtrd.cpp \
Rlatrs.cpp \
Rlaqr4.cpp \
Rlaqr0.cpp \
Csytf2_rook.cpp \
Rlaed4.cpp \
Clatrs.cpp \
Clasyf_rook.cpp \
Clasyf_rk.cpp \
Csytf2_rk.cpp \
Rlaqr5.cpp \
Rtrsyl.cpp \
Rtgsy2.cpp \
Clahef_rk.cpp \
Clahef_rook.cpp \
Rtrevc3.cpp
