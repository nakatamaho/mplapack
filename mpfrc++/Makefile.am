if ENABLE_MPFR

lib_LTLIBRARIES = libmpfrcxx.la

libmpfrcxx_la_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/$(GMP_INCLUDEDIR) -I$(top_srcdir)/$(MPFR_INCLUDEDIR) -I$(top_srcdir)/$(MPC_INCLUDEDIR) -I$(top_srcdir)/$(QD_INCLUDEDIR)
libmpfrcxx_la_SOURCES = mpreal.cpp mpcomplex.cpp

endif
