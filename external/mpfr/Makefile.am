MPFRVERSION=4.1.0
MAKEFLAGS =

#GMPINCLUDEDIR=$(abs_builddir)/../i/GMP/include/
#GMPLIBDIR=$(abs_builddir)/../i/GMP/lib/

CONFIGURE_ARGS=--prefix=$(abs_builddir)/../i/MPFR/ --host=$(host) --with-gmp-include=@GMP_INCLUDEDIR@ --with-gmp-lib=@GMP_LIBDIR@
CONFIGURE_ARGS2=--prefix=$(prefix) --host=$(host) --with-gmp-include=@GMP_INCLUDEDIR@ --with-gmp-lib=@GMP_LIBDIR@

EXTRACT_STAMP=.mpfr_extract_done
PATCH_STAMP=.mpfr_patch_done
CONFIGURE_STAMP=.mpfr_configure_done
BUILD_STAMP=.mpfr_build_done
INSTALL_STAMP=.mpfr_install_done
INSTALL_INTERNAL_STAMP=.mpfr_install_internal_done
CHECK_STAMP=.mpfr_check_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal  ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfj ../../download/mpfr-$(MPFRVERSION).tar.bz2
	cd $(abs_builddir)/work/install; tar xvfj ../../download/mpfr-$(MPFRVERSION).tar.bz2
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
#	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; patch -N -Z -p1 < ../../../download/allpatches
#	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; patch -N -Z -p1 < ../../../download/allpatches
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; env "CC=$(CC) -D__USE_MINGW_ANSI_STDIO" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS) --enable-static --disable-shared
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ;  env "CC=$(CC) -D__USE_MINGW_ANSI_STDIO" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2) --enable-static --disable-shared
else
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath,$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ;  env "CC=$(CC)" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath,$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; VERBOSE=1 ${MAKE} check LOG_COMPILER=wine64
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)  ; VERBOSE=1 ${MAKE} check LOG_COMPILER=wine64
else
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; ${MAKE} check
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE} check
endif

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	ORG=`otool -L $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib  | grep gmp | awk '{print $$1}'` ;\
	_FIX=`otool -L $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib | grep gmp | awk '{print $$1}' | sed 's|/| |g' | awk '{print $$NF}'` ; FIX=`echo $(prefix)/lib/$$_FIX` ;\
	echo "install_name_tool -change $$ORG $$FIX $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib" ;\
	install_name_tool -change $$ORG $$FIX $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

distdir:

clean:
	rm -rf work

distclean:
	rm -rf work
