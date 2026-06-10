MPFRVERSION=4.2.2
MPFR_DISTFILE=mpfr/download/mpfr-$(MPFRVERSION).tar.bz2
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

.PHONY: verify-distfile
verify-distfile:
	@distfile="$(MPFR_DISTFILE)"; \
	cd "$(abs_srcdir)/.." && \
	entry=`awk -v f="$$distfile" '$$2 == f { print $$0 }' distfiles.sha256`; \
	if test -z "$$entry"; then \
		echo "error: no SHA256 checksum recorded for $$distfile" >&2; \
		exit 1; \
	fi; \
	if command -v sha256sum >/dev/null 2>&1; then \
		printf '%s\n' "$$entry" | sha256sum -c -; \
	elif command -v shasum >/dev/null 2>&1; then \
		printf '%s\n' "$$entry" | shasum -a 256 -c -; \
	else \
		echo "error: no SHA256 checksum tool found (sha256sum or shasum)" >&2; \
		exit 1; \
	fi

work/$(EXTRACT_STAMP):
	@$(MAKE) verify-distfile MPFRVERSION="$(MPFRVERSION)" MPFR_DISTFILE="$(MPFR_DISTFILE)"
	rm -rf work ; mkdir -p work/internal  ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfj $(abs_srcdir)/download/mpfr-$(MPFRVERSION).tar.bz2
	cd $(abs_builddir)/work/install; tar xvfj $(abs_srcdir)/download/mpfr-$(MPFRVERSION).tar.bz2
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	@# cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; patch -N -Z -p1 < ../../../download/allpatches
	@# cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; patch -N -Z -p1 < ../../../download/allpatches
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; env "CC=$(CC) -D__USE_MINGW_ANSI_STDIO" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS) --enable-shared --disable-static
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ;  env "CC=$(CC) -D__USE_MINGW_ANSI_STDIO" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2) --enable-shared --disable-static
else
	cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath,$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ;  env "CC=$(CC)" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath,$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; ${MAKE}
	+cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)
if IS_MINGW
	+cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; VERBOSE=1 ${MAKE} check LOG_COMPILER=wine
	+cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)  ; VERBOSE=1 ${MAKE} check LOG_COMPILER=wine
else
	+cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; ${MAKE} check
	+cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE} check
endif

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	+cd $(abs_builddir)/work/internal/mpfr-$(MPFRVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
if IS_MACOS
if !SYSTEM_GMP
	ORG=`otool -L $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib  | grep gmp | awk '{print $$1}'` ;\
	_FIX=`otool -L $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib | grep gmp | awk '{print $$1}' | sed 's|/| |g' | awk '{print $$NF}'` ; FIX=`echo $(prefix)/lib/$$_FIX` ;\
	echo "install_name_tool -change $$ORG $$FIX $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib" ;\
	install_name_tool -change $$ORG $$FIX $(abs_builddir)/work/install/mpfr-$(MPFRVERSION)/src/.libs/libmpfr.dylib
endif
endif
	+cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

uninstall:
	-+cd $(abs_builddir)/work/install/mpfr-$(MPFRVERSION) ; ${MAKE} uninstall

distdir:

clean-local:
	rm -rf work

distclean-local:
	rm -rf work
