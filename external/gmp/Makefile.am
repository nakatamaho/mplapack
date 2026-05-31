GMPVERSION=6.3.0
GMP_DISTFILE=gmp/download/gmp-$(GMPVERSION).tar.xz
MAKEFLAGS =

if IS_MACOS
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix)
else
if IS_MINGW
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/ --build=x86_64-linux-gnu --host=$(host) --enable-shared --disable-static
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix) --build=x86_64-linux-gnu --host=$(host) --enable-shared --disable-static
else
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/ --host=$(host)
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix) --host=$(host)
endif
endif

EXTRACT_STAMP=.gmp_extract_done
PATCH_STAMP=.gmp_patch_done
CONFIGURE_STAMP=.gmp_configure_done
BUILD_STAMP=.gmp_build_done
CHECK_STAMP=.gmp_check_done
INSTALL_STAMP=.gmp_install_done
INSTALL_INTERNAL_STAMP=.gmp_install_internal_done

all: extract config build check install_internal

.PHONY: verify-distfile
verify-distfile:
	@distfile="$(GMP_DISTFILE)"; \
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
	@$(MAKE) verify-distfile GMPVERSION="$(GMPVERSION)" GMP_DISTFILE="$(GMP_DISTFILE)"
	rm -rf work ; mkdir -p work/internal ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfJ $(abs_srcdir)/download/gmp-$(GMPVERSION).tar.xz
	cd $(abs_builddir)/work/install; tar xvfJ $(abs_srcdir)/download/gmp-$(GMPVERSION).tar.xz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
#https://trac.macports.org/ticket/70859
if IS_MACOS
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) && patch -p0 < $(abs_srcdir)/patches/gmp-6.3.0_libtool-2.4.7.patch
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION)  && patch -p0 < $(abs_srcdir)/patches/gmp-6.3.0_libtool-2.4.7.patch
endif
#https://github.com/compiler-explorer/gcc-cross-builder/blob/main/build/patches/crosstool-ng/0001-gcc-15-fix.patch
	if test "@GMP_NEEDS_GCC15_PATCH@" = yes; then \
		cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) && patch -p0 < $(abs_srcdir)/patches/gcc-15.patch && \
		cd $(abs_builddir)/work/install/gmp-$(GMPVERSION)  && patch -p0 < $(abs_srcdir)/patches/gcc-15.patch; \
	fi
	touch $(abs_builddir)/work/$(PATCH_STAMP)


work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
if IS_MACOS
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CFLAGS=$(CFLAGS)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CFLAGS=$(CFLAGS)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
else
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CFLAGS=$(CFLAGS)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CFLAGS=$(CFLAGS)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE}
	+cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if IS_MINGW
#https://gmplib.org/list-archives/gmp-bugs/2020-January/004720.html
#https://gmplib.org/list-archives/gmp-bugs/2020-January/004721.html
	+cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} check LOG_COMPILER=wine
	+cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} check LOG_COMPILER=wine
else
	+cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} check
	+cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} check
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	+cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	+cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

uninstall:
	-+cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} uninstall

distdir:

clean-local:
	rm -rf work

distclean-local:
	rm -rf work
