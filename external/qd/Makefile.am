QDVERSION=2.3.24
QD_DISTFILE=qd/download/qd-$(QDVERSION).tar.gz
MAKEFLAGS =

# Ensure dist tarball contains the vendor tarball and our patch set
EXTRA_DIST = \
	download/qd-$(QDVERSION).tar.gz \
	patches/patch-overflow \
	patches/patch-no-doc \
	patches/patch-mingw

if IS_MINGW
CONFIGURE_ARGS=--prefix=$(abs_builddir)/../i/QD/ --host=$(host) --enable-fortran=no --enable-ieee-add=yes --enable-sloppy-mul=no --enable-sloppy-div=no --enable-shared=yes --enable-static=no $(QDLIB_FLAGS)
CONFIGURE_ARGS2=--prefix=$(prefix) --host=$(host) --enable-fortran=no --enable-ieee-add=yes --enable-sloppy-mul=no --enable-sloppy-div=no --enable-shared=yes --enable-static=no $(QDLIB_FLAGS)
else
CONFIGURE_ARGS=--prefix=$(abs_builddir)/../i/QD/ --host=$(host) --enable-fortran=no --enable-ieee-add=yes --enable-sloppy-mul=no --enable-sloppy-div=no --enable-shared $(QDLIB_FLAGS)
CONFIGURE_ARGS2=--prefix=$(prefix) --host=$(host) --enable-fortran=no --enable-ieee-add=yes --enable-sloppy-mul=no --enable-sloppy-div=no --enable-shared $(QDLIB_FLAGS)
endif

EXTRACT_STAMP=.qd_extract_done
PATCH_STAMP=.qd_patch_done
CONFIGURE_STAMP=.qd_configure_done
BUILD_STAMP=.qd_build_done
INSTALL_STAMP=.qd_install_done
INSTALL_INTERNAL_STAMP=.qd_install_internal_done
CHECK_STAMP=.qd_check_done

all: extract patch config build check install_internal

.PHONY: verify-distfile
verify-distfile:
	@distfile="$(QD_DISTFILE)"; \
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
	@$(MAKE) verify-distfile QDVERSION="$(QDVERSION)" QD_DISTFILE="$(QD_DISTFILE)"
	rm -rf work ; mkdir -p work/internal  ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfz $(abs_srcdir)/download/qd-$(QDVERSION).tar.gz
	cd $(abs_builddir)/work/install; tar xvfz $(abs_srcdir)/download/qd-$(QDVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/QD-$(QDVERSION); patch -p0 < $(abs_srcdir)/patches/patch-mingw ; autoreconf -f -i ; cp config/* .
	cd $(abs_builddir)/work/install/QD-$(QDVERSION); patch -p0 < $(abs_srcdir)/patches/patch-mingw ; autoreconf -f -i ; cp config/* .
endif
	cd $(abs_builddir)/work/internal/QD-$(QDVERSION); patch -p0 < $(abs_srcdir)/patches/patch-overflow
	cd $(abs_builddir)/work/install/QD-$(QDVERSION); patch -p0 < $(abs_srcdir)/patches/patch-overflow
	cd $(abs_builddir)/work/internal/QD-$(QDVERSION); patch -p0 < $(abs_srcdir)/patches/patch-no-doc
	cd $(abs_builddir)/work/install/QD-$(QDVERSION); patch -p0 < $(abs_srcdir)/patches/patch-no-doc

	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
if IS_INTELCC
	QD_CFLAGS="$${CFLAGS// -xHost/} -O2 -fp-model strict" ; \
	QD_CXXFLAGS="$${CXXFLAGS// -xHost/} -O2 -fp-model strict" ; \
	cd $(abs_builddir)/work/internal/QD-$(QDVERSION) ; autoreconf -f -i ; \
	env "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$$QD_CFLAGS" "CXXFLAGS=$$QD_CXXFLAGS" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS) ; \
	cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; autoreconf -f -i ; \
	env "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$$QD_CFLAGS" "CXXFLAGS=$$QD_CXXFLAGS" "FC=$(FC)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath=$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
else
	cd $(abs_builddir)/work/internal/QD-$(QDVERSION) ; autoreconf -f -i ; env "CC=$(CC)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; autoreconf -f -i ; env "CC=$(CC)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@" "FC=$(FC)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath,$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cd $(abs_builddir)/work/internal/QD-$(QDVERSION) ; ${MAKE}
	+cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if IS_MINGW
	+cd $(abs_builddir)/work/internal/QD-$(QDVERSION) ; ${MAKE} check LOG_COMPILER=wine
	+cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; ${MAKE} check LOG_COMPILER=wine
else
	+cd $(abs_builddir)/work/internal/QD-$(QDVERSION) ; ${MAKE} check
	+cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; ${MAKE} check
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	+cd $(abs_builddir)/work/internal/QD-$(QDVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	+cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

uninstall:
	-+cd $(abs_builddir)/work/install/QD-$(QDVERSION) ; ${MAKE} uninstall
clean-local:
	rm -rf work
distclean-local:
	rm -rf work
