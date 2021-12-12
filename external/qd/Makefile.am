QDVERSION=2.3.22
MAKEFLAGS =

CONFIGURE_ARGS=--prefix=$(abs_builddir)/../i/QD/ --host=$(host) --enable-fortran=no --enable-ieee-add=yes --enable-shared --enable-fma=yes $(QDLIB_FLAGS)
CONFIGURE_ARGS2=--prefix=$(prefix) --host=$(host) --enable-fortran=no --enable-ieee-add=yes --enable-shared --enable-fma=yes $(QDLIB_FLAGS)

EXTRACT_STAMP=.qd_extract_done
PATCH_STAMP=.qd_patch_done
CONFIGURE_STAMP=.qd_configure_done
BUILD_STAMP=.qd_build_done
INSTALL_STAMP=.qd_install_done
INSTALL_INTERNAL_STAMP=.qd_install_internal_done
CHECK_STAMP=.qd_check_done

all: extract patch config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal  ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfz ../../download/qd-$(QDVERSION).tar.gz
	cd $(abs_builddir)/work/install; tar xvfz ../../download/qd-$(QDVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-mingw
	cd $(abs_builddir)/work/install/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-mingw
endif
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-configure
	cd $(abs_builddir)/work/install/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-configure
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-sqrtoverflow
	cd $(abs_builddir)/work/install/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-sqrtoverflow

	cd $(abs_builddir)/work/install/qd-$(QDVERSION); cp $(abs_builddir)/patches/config.sub config
	cd $(abs_builddir)/work/install/qd-$(QDVERSION); cp $(abs_builddir)/patches/config.guess config
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); cp $(abs_builddir)/patches/config.sub config
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); cp $(abs_builddir)/patches/config.guess config
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
if IS_INTELCC
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; env "CC=$(CC) -fp-model precise" "CXX=$(CXX) -fp-model precise" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; env "CC=$(CC) -fp-model precise" "CXX=$(CXX) -fp-model precise" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
else
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ;  env "CC=$(CC)" "CXX=$(CXX)" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; TESTS_ENVIRONMENT=wine64 ${MAKE} check
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ;  TESTS_ENVIRONMENT=wine64 ${MAKE} check
else
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE} check
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; ${MAKE} check
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; ${MAKE} install
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
