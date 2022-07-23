QDVERSION=2.3.23
MAKEFLAGS =

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

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal  ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfz ../../download/qd-$(QDVERSION).tar.gz
	cd $(abs_builddir)/work/install; tar xvfz ../../download/qd-$(QDVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-mingw ; autoreconf -f -i ; cp config/* .
	cd $(abs_builddir)/work/install/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-mingw ; autoreconf -f -i ; cp config/* .
endif
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-overflow
	cd $(abs_builddir)/work/install/qd-$(QDVERSION); patch -p0 < $(abs_builddir)/patches/patch-overflow

	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
if IS_INTELCC
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; env "CC=$(CC) -fp-model precise" "CXX=$(CXX) -fp-model precise" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; env "CC=$(CC) -fp-model precise" "CXX=$(CXX) -fp-model precise" "FC=$(FC)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath=$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
else
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "FC=$(FC)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ;  env "CC=$(CC)" "CXX=$(CXX)" "FC=$(FC)" "LDFLAGS=$(LDFLAGS) -Wl,-rpath,$(libdir)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE} check LOG_COMPILER=wine64
	cd $(abs_builddir)/work/install/qd-$(QDVERSION) ; ${MAKE} check LOG_COMPILER=wine64
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
