GMPVERSION=6.1.2
MAKEFLAGS =

if IS_MACOSX
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix)
else
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/ --host=$(host)
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix) --host=$(host)
endif

EXTRACT_STAMP=.gmp_extract_done
PATCH_STAMP=.gmp_patch_done
CONFIGURE_STAMP=.gmp_configure_done
BUILD_STAMP=.gmp_build_done
INSTALL_STAMP=.gmp_install_done
INSTALL_INTERNAL_STAMP=.gmp_install_internal_done
CHECK_STAMP=.gmp_check_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal ; mkdir -p work/install
	cd $(abs_builddir)/work/internal; tar xvfj ../../download/gmp-$(GMPVERSION).tar.bz2
	cd $(abs_builddir)/work/install; tar xvfj ../../download/gmp-$(GMPVERSION).tar.bz2
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)
if IS_MINGW
#	somehow it fails
#	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; TESTS_ENVIRONMENT=wine ${MAKE} check
#	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; TESTS_ENVIRONMENT=wine ${MAKE} check
else
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} check
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} check
endif

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} install
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

