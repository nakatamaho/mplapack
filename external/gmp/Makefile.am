GMPVERSION=6.2.1
MAKEFLAGS =

if IS_MACOSX
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
#https://gmplib.org/list-archives/gmp-bugs/2020-January/004718.html
if IS_MINGW
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; env "CC_FOR_BUILD=x86_64-linux-gnu-gcc" "CC=$(CC) -D__USE_MINGW_ANSI_STDIO" "CXX=$(CXX)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS) #; sed -i 's/#define HAVE_VSNPRINTF 1//g' config.h
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; env "CC_FOR_BUILD=x86_64-linux-gnu-gcc" "CC=$(CC) -D__USE_MINGW_ANSI_STDIO" "CXX=$(CXX)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2) #; sed -i 's/#define HAVE_VSNPRINTF 1//g' config.h
else
#https://gmplib.org/list-archives/gmp-bugs/2014-September/003526.html
if IS_INTELCC
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; env "CC=$(CC) -fp-model precise" "CXX=$(CXX) -fp-model precise" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; env "CC=$(CC) -fp-model precise" "CXX=$(CXX) -fp-model precise" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
else
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS)
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS2)
endif
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)
if IS_MINGW
#https://gmplib.org/list-archives/gmp-bugs/2020-January/004720.html
#https://gmplib.org/list-archives/gmp-bugs/2020-January/004721.html
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} check LOG_COMPILER=wine64
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION)  ; ${MAKE} check LOG_COMPILER=wine64
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

