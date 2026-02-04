GMPVERSION=6.3.0
MAKEFLAGS =

if IS_MACOS
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix)
GMP_CC_EXTRA =
GMP_CXX_EXTRA =
GMP_USE_LD_CLASSIC=
GMP_DISABLE_ASM =
GMP_CC_FOR_BUILD =
else
if IS_MINGW
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/ --build=x86_64-linux-gnu --host=$(host) --enable-shared --disable-static
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix) --build=x86_64-linux-gnu --host=$(host) --enable-shared --disable-static
GMP_CC_EXTRA=-D__USE_MINGW_ANSI_STDIO
GMP_CXX_EXTRA=
GMP_USE_LD_CLASSIC=
GMP_DISABLE_ASM=
GMP_CC_FOR_BUILD=x86_64-linux-gnu-gcc
else
CONFIGURE_ARGS=--enable-cxx --prefix=$(abs_builddir)/../i/GMP/ --host=$(host)
CONFIGURE_ARGS2=--enable-cxx --prefix=$(prefix) --host=$(host)
GMP_CC_EXTRA=
GMP_CXX_EXTRA=
GMP_USE_LD_CLASSIC=
GMP_DISABLE_ASM=
GMP_CC_FOR_BUILD=
endif
endif

if IS_INTELCC
GMP_CC_EXTRA += -fp-model precise
GMP_CXX_EXTRA += -fp-model precise
endif

EXTRACT_STAMP=.gmp_extract_done
PATCH_STAMP=.gmp_patch_done
CONFIGURE_STAMP=.gmp_configure_done
BUILD_STAMP=.gmp_build_done
CHECK_STAMP=.gmp_check_done
INSTALL_STAMP=.gmp_install_done
INSTALL_INTERNAL_STAMP=.gmp_install_internal_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal ; mkdir -p work/install
	cd $(abs_builddir)/work/internal ; tar xvf ../../download/gmp-$(GMPVERSION).tar.xz
	cd $(abs_builddir)/work/install ; tar xvf ../../download/gmp-$(GMPVERSION).tar.xz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

# References:
# mingw: https://gmplib.org/list-archives/gmp-bugs/2020-January/004718.html
# mingw: https://gmplib.org/list-archives/gmp-bugs/2014-September/003526.html
# macOS: https://trac.macports.org/ticket/70859
# macOS: https://gmplib.org/list-archives/gmp-bugs/2024-June/005505.html
# macOS: https://gmplib.org/list-archives/gmp-bugs/2024-July/005506.html
# macOS: https://gmplib.org/list-archives/gmp-bugs/2008-February/000930.html
# macOS: https://gmplib.org/list-archives/gmp-bugs/2008-February/000928.html
work/$(CONFIGURE_STAMP): work/$(EXTRACT_STAMP)
	@set -e; \
	cc="$(CC) $(GMP_CC_EXTRA)"; \
	cxx="$(CXX) $(GMP_CXX_EXTRA)"; \
	ldflags="$(LDFLAGS)"; \
	f77="$(F77)"; \
	nm="$(NM)"; \
	ranlib="$(RANLIB)"; \
	conf_extra=""; \
	cc_for_build="$(GMP_CC_FOR_BUILD)"; \
	if [ "$(GMP_USE_LD_CLASSIC)" = "1" ]; then \
	  ldflags="$$ldflags -Wl,-ld_classic"; \
	fi; \
	if [ "$(GMP_DISABLE_ASM)" = "1" ]; then \
	  conf_extra="$$conf_extra --disable-assembly"; \
	fi; \
	export CC="$$cc"; \
	export CXX="$$cxx"; \
	export F77="$$f77"; \
	export LDFLAGS="$$ldflags"; \
	export NM="$$nm"; \
	export RANLIB="$$ranlib"; \
	if [ -n "$$cc_for_build" ]; then \
	  export CC_FOR_BUILD="$$cc_for_build"; \
	fi; \
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ./configure $(CONFIGURE_ARGS) $$conf_extra; \
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ./configure $(CONFIGURE_ARGS2) $$conf_extra
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

# mingw check with wine:
# https://gmplib.org/list-archives/gmp-bugs/2020-January/004720.html
# https://gmplib.org/list-archives/gmp-bugs/2020-January/004721.html
work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} check LOG_COMPILER=wine
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} check LOG_COMPILER=wine
else
	cd $(abs_builddir)/work/internal/gmp-$(GMPVERSION) ; ${MAKE} check
	cd $(abs_builddir)/work/install/gmp-$(GMPVERSION) ; ${MAKE} check
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

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
