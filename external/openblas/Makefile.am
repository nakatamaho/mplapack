OPENBLASVERSION=0.3.21
MAKEFLAGS =

EXTRACT_STAMP=.openblas_extract_done
PATCH_STAMP=.openblas_patch_done
BUILD_STAMP=.openblas_build_done
INSTALL_STAMP=.openblas_install_done
INSTALL_INTERNAL_STAMP=.openblas_install_internal_done
CHECK_STAMP=.openblas_check_done

# openblas does't allow to use ccache
_CC=echo $$(echo "${CC}" | awk '{print $$NF}' )
_FC=echo $$(echo "${FC}" | awk '{print $$NF}' )

all: extract build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/install
	cd $(abs_builddir)/work/install; tar xvfz ../../download/OpenBLAS-$(OPENBLASVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; patch -p0 < $(abs_builddir)/patches/patch-mingw
endif
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(BUILD_STAMP): work/$(PATCH_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} BINARY=64 CC=`$(_CC)` FC=`$(_FC)` DYNAMIC_ARCH=1 HOSTCC=gcc
else
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} BINARY=64 CC=`$(_CC)` FC=`$(_FC)` DYNAMIC_ARCH=1
endif
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if !IS_MINGW
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} -j lapack-test
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} PREFIX=$(abs_builddir)/../i/OPENBLAS/ install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} PREFIX=$(prefix) install
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

