OPENBLASVERSION=0.3.32
OPENBLAS_DISTFILE=openblas/download/OpenBLAS-$(OPENBLASVERSION).tar.gz
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

.PHONY: verify-distfile
verify-distfile:
	@distfile="$(OPENBLAS_DISTFILE)"; \
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
	@$(MAKE) verify-distfile OPENBLASVERSION="$(OPENBLASVERSION)" OPENBLAS_DISTFILE="$(OPENBLAS_DISTFILE)"
	rm -rf work ; mkdir -p work/install
	cd $(abs_builddir)/work/install; tar xvfz $(abs_srcdir)/download/OpenBLAS-$(OPENBLASVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; patch -p1 < $(abs_srcdir)/patches/patch-mingw-laed3
endif
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(BUILD_STAMP): work/$(PATCH_STAMP)
if IS_MINGW
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} BINARY=64 CC=`$(_CC)` FC=`$(_FC)` DYNAMIC_ARCH=1 HOSTCC=gcc
else
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} BINARY=64 CC=`$(_CC)` FC=`$(_FC)` DYNAMIC_ARCH=1
endif
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if !IS_MINGW
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} -j lapack-test
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
if IS_MINGW
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} CC=`$(_CC)` FC=`$(_FC)` DYNAMIC_ARCH=1 HOSTCC=gcc PREFIX=$(abs_builddir)/../i/OPENBLAS/ install
else
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} PREFIX=$(abs_builddir)/../i/OPENBLAS/ install
endif
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
if IS_MINGW
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} CC=`$(_CC)` FC=`$(_FC)` DYNAMIC_ARCH=1 HOSTCC=gcc PREFIX=$(prefix) install
else
	+cd $(abs_builddir)/work/install/OpenBLAS-$(OPENBLASVERSION) ; ${MAKE} PREFIX=$(prefix) install
endif
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

uninstall:
	-rm -f  $(DESTDIR)$(prefix)/lib/libopenblas*.a
	-rm -f  $(DESTDIR)$(prefix)/lib/libopenblas*.so*
	-rm -f  $(DESTDIR)$(prefix)/lib/libopenblas*.dylib
	-rm -f  $(DESTDIR)$(prefix)/lib/pkgconfig/openblas.pc
	-rm -rf $(DESTDIR)$(prefix)/lib/cmake/openblas
	-rm -f  $(DESTDIR)$(prefix)/include/cblas.h
	-rm -f  $(DESTDIR)$(prefix)/include/lapacke.h
	-rm -f  $(DESTDIR)$(prefix)/include/lapacke_config.h
	-rm -f  $(DESTDIR)$(prefix)/include/lapacke_mangling.h
	-rm -f  $(DESTDIR)$(prefix)/include/lapacke_utils.h
	-rm -f  $(DESTDIR)$(prefix)/include/lapack.h
	-rm -f  $(DESTDIR)$(prefix)/include/openblas_config.h
	-rm -f  $(DESTDIR)$(prefix)/include/f77blas.h
	-rm -f  $(DESTDIR)$(prefix)/bin/libopenblas*.dll

distdir:

clean-local:
	rm -rf work
distclean-local:
	rm -rf work
