DLFCN_WIN32VERSION=1.4.2
DLFCN_WIN32_DISTFILE=dlfcn-win32/download/dlfcn-win32-$(DLFCN_WIN32VERSION).tar.gz
MAKEFLAGS =

EXTRACT_STAMP=.dlfcn-win32_extract_done
PATCH_STAMP=.dlfcn-win32_patch_done
CONFIGURE_STAMP=.dlfcn-win32_configure_done
BUILD_STAMP=.dlfcn-win32_build_done
INSTALL_STAMP=.dlfcn-win32_install_done
CHECK_STAMP=.dlfcn-win32_check_done

all: extract patch config build check install

.PHONY: verify-distfile
verify-distfile:
	@distfile="$(DLFCN_WIN32_DISTFILE)"; \
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
	@$(MAKE) verify-distfile DLFCN_WIN32VERSION="$(DLFCN_WIN32VERSION)" DLFCN_WIN32_DISTFILE="$(DLFCN_WIN32_DISTFILE)"
	rm -rf work ; mkdir -p work/install
	cd $(abs_builddir)/work/install; tar xvfz $(abs_srcdir)/download/dlfcn-win32-$(DLFCN_WIN32VERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION); ./configure --enable-wine=wine --cross-prefix=x86_64-w64-mingw32- --prefix=$(abs_builddir)/../i/DLFCN_WIN32/
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION); ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)


work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	+cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install: work/$(INSTALL_STAMP)

uninstall:
	-+cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION) ; ${MAKE} uninstall

distdir:

clean-local:
	rm -rf work
distclean-local:
	rm -rf work
