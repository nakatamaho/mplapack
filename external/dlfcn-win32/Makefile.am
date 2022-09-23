DLFCN_WIN32VERSION=1.3.1
MAKEFLAGS =

EXTRACT_STAMP=.dlfcn-win32_extract_done
PATCH_STAMP=.dlfcn-win32_patch_done
CONFIGURE_STAMP=.dlfcn-win32_configure_done
BUILD_STAMP=.dlfcn-win32_build_done
INSTALL_STAMP=.dlfcn-win32_install_done
CHECK_STAMP=.dlfcn-win32_check_done

all: extract patch config build check install

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/install
	cd $(abs_builddir)/work/install; tar xvfz ../../download/dlfcn-win32-$(DLFCN_WIN32VERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION); ./configure --enable-wine=wine64 --cross-prefix=x86_64-w64-mingw32- --prefix=$(abs_builddir)/../i/DLFCN_WIN32/
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION); make
	cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION)
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/install/dlfcn-win32-$(DLFCN_WIN32VERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install: work/$(INSTALL_STAMP)

distdir:

clean:
	rm -rf work
distclean:
	rm -rf work
