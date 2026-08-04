QDVERSION=1.2.0
QD_DISTFILE=qd/download/libQD3-$(QDVERSION).tar.gz
QD_SRCDIR=libQD3-$(QDVERSION)
QD_INTERNAL_BUILDDIR=$(abs_builddir)/work/internal-build
QD_INSTALL_BUILDDIR=$(abs_builddir)/work/install-build
QD_INTERNAL_LIBDIR=$(abs_builddir)/../i/QD/lib
QD_INSTALL_LIBDIR=$(libdir)
QD_PATCHES= \
	patches/patch-std-complex-interop \
	patches/patch-int-assignment-overloads \
	patches/patch-integer-conversion-helpers \
	patches/patch-cmake-uninstall-target \
	patches/patch-fpu-fix-80bit \
	patches/patch-qd-real-integer-arithmetic-overloads \
	patches/patch-real-integer-comparison-overloads \
	patches/patch-dd-td-edd-real-integer-arithmetic-overloads
MAKEFLAGS =

# Ensure dist tarball contains the vendor tarball.
EXTRA_DIST = \
	download/libQD3-$(QDVERSION).tar.gz \
	$(QD_PATCHES)

if IS_MINGW
QD_BUILD_STATIC=OFF
qd_internal_rpath_ldflags=
qd_install_rpath_ldflags=
else
QD_BUILD_STATIC=ON
qd_internal_rpath_ldflags=-Wl,-rpath,$(QD_INTERNAL_LIBDIR)
qd_install_rpath_ldflags=-Wl,-rpath,$(QD_INSTALL_LIBDIR)
endif

CMAKE_ARGS= \
	-DCMAKE_INSTALL_PREFIX=$(abs_builddir)/../i/QD/ \
	-DCMAKE_INSTALL_RPATH=$(QD_INTERNAL_LIBDIR) \
	-DBUILD_TESTING=ON \
	-DBUILD_SHARED_LIBS=ON \
	-DQD_BUILD_STATIC=$(QD_BUILD_STATIC) \
	-DQD_BUILD_FORTRAN=OFF \
	-DQD_ENABLE_IEEE_ADD=ON \
	-DQD_ENABLE_SLOPPY_MUL=OFF \
	-DQD_ENABLE_SLOPPY_DIV=OFF \
	$(QDLIB_FLAGS)
CMAKE_ARGS2= \
	-DCMAKE_INSTALL_PREFIX=$(prefix) \
	-DCMAKE_INSTALL_RPATH=$(QD_INSTALL_LIBDIR) \
	-DBUILD_TESTING=ON \
	-DBUILD_SHARED_LIBS=ON \
	-DQD_BUILD_STATIC=$(QD_BUILD_STATIC) \
	-DQD_BUILD_FORTRAN=OFF \
	-DQD_ENABLE_IEEE_ADD=ON \
	-DQD_ENABLE_SLOPPY_MUL=OFF \
	-DQD_ENABLE_SLOPPY_DIV=OFF \
	$(QDLIB_FLAGS)

EXTRACT_STAMP=.qd_extract_done
CONFIGURE_STAMP=.qd_configure_done
BUILD_STAMP=.qd_build_done
INSTALL_STAMP=.qd_install_done
INSTALL_INTERNAL_STAMP=.qd_install_internal_done
CHECK_STAMP=.qd_check_done

all: extract config build check install_internal

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
	cd $(abs_builddir)/work/internal; tar xvfz $(abs_srcdir)/download/libQD3-$(QDVERSION).tar.gz
	cd $(abs_builddir)/work/install; tar xvfz $(abs_srcdir)/download/libQD3-$(QDVERSION).tar.gz
	for patchfile in $(QD_PATCHES); do \
		cd $(abs_builddir)/work/internal/$(QD_SRCDIR) && patch -p1 < $(abs_srcdir)/$$patchfile; \
		cd $(abs_builddir)/work/install/$(QD_SRCDIR) && patch -p1 < $(abs_srcdir)/$$patchfile; \
	done
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(CONFIGURE_STAMP): work/$(EXTRACT_STAMP)
if IS_INTELCC
	QD_CFLAGS="$${CFLAGS// -xHost/} -O2 -fp-model strict" ; \
	QD_CXXFLAGS="$${CXXFLAGS// -xHost/} -O2 -fp-model strict" ; \
	env "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$$QD_CFLAGS" "CXXFLAGS=$$QD_CXXFLAGS" "LDFLAGS=$(LDFLAGS) $(qd_internal_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/internal/$(QD_SRCDIR) -B $(QD_INTERNAL_BUILDDIR) $(CMAKE_ARGS) ; \
	env "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$$QD_CFLAGS" "CXXFLAGS=$$QD_CXXFLAGS" "LDFLAGS=$(LDFLAGS) $(qd_install_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/install/$(QD_SRCDIR) -B $(QD_INSTALL_BUILDDIR) $(CMAKE_ARGS2)
else
	env "CC=$(CC)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@" "LDFLAGS=$(LDFLAGS) $(qd_internal_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/internal/$(QD_SRCDIR) -B $(QD_INTERNAL_BUILDDIR) $(CMAKE_ARGS)
	env "CC=$(CC)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@" "LDFLAGS=$(LDFLAGS) $(qd_install_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/install/$(QD_SRCDIR) -B $(QD_INSTALL_BUILDDIR) $(CMAKE_ARGS2)
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cmake --build $(QD_INTERNAL_BUILDDIR)
	+cmake --build $(QD_INSTALL_BUILDDIR)
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	+ctest --test-dir $(QD_INTERNAL_BUILDDIR)
	+ctest --test-dir $(QD_INSTALL_BUILDDIR)
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	+cmake --install $(QD_INTERNAL_BUILDDIR)
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	+cmake --install $(QD_INSTALL_BUILDDIR)
	touch $(abs_builddir)/work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

uninstall:
	-+cmake --build $(QD_INSTALL_BUILDDIR) --target uninstall
clean-local:
	rm -rf work
distclean-local:
	rm -rf work
