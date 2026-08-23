QDVERSION=1.4.0
QD_DISTFILE=qd/download/libQD3-$(QDVERSION).tar.gz
QD_SRCDIR=libQD3-$(QDVERSION)
QD_INTERNAL_BUILDDIR=$(abs_builddir)/work/internal-build
QD_INSTALL_BUILDDIR=$(abs_builddir)/work/install-build
QD_INTERNAL_LIBDIR=$(abs_builddir)/../i/QD/lib
QD_INSTALL_LIBDIR=$(libdir)
MAKEFLAGS =

# Ensure dist tarball contains the vendor tarball.
EXTRA_DIST = \
	download/libQD3-$(QDVERSION).tar.gz

if IS_MINGW
QD_BUILD_STATIC=OFF
qd_edd_real_arg=-DQD_ENABLE_EDD_REAL=OFF
qd_internal_rpath_ldflags=
qd_install_rpath_ldflags=
else
QD_BUILD_STATIC=ON
qd_edd_real_arg=
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
	$(qd_edd_real_arg) \
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
	$(qd_edd_real_arg) \
	$(QDLIB_FLAGS)

EXTRACT_STAMP=.qd_extract_done
CONFIGURE_STAMP=.qd_configure_done
BUILD_STAMP=.qd_build_done
INSTALL_INTERNAL_STAMP=.qd_install_internal_done
CHECK_STAMP=.qd_check_done

all: extract config build check install_internal

.PHONY: verify-distfile install uninstall
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
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(CONFIGURE_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	@set -e; \
	cmake_cc=; cmake_cxx=; \
	cmake_cc_launcher=; cmake_cxx_launcher=; \
	cmake_cflags_extra=; cmake_cxxflags_extra=; \
	for word in $(CC); do \
		case "$$word" in \
			ccache|*/ccache) cmake_cc_launcher=ccache ;; \
			-*) cmake_cflags_extra="$$cmake_cflags_extra $$word" ;; \
			*) if test -z "$$cmake_cc"; then cmake_cc="$$word"; else cmake_cflags_extra="$$cmake_cflags_extra $$word"; fi ;; \
		esac; \
	done; \
	for word in $(CXX); do \
		case "$$word" in \
			ccache|*/ccache) cmake_cxx_launcher=ccache ;; \
			-*) cmake_cxxflags_extra="$$cmake_cxxflags_extra $$word" ;; \
			*) if test -z "$$cmake_cxx"; then cmake_cxx="$$word"; else cmake_cxxflags_extra="$$cmake_cxxflags_extra $$word"; fi ;; \
		esac; \
	done; \
	if test -z "$$cmake_cc"; then cmake_cc=x86_64-w64-mingw32-gcc-posix; fi; \
	if test -z "$$cmake_cxx"; then cmake_cxx=x86_64-w64-mingw32-g++-posix; fi; \
	cmake_system_processor=`"$$cmake_cc" -dumpmachine 2>/dev/null | sed 's/-.*//'`; \
	if test -z "$$cmake_system_processor"; then cmake_system_processor=x86_64; fi; \
	cmake_cc_launcher_arg=; \
	if test -n "$$cmake_cc_launcher"; then cmake_cc_launcher_arg="-DCMAKE_C_COMPILER_LAUNCHER=$$cmake_cc_launcher"; fi; \
	cmake_cxx_launcher_arg=; \
	if test -n "$$cmake_cxx_launcher"; then cmake_cxx_launcher_arg="-DCMAKE_CXX_COMPILER_LAUNCHER=$$cmake_cxx_launcher"; fi; \
	cmake_path=; \
	old_ifs=$$IFS; \
	IFS=:; \
	for dir in $$PATH; do \
		case "$$dir" in \
			*/ccache) ;; \
			*) if test -z "$$cmake_path"; then cmake_path="$$dir"; else cmake_path="$$cmake_path:$$dir"; fi ;; \
		esac; \
	done; \
	IFS=$$old_ifs; \
	QD_CFLAGS="$(CFLAGS)$$cmake_cflags_extra"; \
	QD_CXXFLAGS="$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@$$cmake_cxxflags_extra"; \
	env PATH="$$cmake_path" CC="$$cmake_cc" CXX="$$cmake_cxx" CFLAGS="$$QD_CFLAGS" CXXFLAGS="$$QD_CXXFLAGS" LDFLAGS="$(LDFLAGS) $(qd_internal_rpath_ldflags)" NM="$(NM)" RANLIB="$(RANLIB)" cmake \
		-S $(abs_builddir)/work/internal/$(QD_SRCDIR) \
		-B $(QD_INTERNAL_BUILDDIR) \
		$(CMAKE_ARGS) \
		-DCMAKE_SYSTEM_NAME=Windows \
		-DCMAKE_SYSTEM_PROCESSOR="$$cmake_system_processor" \
		-DCMAKE_C_COMPILER="$$cmake_cc" \
		-DCMAKE_CXX_COMPILER="$$cmake_cxx" \
		$$cmake_cc_launcher_arg \
		$$cmake_cxx_launcher_arg \
		-DCMAKE_CROSSCOMPILING_EMULATOR=wine; \
	env PATH="$$cmake_path" CC="$$cmake_cc" CXX="$$cmake_cxx" CFLAGS="$$QD_CFLAGS" CXXFLAGS="$$QD_CXXFLAGS" LDFLAGS="$(LDFLAGS) $(qd_install_rpath_ldflags)" NM="$(NM)" RANLIB="$(RANLIB)" cmake \
		-S $(abs_builddir)/work/install/$(QD_SRCDIR) \
		-B $(QD_INSTALL_BUILDDIR) \
		$(CMAKE_ARGS2) \
		-DCMAKE_SYSTEM_NAME=Windows \
		-DCMAKE_SYSTEM_PROCESSOR="$$cmake_system_processor" \
		-DCMAKE_C_COMPILER="$$cmake_cc" \
		-DCMAKE_CXX_COMPILER="$$cmake_cxx" \
		$$cmake_cc_launcher_arg \
		$$cmake_cxx_launcher_arg \
		-DCMAKE_CROSSCOMPILING_EMULATOR=wine
else
if IS_INTELCC
	QD_CFLAGS="$${CFLAGS// -xHost/} -O2 -fp-model strict" ; \
	QD_CXXFLAGS="$${CXXFLAGS// -xHost/} -O2 -fp-model strict" ; \
	env "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$$QD_CFLAGS" "CXXFLAGS=$$QD_CXXFLAGS" "LDFLAGS=$(LDFLAGS) $(qd_internal_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/internal/$(QD_SRCDIR) -B $(QD_INTERNAL_BUILDDIR) $(CMAKE_ARGS) ; \
	env "CC=$(CC)" "CXX=$(CXX)" "CFLAGS=$$QD_CFLAGS" "CXXFLAGS=$$QD_CXXFLAGS" "LDFLAGS=$(LDFLAGS) $(qd_install_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/install/$(QD_SRCDIR) -B $(QD_INSTALL_BUILDDIR) $(CMAKE_ARGS2)
else
	env "CC=$(CC)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@" "LDFLAGS=$(LDFLAGS) $(qd_internal_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/internal/$(QD_SRCDIR) -B $(QD_INTERNAL_BUILDDIR) $(CMAKE_ARGS)
	env "CC=$(CC)" "CXX=$(CXX)" "CXXFLAGS=$(CXXFLAGS) @QD_EXTRA_CXXFLAGS@" "LDFLAGS=$(LDFLAGS) $(qd_install_rpath_ldflags)" "NM=$(NM)" "RANLIB=$(RANLIB)" cmake -S $(abs_builddir)/work/install/$(QD_SRCDIR) -B $(QD_INSTALL_BUILDDIR) $(CMAKE_ARGS2)
endif
endif
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cmake --build $(QD_INTERNAL_BUILDDIR)
	+cmake --build $(QD_INSTALL_BUILDDIR)
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
if IS_MINGW
	+env WINEPATH="$(QD_INTERNAL_BUILDDIR);$${WINEPATH}" ctest --test-dir $(QD_INTERNAL_BUILDDIR) --output-on-failure
	+env WINEPATH="$(QD_INSTALL_BUILDDIR);$${WINEPATH}" ctest --test-dir $(QD_INSTALL_BUILDDIR) --output-on-failure
else
	+ctest --test-dir $(QD_INTERNAL_BUILDDIR) --output-on-failure
	+ctest --test-dir $(QD_INSTALL_BUILDDIR) --output-on-failure
endif
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	+cmake --install $(QD_INTERNAL_BUILDDIR)
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(CHECK_STAMP)
	+cmake -S $(abs_builddir)/work/install/$(QD_SRCDIR) -B $(QD_INSTALL_BUILDDIR) $(CMAKE_ARGS2)
	+cmake --build $(QD_INSTALL_BUILDDIR)
	+env DESTDIR="$(DESTDIR)" cmake --install $(QD_INSTALL_BUILDDIR)

uninstall:
	-+cmake -S $(abs_builddir)/work/install/$(QD_SRCDIR) -B $(QD_INSTALL_BUILDDIR) $(CMAKE_ARGS2)
	-+cmake --build $(QD_INSTALL_BUILDDIR) --target uninstall
clean-local:
	rm -rf work
distclean-local:
	rm -rf work
