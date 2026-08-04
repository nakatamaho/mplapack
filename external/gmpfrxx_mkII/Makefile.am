GMPFRXX_MKII_VERSION = 1.2.0
GMPFRXX_MKII_DISTFILE = gmpfrxx_mkII.$(GMPFRXX_MKII_VERSION).tar.xz
GMPFRXX_MKII_SOURCE_DIR = gmpfrxx_mkII.$(GMPFRXX_MKII_VERSION)
GMPFRXX_MKII_INTERNAL_PREFIX = $(abs_builddir)/../i/GMPFRXX_MKII
GMPFRXX_MKII_WORK = $(abs_builddir)/work
GMPFRXX_MKII_CHECK_BUILD = $(GMPFRXX_MKII_WORK)/check
GMPFRXX_MKII_OPTIONS = -DBUILD_TESTING=ON -DGMPFRXX_MKII_DEPS_AUTO_FETCH=OFF -DGMPFRXX_MKII_BUILD_EXAMPLES=OFF -DGMPFRXX_MKII_BUILD_BENCHMARKS=OFF -DGMPFRXX_MKII_COMPONENTS=GMP,MPFR,MPC
GMPFRXX_MKII_DEPS = -DCMAKE_PREFIX_PATH="$(abs_builddir)/../i/GMP;$(abs_builddir)/../i/MPFR;$(abs_builddir)/../i/MPC"
if ENABLE_GMP
GMPFRXX_MKII_UNINSTALL_ENABLED = yes
else
if ENABLE_MPFR
GMPFRXX_MKII_UNINSTALL_ENABLED = yes
else
GMPFRXX_MKII_UNINSTALL_ENABLED = no
endif
endif

all: $(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done
$(GMPFRXX_MKII_WORK)/.verify_done:
	@set -e; \
	distfile="gmpfrxx_mkII/download/$(GMPFRXX_MKII_DISTFILE)"; \
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
	fi; \
	mkdir -p $(GMPFRXX_MKII_WORK); \
	touch $@
$(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done: $(GMPFRXX_MKII_WORK)/.verify_done
	rm -rf $(GMPFRXX_MKII_INTERNAL_PREFIX)
	mkdir -p $(GMPFRXX_MKII_INTERNAL_PREFIX)
	tar -xJf $(abs_srcdir)/download/$(GMPFRXX_MKII_DISTFILE) -C $(GMPFRXX_MKII_INTERNAL_PREFIX) --strip-components=1
	cd $(GMPFRXX_MKII_INTERNAL_PREFIX) && patch -p0 < $(abs_srcdir)/patches/patch-binary128-test-mpfr-float128
	touch $@
check-local: $(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done
	@set -e; \
	cmake_cxx=; \
	cmake_cxx_launcher=; \
	for word in $(CXX); do \
		case "$$word" in \
			ccache|*/ccache) cmake_cxx_launcher=ccache ;; \
			-*) ;; \
			*) cmake_cxx="$$word"; break ;; \
		esac; \
	done; \
	if test -z "$$cmake_cxx"; then cmake_cxx=c++; fi; \
	cmake_cxx_launcher_arg=; \
	if test -n "$$cmake_cxx_launcher"; then \
		cmake_cxx_launcher_arg="-DCMAKE_CXX_COMPILER_LAUNCHER=$$cmake_cxx_launcher"; \
	fi; \
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
	env PATH="$$cmake_path" CXX="$$cmake_cxx" CXXFLAGS="$(CXXFLAGS)" LDFLAGS="$(LDFLAGS)" cmake \
		-S $(GMPFRXX_MKII_INTERNAL_PREFIX) \
		-B $(GMPFRXX_MKII_CHECK_BUILD) \
		-DCMAKE_CXX_COMPILER="$$cmake_cxx" \
		-DCMAKE_CXX_FLAGS="$(CXXFLAGS)" \
		-DCMAKE_EXE_LINKER_FLAGS="$(LDFLAGS)" \
		$$cmake_cxx_launcher_arg \
		$(GMPFRXX_MKII_OPTIONS) \
		$(GMPFRXX_MKII_DEPS)
	+cmake --build $(GMPFRXX_MKII_CHECK_BUILD)
	ctest --test-dir $(GMPFRXX_MKII_CHECK_BUILD) --output-on-failure
clean-local distclean-local:
	rm -rf $(GMPFRXX_MKII_WORK) $(GMPFRXX_MKII_INTERNAL_PREFIX)
install-data-local: $(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done
	$(MKDIR_P) $(DESTDIR)$(includedir)/mplapack
	cp -R $(GMPFRXX_MKII_INTERNAL_PREFIX)/include/. $(DESTDIR)$(includedir)/mplapack/

uninstall-local:
	@if test x"$(GMPFRXX_MKII_UNINSTALL_ENABLED)" = xyes; then \
		rm -rf $(DESTDIR)$(includedir)/mplapack/gmpfrxx_mkII; \
		rm -f $(DESTDIR)$(includedir)/mplapack/gmpfrxx_mkII.h; \
		rm -f $(DESTDIR)$(includedir)/mplapack/gmpxx_mkII.h; \
		rm -f $(DESTDIR)$(includedir)/mplapack/mpfrxx_mkII.h; \
		rm -f $(DESTDIR)$(includedir)/mplapack/mpcxx_mkII.h; \
	fi
EXTRA_DIST = \
	download/$(GMPFRXX_MKII_DISTFILE) \
	patches/patch-binary128-test-mpfr-float128
