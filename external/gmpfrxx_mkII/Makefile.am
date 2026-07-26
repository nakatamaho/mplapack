GMPFRXX_MKII_VERSION = 1.1.0
GMPFRXX_MKII_DISTFILE = gmpfrxx_mkII.$(GMPFRXX_MKII_VERSION).tar.xz
GMPFRXX_MKII_SOURCE_DIR = gmpfrxx_mkII.$(GMPFRXX_MKII_VERSION)
GMPFRXX_MKII_INTERNAL_PREFIX = $(abs_builddir)/../i/GMPFRXX_MKII
GMPFRXX_MKII_WORK = $(abs_builddir)/work
GMPFRXX_MKII_BUILD_JOBS = 32
GMPFRXX_MKII_CHECK_BUILD = $(GMPFRXX_MKII_WORK)/check
GMPFRXX_MKII_OPTIONS = -DBUILD_TESTING=ON -DGMPFRXX_MKII_DEPS_AUTO_FETCH=OFF -DGMPFRXX_MKII_BUILD_EXAMPLES=OFF -DGMPFRXX_MKII_BUILD_BENCHMARKS=OFF -DGMPFRXX_MKII_COMPONENTS=GMP,MPFR,MPC
GMPFRXX_MKII_DEPS = -DCMAKE_PREFIX_PATH="$(abs_builddir)/../i/GMP;$(abs_builddir)/../i/MPFR;$(abs_builddir)/../i/MPC" -DGMP_INCLUDE_DIR="@GMP_INCLUDEDIR@" -DGMP_LIBRARY="@GMP_LIBDIR@/libgmp.so" -DMPFR_INCLUDE_DIR="@MPFR_INCLUDEDIR@" -DMPFR_LIBRARY="@MPFR_LIBDIR@/libmpfr.so" -DMPC_INCLUDE_DIR="@MPC_INCLUDEDIR@" -DMPC_LIBRARY="@MPC_LIBDIR@/libmpc.so"

all: $(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done
$(GMPFRXX_MKII_WORK)/.verify_done:
	@set -e; cd $(abs_srcdir); expected=`awk '$$2 == "gmpfrxx_mkII/download/$(GMPFRXX_MKII_DISTFILE)" {print $$1}' ../distfiles.sha256`; test -n "$$expected"; actual=`sha256sum download/$(GMPFRXX_MKII_DISTFILE) | awk '{print $$1}'`; test "$$expected" = "$$actual"; mkdir -p $(GMPFRXX_MKII_WORK); touch $@
$(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done: $(GMPFRXX_MKII_WORK)/.verify_done
	rm -rf $(GMPFRXX_MKII_INTERNAL_PREFIX)
	mkdir -p $(GMPFRXX_MKII_INTERNAL_PREFIX)
	tar -xJf $(abs_srcdir)/download/$(GMPFRXX_MKII_DISTFILE) -C $(GMPFRXX_MKII_INTERNAL_PREFIX) --strip-components=1
	touch $@
check-local: $(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done
	cmake -S $(GMPFRXX_MKII_INTERNAL_PREFIX) -B $(GMPFRXX_MKII_CHECK_BUILD) $(GMPFRXX_MKII_OPTIONS) $(GMPFRXX_MKII_DEPS)
	cmake --build $(GMPFRXX_MKII_CHECK_BUILD) -j$(GMPFRXX_MKII_BUILD_JOBS)
	ctest --test-dir $(GMPFRXX_MKII_CHECK_BUILD) --output-on-failure
clean-local distclean-local:
	rm -rf $(GMPFRXX_MKII_WORK) $(GMPFRXX_MKII_INTERNAL_PREFIX)
install-data-local: $(GMPFRXX_MKII_INTERNAL_PREFIX)/.extract_done
	$(MKDIR_P) $(DESTDIR)$(includedir)/mplapack
	cp -R $(GMPFRXX_MKII_INTERNAL_PREFIX)/include/. $(DESTDIR)$(includedir)/mplapack/
EXTRA_DIST = download/$(GMPFRXX_MKII_DISTFILE)
