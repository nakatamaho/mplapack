GMPFRXX_MKII_VERSION = 1.1.0
GMPFRXX_MKII_DISTFILE = gmpfrxx_mkII.$(GMPFRXX_MKII_VERSION).tar.xz
GMPFRXX_MKII_SOURCE_DIR = gmpfrxx_mkII.$(GMPFRXX_MKII_VERSION)
GMPFRXX_MKII_INTERNAL_PREFIX = $(abs_builddir)/../i/GMPFRXX_MKII
GMPFRXX_MKII_FINAL_PREFIX = $(prefix)
GMPFRXX_MKII_WORK = $(abs_builddir)/work
GMPFRXX_MKII_INTERNAL_BUILD = $(GMPFRXX_MKII_WORK)/internal
GMPFRXX_MKII_FINAL_BUILD = $(GMPFRXX_MKII_WORK)/install
GMPFRXX_MKII_BUILD_JOBS = 32
GMPFRXX_MKII_OPTIONS = -DGMPFRXX_MKII_DEPS_AUTO_FETCH=OFF -DGMPFRXX_MKII_BUILD_EXAMPLES=OFF -DGMPFRXX_MKII_BUILD_BENCHMARKS=OFF -DGMPFRXX_MKII_COMPONENTS=GMP,MPFR,MPC
GMPFRXX_MKII_DEPS = -DCMAKE_PREFIX_PATH="$(abs_builddir)/../i/GMP;$(abs_builddir)/../i/MPFR;$(abs_builddir)/../i/MPC" -DGMP_INCLUDE_DIR="@GMP_INCLUDEDIR@" -DGMP_LIBRARY="@GMP_LIBDIR@/libgmp.so" -DMPFR_INCLUDE_DIR="@MPFR_INCLUDEDIR@" -DMPFR_LIBRARY="@MPFR_LIBDIR@/libmpfr.so" -DMPC_INCLUDE_DIR="@MPC_INCLUDEDIR@" -DMPC_LIBRARY="@MPC_LIBDIR@/libmpc.so"

all: $(GMPFRXX_MKII_WORK)/.install_internal_done
$(GMPFRXX_MKII_WORK)/.verify_done:
	@set -e; cd $(abs_srcdir); expected=`awk '$$2 == "gmpfrxx_mkII/download/$(GMPFRXX_MKII_DISTFILE)" {print $$1}' ../distfiles.sha256`; test -n "$$expected"; actual=`sha256sum download/$(GMPFRXX_MKII_DISTFILE) | awk '{print $$1}'`; test "$$expected" = "$$actual"; mkdir -p $(GMPFRXX_MKII_WORK); touch $@
$(GMPFRXX_MKII_WORK)/.extract_done: $(GMPFRXX_MKII_WORK)/.verify_done
	rm -rf $(GMPFRXX_MKII_WORK)/$(GMPFRXX_MKII_SOURCE_DIR); mkdir -p $(GMPFRXX_MKII_WORK); tar -xJf $(abs_srcdir)/download/$(GMPFRXX_MKII_DISTFILE) -C $(GMPFRXX_MKII_WORK); touch $@
$(GMPFRXX_MKII_WORK)/.configure_internal_done: $(GMPFRXX_MKII_WORK)/.extract_done
	cmake -S $(GMPFRXX_MKII_WORK)/$(GMPFRXX_MKII_SOURCE_DIR) -B $(GMPFRXX_MKII_INTERNAL_BUILD) -DCMAKE_INSTALL_PREFIX=$(GMPFRXX_MKII_INTERNAL_PREFIX) $(GMPFRXX_MKII_OPTIONS) $(GMPFRXX_MKII_DEPS) && touch $@
$(GMPFRXX_MKII_WORK)/.build_internal_done: $(GMPFRXX_MKII_WORK)/.configure_internal_done
	cmake --build $(GMPFRXX_MKII_INTERNAL_BUILD) -j$(GMPFRXX_MKII_BUILD_JOBS) && touch $@
$(GMPFRXX_MKII_WORK)/.install_internal_done: $(GMPFRXX_MKII_WORK)/.build_internal_done
	cmake --install $(GMPFRXX_MKII_INTERNAL_BUILD) && touch $@
install-exec-local:
	@if test "$(GMPFRXX_MKII_SOURCE)" = bundled; then $(MAKE) $(GMPFRXX_MKII_WORK)/.configure_final_done; $(MAKE) $(GMPFRXX_MKII_WORK)/.build_final_done; cmake --install $(GMPFRXX_MKII_FINAL_BUILD); fi
$(GMPFRXX_MKII_WORK)/.configure_final_done: $(GMPFRXX_MKII_WORK)/.extract_done
	cmake -S $(GMPFRXX_MKII_WORK)/$(GMPFRXX_MKII_SOURCE_DIR) -B $(GMPFRXX_MKII_FINAL_BUILD) -DCMAKE_INSTALL_PREFIX=$(GMPFRXX_MKII_FINAL_PREFIX) $(GMPFRXX_MKII_OPTIONS) $(GMPFRXX_MKII_DEPS) && touch $@
$(GMPFRXX_MKII_WORK)/.build_final_done: $(GMPFRXX_MKII_WORK)/.configure_final_done
	cmake --build $(GMPFRXX_MKII_FINAL_BUILD) -j$(GMPFRXX_MKII_BUILD_JOBS) && touch $@
clean-local distclean-local:
	rm -rf $(GMPFRXX_MKII_WORK)
EXTRA_DIST = download/$(GMPFRXX_MKII_DISTFILE)
