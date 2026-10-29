# mpc_cuda (https://github.com/tkouya/mpc_cuda): GMP/MPFR/MPC for CUDA kernels.
# Used by libmplapack_mpfr_opt_cuda (--enable-mpfr-cuda) when no installed
# mpc_cuda is given with --with-mpc-cuda-includedir.
#
# The distfile is `git archive` of upstream commit abd9407 with a generated
# configure.  mpc_cuda builds from the GMP/MPFR/MPC sources, taken from the
# distfiles of ../gmp, ../mpfr and ../mpc.  It needs CUDA 13 (older nvcc
# rejects the early-clobber asm constraints) and builds for one GPU
# architecture, sm_$(MPFR_CUDA_ARCH), to which MPC_CUDA_EXTRA_GENCODE adds the
# other architectures of --with-mpfr-cuda-arch (e.g. 80,90: sm_80 and sm_90).
MPC_CUDA_VERSION=0.0.1-abd9407
MPC_CUDA_DISTFILE=mpc_cuda/download/mpc_cuda-$(MPC_CUDA_VERSION).tar.xz
GMP_DISTFILE=gmp/download/gmp-6.3.0.tar.xz
MPFR_DISTFILE=mpfr/download/mpfr-4.2.2.tar.bz2
MPC_DISTFILE=mpc/download/mpc-1.4.1.tar.xz
MAKEFLAGS =

SRCDIR=$(abs_builddir)/work/mpc_cuda-$(MPC_CUDA_VERSION)
PREFIX=$(abs_builddir)/../i/MPC_CUDA

EXTRACT_STAMP=.mpc_cuda_extract_done
CONFIGURE_STAMP=.mpc_cuda_configure_done
BUILD_STAMP=.mpc_cuda_build_done
INSTALL_INTERNAL_STAMP=.mpc_cuda_install_internal_done
PATCHES=$(srcdir)/patches/patch-cu_freal-correctly-rounded-div-sqrt
EXTRA_DIST=$(PATCHES) patches/README.md

all: extract config build install_internal

.PHONY: verify-distfiles
verify-distfiles:
	@cd "$(abs_srcdir)/.." && \
	for distfile in $(MPC_CUDA_DISTFILE) $(GMP_DISTFILE) $(MPFR_DISTFILE) $(MPC_DISTFILE); do \
		entry=`awk -v f="$$distfile" '$$2 == f { print $$0 }' distfiles.sha256`; \
		if test -z "$$entry"; then \
			echo "error: no SHA256 checksum recorded for $$distfile" >&2; \
			exit 1; \
		fi; \
		if command -v sha256sum >/dev/null 2>&1; then \
			printf '%s\n' "$$entry" | sha256sum -c - || exit 1; \
		elif command -v shasum >/dev/null 2>&1; then \
			printf '%s\n' "$$entry" | shasum -a 256 -c - || exit 1; \
		else \
			echo "error: no SHA256 checksum tool found (sha256sum or shasum)" >&2; \
			exit 1; \
		fi; \
	done

work/$(EXTRACT_STAMP): $(PATCHES)
	@$(MAKE) verify-distfiles
	rm -rf work ; mkdir -p work
	cd $(abs_builddir)/work; tar xfJ $(abs_srcdir)/download/mpc_cuda-$(MPC_CUDA_VERSION).tar.xz
	cd $(SRCDIR); tar xfJ $(abs_srcdir)/../$(GMP_DISTFILE)
	cd $(SRCDIR); tar xfj $(abs_srcdir)/../$(MPFR_DISTFILE)
	cd $(SRCDIR); tar xfJ $(abs_srcdir)/../$(MPC_DISTFILE)
	cd $(SRCDIR); patch -p1 --batch < $(abs_srcdir)/patches/patch-cu_freal-correctly-rounded-div-sqrt
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

# glibc 2.42 exposes C23 rsqrt declarations that conflict with CUDA 13.1;
# X/Open keeps the POSIX declarations (as NVCC_GLIBC_FEATURE_FLAGS of
# mpblas/optimized/mpfr/cuda).  mpc_cuda takes them as part of NVCC; it then
# cannot derive the CUDA library directory from NVCC, so it is given.
NVCC_GLIBC_FEATURE_FLAGS = -U_GNU_SOURCE -D_XOPEN_SOURCE=700

work/$(CONFIGURE_STAMP): work/$(EXTRACT_STAMP)
	cd $(SRCDIR) ; ./configure --prefix=$(PREFIX) --with-cuda-arch=sm_$(MPFR_CUDA_ARCH) "NVCC=$(NVCC) $(NVCC_GLIBC_FEATURE_FLAGS) $(MPC_CUDA_EXTRA_GENCODE)"
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

# mpc_cuda must not be built in parallel: its steps share one build directory.
# Its tools/build_lib.sh prunes objects with sort | comm, which disagree in
# some locales (ja_JP.UTF-8: most of MPFR, cu_mpfr_add etc., was dropped), so
# it runs in the C locale.
work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(SRCDIR) ; env MAKEFLAGS= LC_ALL=C CUDA_LIBDIR=$(CUDATOOLKITHOME)/lib64 $(MAKE) -j1
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(BUILD_STAMP)
	cd $(SRCDIR) ; env MAKEFLAGS= $(MAKE) -j1 install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)
config: work/$(CONFIGURE_STAMP)
build: work/$(BUILD_STAMP)
install_internal: work/$(INSTALL_INTERNAL_STAMP)

# libmpc_cuda.a is linked into libmplapack_mpfr_opt_cuda; nothing is installed.
install:

uninstall:

check:

distdir:

clean-local:
	rm -rf work

distclean-local:
	rm -rf work
