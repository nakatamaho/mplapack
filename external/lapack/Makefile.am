LAPACK_VERSION    := $(or $(LAPACK_VERSION),3.12.1)
LAPACK_GIT_URL    := $(or $(LAPACK_GIT_URL),https://github.com/Reference-LAPACK/lapack.git)
LAPACK_GIT_REF    := $(or $(LAPACK_GIT_REF),master)
LAPACK_DISTFILE   := lapack/download/lapack-$(LAPACK_VERSION).tar.gz
MAKEFLAGS =

EXTRACT_STAMP=.lapack_extract_$(LAPACK_VERSION)_done
PATCH_STAMP=.lapack_patch_$(LAPACK_VERSION)_done
CONFIGURE_STAMP=.lapack_configure_$(LAPACK_VERSION)_done
BUILD_STAMP=.lapack_build_$(LAPACK_VERSION)_done
INSTALL_INTERNAL_STAMP=.lapack_install_internal_$(LAPACK_VERSION)_done
CHECK_STAMP=.lapack_check_$(LAPACK_VERSION)_done

all: extract config build check install_internal

LAPACK_SOURCE_DIR = lapack-$(LAPACK_VERSION)

.PHONY: verify-distfile
verify-distfile:
	@if test "x$(LAPACK_VERSION)" = "xcurrent"; then \
		echo "Skipping SHA256 check for LAPACK current git checkout"; \
	else \
		distfile="$(LAPACK_DISTFILE)"; \
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
	fi

work/$(EXTRACT_STAMP):
	@$(MAKE) verify-distfile LAPACK_VERSION="$(LAPACK_VERSION)" LAPACK_DISTFILE="$(LAPACK_DISTFILE)"
	rm -rf work ; mkdir -p work/internal
	@if test "x$(LAPACK_VERSION)" = "xcurrent"; then \
		git clone "$(LAPACK_GIT_URL)" "$(abs_builddir)/work/internal/$(LAPACK_SOURCE_DIR)" && \
		cd "$(abs_builddir)/work/internal/$(LAPACK_SOURCE_DIR)" && \
		git checkout "$(LAPACK_GIT_REF)" && \
		git rev-parse HEAD > .lapack_git_head; \
	else \
		cd "$(abs_builddir)/work/internal" && \
		tar xvfz "$(abs_srcdir)/download/lapack-$(LAPACK_VERSION).tar.gz"; \
	fi
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

LAPACK_WORK_DIR = $(abs_builddir)/work/internal/$(LAPACK_SOURCE_DIR)
LAPACK_PATCH_DIR = $(abs_srcdir)/patches/$(LAPACK_VERSION)
LAPACK_VERSION_KEY = $(subst .,_,$(LAPACK_VERSION))

PATCHES_current =

PATCHES_3_9_1 = \
	patch-zchkhe_aa_2stage.f \
	patch-ddrvst2stg.f \
	patch-ddrvst.f \
	patch-dlatm4.f \
	patch-zlqt04.f \
	patch-zlqt05.f \
	patch-ztsqr01.f \
	patch-iseed

PATCHES_3_12_1 = \
	patch-iseed \
	patch-xLAQP2RK.f \
	patch-xxxGVD.f \
	patch-stebz-error-propagation \
	patch-redundant-dble-real \
	patch-trsyl3-ldswork-overwrite \
	patch-xpgvx-info-to-iinfo \
	patch-zunbdb3-conjugate-x21 \
	patch-fix-testparameters \
	patch-xERRTSQR.f \
	patch-tsqr01-cq-residual \
	patch-cabs1_unification \
	patch-cabsmax

LAPACK_PATCHES = $(PATCHES_$(LAPACK_VERSION_KEY))

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	@if test -f "$(LAPACK_PATCH_DIR)/patch-print"; then \
		cd "$(LAPACK_WORK_DIR)" && patch -p1 < "$(LAPACK_PATCH_DIR)/patch-print" && \
		cp "$(abs_srcdir)/patches/print.c" "$(LAPACK_WORK_DIR)/BLAS/SRC"; \
	fi
if IS_MINGW
	@if test -f "$(LAPACK_PATCH_DIR)/patch-mingw"; then \
		cd "$(LAPACK_WORK_DIR)" && patch -p1 < "$(LAPACK_PATCH_DIR)/patch-mingw"; \
	fi
endif
	@if test -n "$(LAPACK_PATCHES)"; then \
		$(foreach patch,$(LAPACK_PATCHES), \
			cd "$(LAPACK_WORK_DIR)" && patch -p1 < "$(LAPACK_PATCH_DIR)/$(patch)" && ) \
		:; \
	fi
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	echo "SHELL = /bin/sh"      >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "CC = $(CC)"           >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "CFLAGS = $(CFLAGS)"   >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "FC = $(FC)"           >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
if IS_INTELCC
	echo "FFLAGS = $(FCFLAGS) -recursive -fp-model strict -assume protect-parens" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
else
	echo "FFLAGS = $(FCFLAGS)" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
endif
	echo "FFLAGS_DRV = \$$(FFLAGS)" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "FFLAGS_NOOPT = -O0" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
if IS_MINGW
	echo "LDFLAGS = -Wl,--stack,33554432" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
else
	echo "LDFLAGS =" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
endif
	echo "AR = $(AR)"        >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "ARFLAGS = cr"      >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "RANLIB = $(RANLIB)"   >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
if IS_INTELCC
	echo "TIMER = INT_CPU_TIME"      >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
else
	echo "TIMER = INT_ETIME"      >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
endif
	echo "TOPSRCDIR = $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "BLASLIB      = \$$(TOPSRCDIR)/libblas.a"     >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "LAPACKLIB    = \$$(TOPSRCDIR)/liblapack.a"   >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "LAPACKELIB   = \$$(TOPSRCDIR)/liblapacke.a"  >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "TMGLIB       = \$$(TOPSRCDIR)/libtmglib.a"   >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "LAPACK_WITH_TMG = Yes"                       >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	+cd $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/BLAS/SRC ; ${MAKE}
	+cd $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION) ; ${MAKE}
	+cd $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/LAPACKE ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	rm -rf $(abs_builddir)/../i/LAPACK
	mkdir -p $(abs_builddir)/../i/LAPACK/lib
	mkdir -p $(abs_builddir)/../i/LAPACK/include
	cp $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/liblapack.*  $(abs_builddir)/../i/LAPACK/lib 
	cp $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/libblas.*    $(abs_builddir)/../i/LAPACK/lib
	cp $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/libtmglib.*  $(abs_builddir)/../i/LAPACK/lib
	cp $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/liblapacke.* $(abs_builddir)/../i/LAPACK/lib
	cp $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/LAPACKE/include/*.h $(abs_builddir)/../i/LAPACK/include
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

distdir:

clean-local:
	rm -rf work
distclean-local:
	rm -rf work
