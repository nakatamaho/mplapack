LAPACK_VERSION ?= 3.12.1
MAKEFLAGS =

EXTRACT_STAMP=.lapack_extract_$(LAPACK_VERSION)_done
PATCH_STAMP=.lapack_patch_$(LAPACK_VERSION)_done
CONFIGURE_STAMP=.lapack_configure_$(LAPACK_VERSION)_done
BUILD_STAMP=.lapack_build_$(LAPACK_VERSION)_done
INSTALL_INTERNAL_STAMP=.lapack_install_internal_$(LAPACK_VERSION)_done
CHECK_STAMP=.lapack_check_$(LAPACK_VERSION)_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal
	cd $(abs_builddir)/work/internal; tar xvfz $(abs_srcdir)/download/lapack-$(LAPACK_VERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

LAPACK_WORK_DIR = $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)
LAPACK_PATCH_DIR = $(abs_srcdir)/patches/$(LAPACK_VERSION)

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
	patch-fix-testparameters \
	patch-xERRTSQR.f \
	patch-cabs1_unification \
	patch-cabsmax

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	cd $(LAPACK_WORK_DIR) && patch -p0 < $(LAPACK_PATCH_DIR)/patch-print
if IS_MINGW
	cd $(LAPACK_WORK_DIR) && patch -p0 < $(LAPACK_PATCH_DIR)/patch-mingw
endif
	@if test "x$(LAPACK_VERSION)" = "x3.9.1"; then \
		$(foreach patch,$(PATCHES_3_9_1), \
			cd $(LAPACK_WORK_DIR) && patch -p0 < $(LAPACK_PATCH_DIR)/$(patch) && ) \
		cp $(abs_srcdir)/patches/print.c $(LAPACK_WORK_DIR)/BLAS/SRC; \
	fi
	@if test "x$(LAPACK_VERSION)" = "x3.12.1"; then \
		$(foreach patch,$(PATCHES_3_12_1), \
			cd $(LAPACK_WORK_DIR) && patch -p1 < $(LAPACK_PATCH_DIR)/$(patch) && ) \
		cp $(abs_srcdir)/patches/print.c $(LAPACK_WORK_DIR)/BLAS/SRC; \
	fi
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	echo "SHELL = /bin/sh"      >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "CC = $(CC)"           >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "CFLAGS = $(CFLAGS)"   >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "FC = $(FC)"           >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
	echo "FFLAGS = $(FCFLAGS)" >> $(abs_builddir)/work/internal/lapack-$(LAPACK_VERSION)/make.inc
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

