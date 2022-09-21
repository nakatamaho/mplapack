LAPACKVERSION=$(LAPACK_VERSION)
MAKEFLAGS =

EXTRACT_STAMP=.lapack_extract_done
PATCH_STAMP=.lapack_patch_done
CONFIGURE_STAMP=.lapack_configure_done
BUILD_STAMP=.lapack_build_done
INSTALL_INTERNAL_STAMP=.lapack_install_internal_done
CHECK_STAMP=.lapack_check_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal
	cd $(abs_builddir)/work/internal; tar xvfz ../../download/lapack-$(LAPACKVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
if IS_MINGW
	cd $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION); patch -p0 < $(abs_builddir)/patches/patch-mingw
endif
if IS_MIC
	cd $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION); patch -p0 < $(abs_builddir)/patches/patch-mic
endif
	cd $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION); patch -p0 < $(abs_builddir)/patches/patch-print
	cp $(abs_builddir)/patches/print.c $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/BLAS/SRC
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	echo "SHELL = /bin/sh"      >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "CC = $(CC)"           >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "CFLAGS = $(CFLAGS)"   >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "FC = $(FC)"           >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
if IS_MINGW
	echo "FFLAGS = -O2 -march=native" >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
else
	echo "FFLAGS = $(FCFLAGS)" >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
endif
	echo "FFLAGS_DRV = \$$(FFLAGS)" >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "FFLAGS_NOOPT = -O0" >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "LDFLAGS ="       >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "AR = $(AR)"         >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "ARFLAGS = cr"      >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "RANLIB = $(RANLIB)"   >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
if IS_INTELCC
	echo "TIMER = INT_CPU_TIME"      >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
else
	echo "TIMER = INT_ETIME"      >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
endif
	echo "TOPSRCDIR = $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)" >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "BLASLIB      = \$$(TOPSRCDIR)/libblas.a"     >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "LAPACKLIB    = \$$(TOPSRCDIR)/liblapack.a"   >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "LAPACKELIB   = \$$(TOPSRCDIR)/liblapacke.a"  >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "TMGLIB       = \$$(TOPSRCDIR)/libtmglib.a"   >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	echo "LAPACK_WITH_TMG = Yes"                       >> $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/make.inc
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/BLAS/SRC ; ${MAKE}
	cd $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION) ; ${MAKE}
	cd $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/LAPACKE ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	rm -rf $(abs_builddir)/../i/LAPACK
	mkdir -p $(abs_builddir)/../i/LAPACK/lib
	mkdir -p $(abs_builddir)/../i/LAPACK/include
	cp $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/liblapack.*  $(abs_builddir)/../i/LAPACK/lib 
	cp $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/libblas.*    $(abs_builddir)/../i/LAPACK/lib
	cp $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/libtmglib.*  $(abs_builddir)/../i/LAPACK/lib
	cp $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/liblapacke.* $(abs_builddir)/../i/LAPACK/lib
	cp $(abs_builddir)/work/internal/lapack-$(LAPACKVERSION)/LAPACKE/include/*.h $(abs_builddir)/../i/LAPACK/include
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

distdir:

clean:
	rm -rf work
distclean:
	rm -rf work

