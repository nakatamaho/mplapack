MAKEFLAGS = -j 1

SUBDIRS = @GMP_DIR@ @MPFR_DIR@ @MPC_DIR@ @QD_DIR@ @LAPACK_DIR@ @OPENBLAS_DIR@
if IS_MINGW
SUBDIRS += dlfcn-win32
endif

all:

clean:
	rm -rf */work
distclean:

