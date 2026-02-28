SUBDIRS = @GMP_DIR@ @MPFR_DIR@ @MPC_DIR@ @QD_DIR@ @LAPACK_DIR@ @OPENBLAS_DIR@
if IS_MINGW
SUBDIRS += dlfcn-win32
endif

all:

clean-local:
	rm -rf */work

distclean-local:
	rm -rf i

uninstall-hook:
	-rm -f $(DESTDIR)$(infodir)/dir
