SUBDIRS = @GMP_DIR@ @MPFR_DIR@ @MPC_DIR@ @QD_DIR@ @LAPACK_DIR@ @OPENBLAS_DIR@
if IS_MINGW
SUBDIRS += dlfcn-win32
endif

EXTRA_DIST = distfiles.sha256

all:

.PHONY: verify-distfiles
verify-distfiles:
	@cd "$(abs_srcdir)" && \
	if command -v sha256sum >/dev/null 2>&1; then \
		sha256sum -c distfiles.sha256; \
	elif command -v shasum >/dev/null 2>&1; then \
		shasum -a 256 -c distfiles.sha256; \
	else \
		echo "error: no SHA256 checksum tool found (sha256sum or shasum)" >&2; \
		exit 1; \
	fi

clean-local:
	rm -rf */work

distclean-local:
	rm -rf i

uninstall-hook:
	-rm -f $(DESTDIR)$(infodir)/dir
