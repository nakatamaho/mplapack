SUBDIRS =

if ENABLE_GMP
SUBDIRS += gmp
endif

if ENABLE_MPFR
SUBDIRS += mpfr
endif

if ENABLE_QD
SUBDIRS += qd
endif

if ENABLE_DD
SUBDIRS += dd
endif

if ENABLE_DOUBLE
SUBDIRS += double
endif

if ENABLE_BINARY128
SUBDIRS += binary128
endif

if ENABLE_BINARY80
SUBDIRS += binary80
endif

CHECK_SUBDIRS = $(SUBDIRS)

# The per-backend Makefile.am files are committed (generated once via
# gen.Makefile.am.sh from the *.part templates).  Ship the generator and its
# templates so the test programs can be regenerated from a dist tarball after
# new *.test.cpp cases are added under common/.
EXTRA_DIST = \
	gen.Makefile.am.sh \
	gmp/Makefile.am.part \
	mpfr/Makefile.am.part \
	qd/Makefile.am.part \
	dd/Makefile.am.part \
	double/Makefile.am.part \
	binary80/Makefile.am.part \
	binary128/Makefile.am.part

.PHONY: check-recursive $(CHECK_SUBDIRS:%=check-subdir-%)

check-recursive: $(CHECK_SUBDIRS:%=check-subdir-%)

check-subdir-%:
	$(AM_V_at)$(MAKE) $(AM_MAKEFLAGS) -C $* check
