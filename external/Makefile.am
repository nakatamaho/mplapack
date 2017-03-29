MAKEFLAGS = -j 1

SUBDIRS =

if NEED_INTERNAL_GMPLIB
SUBDIRS += gmp
endif

if NEED_INTERNAL_MPFRLIB
SUBDIRS += mpfr
endif

if NEED_INTERNAL_MPCLIB
SUBDIRS += mpc
endif

if NEED_INTERNAL_QDLIB
SUBDIRS += qd
endif

# We use BLAS comes from LAPACK.
if ENABLE_DEBUG
SUBDIRS += lapack
endif

all:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i || exit 1; \
	done

clean:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i clean || exit 1; \
	done
	rm -rf i

distclean:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i clean || exit 1; \
	done
	rm -rf i
