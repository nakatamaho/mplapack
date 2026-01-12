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

if ENABLE__FLOAT128
SUBDIRS += _Float128
endif

if ENABLE__FLOAT64X
SUBDIRS += _Float64x
endif

CHECK_SUBDIRS = $(SUBDIRS)

.PHONY: check-recursive $(CHECK_SUBDIRS:%=check-subdir-%)

check-recursive: $(CHECK_SUBDIRS:%=check-subdir-%)

check-subdir-%:
	$(AM_V_at)$(MAKE) $(AM_MAKEFLAGS) -C $* check
