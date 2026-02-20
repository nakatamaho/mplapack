SUBDIRS = compare matgen lin eig

CHECK_SUBDIRS = compare lin eig

.PHONY: check-recursive \
        check-subdir-compare check-subdir-lin check-subdir-eig

check-recursive: check-subdir-compare check-subdir-lin check-subdir-eig

check-subdir-compare:
	$(AM_V_at)$(MAKE) -C compare $(AM_MAKEFLAGS) check

check-subdir-lin:
	$(AM_V_at)$(MAKE) -C lin $(AM_MAKEFLAGS) check

check-subdir-eig:
	$(AM_V_at)$(MAKE) -C eig $(AM_MAKEFLAGS) check

check-local:
	python summarize_mplapack_tests.py .
