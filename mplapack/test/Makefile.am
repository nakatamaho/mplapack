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

# -----------------------------------------------------------------------------
# Output directory: CPU_COMPILER_OS  (or CPU_COMPILER_mingw for MinGW builds)
#
# cpu_os_detection.sh prints "CPU|OS" (e.g. "Ryzen_Threadripper_3970X|ubuntu22_04").
# MPLAPACK_CC_TAG is set by AX_MPLAPACK_CC_TAG in configure.ac (e.g. "gcc_14_2_0").
#
# Native layout : Ryzen_Threadripper_3970X_gcc_14_2_0_ubuntu22_04/<test>.<profile>.out
# MinGW  layout : Ryzen_Threadripper_3970X_gcc_14_2_0_windows/<test>.<profile>.out
# -----------------------------------------------------------------------------

_CPU_OS_RAW := $(shell bash $(abs_top_srcdir)/misc/cpu_os_detection.sh)
_CPU_TAG    := $(word 1,$(subst |, ,$(_CPU_OS_RAW)))
_OS_TAG     := $(word 2,$(subst |, ,$(_CPU_OS_RAW)))

if IS_MINGW
OUTDIR = $(_CPU_TAG)_windows_$(MPLAPACK_CC_TAG)
else
OUTDIR = $(_CPU_TAG)_$(_OS_TAG)_$(MPLAPACK_CC_TAG)
endif

check-local:
	python3 $(top_srcdir)/misc/summarize_mplapack_tests.py $(OUTDIR) eig
	python3 $(top_srcdir)/misc/summarize_mplapack_tests.py $(OUTDIR) lin
	@if test -n "$${MPLAPACK_TEST_RESULTS_STAGING:-}"; then \
		for suite in eig lin; do \
			src="$$suite/results/$(OUTDIR)"; \
			dst="$${MPLAPACK_TEST_RESULTS_STAGING}/$$suite/results/$(OUTDIR)"; \
			if test -d "$$src"; then \
				mkdir -p "$$dst"; \
				if command -v rsync >/dev/null 2>&1; then \
					rsync -a "$$src/" "$$dst/"; \
				else \
					cp -a "$$src/." "$$dst/"; \
				fi; \
				echo "Staged $$suite test results: $$src -> $$dst"; \
			fi; \
		done; \
	fi
