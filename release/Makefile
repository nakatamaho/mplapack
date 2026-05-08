# Makefile
# MPLAPACK Release Engineering

SCRIPT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
PROJECT_ROOT := $(abspath $(SCRIPT_DIR)/..)
CONF_FILE ?= $(SCRIPT_DIR)/build-matrix.conf
LOGDIR ?= $(SCRIPT_DIR)/logs/$(shell date +%Y%m%d_%H%M%S)
HOST_WORK_DIR ?= /work/mplapack-work
HOST_CCACHE_DIR ?= /work/ccache
SUCCESS_DIR ?= $(SCRIPT_DIR)/success

export CONF_FILE LOGDIR PROJECT_ROOT HOST_WORK_DIR HOST_CCACHE_DIR SUCCESS_DIR

.PHONY: all branch tarball dist clean clean-all clean-success list help
.PHONY: intel mingw cuda ubuntu debian fedora rocky suse alpine
.PHONY: amd64 arm64 ppc64le i386 s390x mips64le
.PHONY: branch-ubuntu branch-debian branch-amd64 branch-cuda branch-rocky branch-suse
.PHONY: tarball-ubuntu tarball-debian tarball-amd64 tarball-cuda tarball-rocky tarball-suse
.PHONY: ccache-stats ccache-clear

# === Main Targets ===

all:
	@./build-all.sh

branch:
	@./build-all.sh branch

dist:
	@./build-all.sh dist

tarball:
	@./build-all.sh tarball

tarball-with:
	@test -n "$(TARBALL)" || (echo "Usage: make tarball-with TARBALL=path/to/mplapack-x.y.z.tar.gz" && exit 1)
	@TARBALL=$(TARBALL) ./build-all.sh tarball

# === Distro Filters ===

intel:
	@FILTER_NAME=intel ./build-all.sh branch

mingw:
	@FILTER_NAME=mingw ./build-all.sh branch

cuda:
	@FILTER_NAME=cuda ./build-all.sh branch

ubuntu:
	@FILTER_NAME=ubuntu ./build-all.sh branch

debian:
	@FILTER_NAME=debian ./build-all.sh branch

fedora:
	@FILTER_NAME=fedora ./build-all.sh branch

rocky:
	@FILTER_NAME=rocky ./build-all.sh branch

suse:
	@FILTER_NAME=suse ./build-all.sh branch

alpine:
	@FILTER_NAME=alpine ./build-all.sh branch

# === Arch Filters ===

amd64:
	@FILTER_ARCH=amd64 ./build-all.sh branch

arm64:
	@FILTER_ARCH=arm64 ./build-all.sh branch

i386:
	@FILTER_ARCH=386 ./build-all.sh branch

ppc64le:
	@FILTER_ARCH=ppc64le ./build-all.sh branch

s390x:
	@FILTER_ARCH=s390x ./build-all.sh branch

mips64le:
	@FILTER_ARCH=mips64le ./build-all.sh branch

# === Combined Filters ===

branch-ubuntu:
	@FILTER_NAME=ubuntu ./build-all.sh branch

branch-debian:
	@FILTER_NAME=debian ./build-all.sh branch

branch-amd64:
	@FILTER_ARCH=amd64 ./build-all.sh branch

branch-cuda:
	@FILTER_NAME=cuda ./build-all.sh branch

branch-rocky:
	@FILTER_NAME=rocky ./build-all.sh branch

branch-suse:
	@FILTER_NAME=suse ./build-all.sh branch

tarball-ubuntu:
	@FILTER_NAME=ubuntu ./build-all.sh tarball

tarball-debian:
	@FILTER_NAME=debian ./build-all.sh tarball

tarball-amd64:
	@FILTER_ARCH=amd64 ./build-all.sh tarball

tarball-cuda:
	@FILTER_NAME=cuda ./build-all.sh tarball

tarball-rocky:
	@FILTER_NAME=rocky ./build-all.sh tarball

tarball-suse:
	@FILTER_NAME=suse ./build-all.sh tarball

# === ccache Utilities ===

ccache-stats:
	@echo "ccache directory: $(HOST_CCACHE_DIR)"
	@ccache -d $(HOST_CCACHE_DIR) -s 2>/dev/null || \
		echo "Install ccache on host or check directory permissions"

ccache-clear:
	@echo "Clearing ccache at $(HOST_CCACHE_DIR)..."
	@ccache -d $(HOST_CCACHE_DIR) -C 2>/dev/null || \
		rm -rf $(HOST_CCACHE_DIR)/*
	@echo "Done."

# === Utilities ===

list:
	@echo "=== Build Matrix ==="
	@printf "%-25s %-45s %-40s %s\n" "NAME" "BASE" "ARCHS" "SOURCE"
	@echo "----------------------------------------------------------------------------------------------------------------------"
	@grep -v '^#' $(CONF_FILE) | grep -v '^$$' | \
		awk -F'|' '{printf "%-25s %-45s %-40s %s\n", $$1, $$2, $$3, $$6}'

clean:
	@mkdir -p logs/ $(SUCCESS_DIR)
	@tmp="$$(mktemp)"; \
	trap 'rm -f "$$tmp"' EXIT; \
	find "$(SUCCESS_DIR)" -maxdepth 1 -type l -name '*.ok' | while IFS= read -r link; do \
		target="$$(readlink -f "$$link" 2>/dev/null || true)"; \
		if [ -n "$$target" ] && [ -e "$$target" ]; then \
			dir="$$(dirname "$$target")"; \
			case "$$dir" in \
				"$(abspath $(SCRIPT_DIR)/logs)"/*) printf '%s\n' "$$dir" ;; \
			esac; \
		fi; \
	done | sort -u > "$$tmp"; \
	if [ -s "$$tmp" ]; then \
		echo "Preserving successful log directories:"; \
		sed 's#^#  #g' "$$tmp"; \
		find logs -mindepth 1 -maxdepth 1 -type d | while IFS= read -r d; do \
			rd="$$(readlink -f "$$d" 2>/dev/null || true)"; \
			if grep -Fxq "$$rd" "$$tmp"; then \
				echo "KEEP $$d"; \
			else \
				echo "RM   $$d"; \
				rm -rf "$$d"; \
			fi; \
		done; \
	else \
		echo "No protected success logs; removing logs/"; \
		rm -rf logs/; \
		mkdir -p logs/; \
	fi
	@# Stop and remove any running containers using our labeled images first.
	@# "docker rmi -f" cannot remove images used by *running* containers even with -f.
	-docker ps -q  --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker stop
	-docker ps -aq --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rm -f
	-docker images -q --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rmi -f


clean-success:
	@rm -rf $(SUCCESS_DIR)
	@mkdir -p $(SUCCESS_DIR)

clean-all:
	@rm -rf logs/ $(SUCCESS_DIR)
	@mkdir -p logs/ $(SUCCESS_DIR)
	@# Stop and remove any running containers using our labeled images first.
	@# "docker rmi -f" cannot remove images used by *running* containers even with -f.
	-docker ps -q  --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker stop
	-docker ps -aq --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rm -f
	-docker images -q --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rmi -f

help:
	@echo "MPLAPACK Release Engineering"
	@echo ""
	@echo "Full test cycle:"
	@echo "  make             - run all (branch -> dist -> tarball)"
	@echo ""
	@echo "Individual phases:"
	@echo "  make branch      - branch tests only"
	@echo "  make dist        - make dist only"
	@echo "  make tarball     - tarball tests (runs dist if needed)"
	@echo "  make tarball-with TARBALL=path/to/foo.tar.gz"
	@echo ""
	@echo "Distro filters:"
	@echo "  make ubuntu      - Ubuntu only"
	@echo "  make debian      - Debian only"
	@echo "  make fedora      - Fedora only"
	@echo "  make rocky       - Rocky Linux only"
	@echo "  make suse        - openSUSE only"
	@echo "  make alpine      - Alpine only"
	@echo "  make intel       - Intel compiler only"
	@echo "  make mingw       - MinGW cross-compile only"
	@echo "  make cuda        - CUDA builds only"
	@echo ""
	@echo "Arch filters:"
	@echo "  make amd64       - x86_64"
	@echo "  make arm64       - aarch64"
	@echo "  make ppc64le     - PowerPC 64-bit LE"
	@echo "  make i386        - x86 32-bit"
	@echo "  make s390x       - IBM Z mainframe"
	@echo "  make mips64le    - MIPS 64-bit LE"
	@echo ""
	@echo "Combined filters:"
	@echo "  make branch-ubuntu    - Ubuntu branch tests"
	@echo "  make branch-debian    - Debian branch tests"
	@echo "  make branch-rocky     - Rocky Linux branch tests"
	@echo "  make branch-suse      - openSUSE branch tests"
	@echo "  make branch-cuda      - CUDA branch tests"
	@echo "  make tarball-amd64    - amd64 tarball tests"
	@echo ""
	@echo "ccache:"
	@echo "  make ccache-stats     - show ccache statistics"
	@echo "  make ccache-clear     - clear ccache"
	@echo "  HOST_CCACHE_DIR=/path - override ccache directory (default: /work/ccache)"
	@echo ""
	@echo "Mount options:"
	@echo "  HOST_WORK_DIR=/path   - override source directory (default: project root)"
	@echo ""
	@echo "Utilities:"
	@echo "  make list        - show build matrix"
	@echo "  make clean       - remove unprotected logs and images"
	@echo "  make clean-success - remove success stamps only"
	@echo "  make clean-all   - remove all logs, success stamps, and images"
	@echo ""
	@echo "Architecture notes:"
	@echo "  - i386/s390x/mips64le: Debian only (QEMU, very slow)"
	@echo "  - arm64/ppc64le: QEMU on non-native hosts"
	@echo "  - amd64: Native on x86_64 hosts"
