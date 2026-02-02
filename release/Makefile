# Makefile
# MPLAPACK Release Engineering

SCRIPT_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
PROJECT_ROOT := $(abspath $(SCRIPT_DIR)/..)
CONF_FILE ?= $(SCRIPT_DIR)/build-matrix.conf
LOGDIR ?= $(SCRIPT_DIR)/logs/$(shell date +%Y%m%d_%H%M%S)
HOST_WORK_DIR ?= $(PROJECT_ROOT)
HOST_CCACHE_DIR ?= /work/ccache

export CONF_FILE LOGDIR PROJECT_ROOT HOST_WORK_DIR HOST_CCACHE_DIR

.PHONY: all branch tarball dist clean list help
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
	@FILTER_NAME=intel ./build-all.sh

mingw:
	@FILTER_NAME=mingw ./build-all.sh

cuda:
	@FILTER_NAME=cuda ./build-all.sh

ubuntu:
	@FILTER_NAME=ubuntu ./build-all.sh

debian:
	@FILTER_NAME=debian ./build-all.sh

fedora:
	@FILTER_NAME=fedora ./build-all.sh

rocky:
	@FILTER_NAME=rocky ./build-all.sh

suse:
	@FILTER_NAME=suse ./build-all.sh

alpine:
	@FILTER_NAME=alpine ./build-all.sh

# === Arch Filters ===

amd64:
	@FILTER_ARCH=amd64 ./build-all.sh

arm64:
	@FILTER_ARCH=arm64 ./build-all.sh

i386:
	@FILTER_ARCH=386 ./build-all.sh

ppc64le:
	@FILTER_ARCH=ppc64le ./build-all.sh

s390x:
	@FILTER_ARCH=s390x ./build-all.sh

mips64le:
	@FILTER_ARCH=mips64le ./build-all.sh

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
	@echo ""
	@echo "Architecture notes:"
	@echo "  - i386/s390x/mips64le: Debian only (QEMU, very slow)"
	@echo "  - arm64/ppc64le: QEMU on non-native hosts"
	@echo "  - amd64: Native on x86_64 hosts"

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
	rm -rf logs/
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
	@echo "  make clean       - remove logs and images"
	@echo ""
	@echo "Architecture notes:"
	@echo "  - i386/s390x/mips64le: Debian only (QEMU, very slow)"
	@echo "  - arm64/ppc64le: QEMU on non-native hosts"
	@echo "  - amd64: Native on x86_64 hosts"
