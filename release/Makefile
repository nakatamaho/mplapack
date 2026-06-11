# Makefile
# MPLAPACK Release Engineering

SCRIPT_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
PROJECT_ROOT := $(abspath $(SCRIPT_DIR)/..)
CONF_FILE ?= $(SCRIPT_DIR)/build-matrix.conf
LOGDIR ?= $(SCRIPT_DIR)/logs/$(shell date +%Y%m%d_%H%M%S)
HOST_WORK_DIR ?= /work/mplapack-work
HOST_CCACHE_DIR ?= /work/ccache
SUCCESS_DIR ?= $(SCRIPT_DIR)/success
MACOS_REMOTE_SSH ?= ssh
REMOTE_LINUX_SSH ?= ssh

export CONF_FILE LOGDIR PROJECT_ROOT HOST_WORK_DIR HOST_CCACHE_DIR SUCCESS_DIR

TIER1_TARGETS := tier1-macos-arm64 tier1-macos-amd64 tier1-ubuntu2404-arm64 tier1-ubuntu2604-arm64 tier1-ubuntu2404-amd64 tier1-ubuntu2604-amd64 tier1-ubuntu2404-mingw64-amd64 tier1-ubuntu2604-mingw64-amd64 tier1-ubuntu2404-inteloneapi-amd64 tier1-ubuntu2604-inteloneapi-amd64 tier1-debian12-i386 tier1-debian13-i386
TIER1_MACOS_TARGETS := tier1-macos-arm64 tier1-macos-amd64
TIER1_LINUX_TARGETS := tier1-ubuntu2404-arm64 tier1-ubuntu2604-arm64 tier1-ubuntu2404-amd64 tier1-ubuntu2604-amd64 tier1-ubuntu2404-mingw64-amd64 tier1-ubuntu2604-mingw64-amd64 tier1-ubuntu2404-inteloneapi-amd64 tier1-ubuntu2604-inteloneapi-amd64 tier1-debian12-i386 tier1-debian13-i386
TIER2_TARGETS :=
RUN_REMOTE_TARGETS = CONF_FILE=$(CONF_FILE) LOGDIR=$(LOGDIR) PROJECT_ROOT=$(PROJECT_ROOT) SUCCESS_DIR=$(SUCCESS_DIR) MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-targets-by-host.sh

.PHONY: all branch tarball dist clean clean-remote clean-all clean-success list help
.PHONY: intel mingw cuda ubuntu debian fedora rocky suse alpine
.PHONY: amd64 arm64 ppc64le i386 s390x mips64le
.PHONY: branch-ubuntu branch-debian branch-amd64 branch-cuda branch-rocky branch-suse
.PHONY: tarball-ubuntu tarball-amd64
.PHONY: tier1 tier1-macos tier1-linux tier2 tier1-macos-arm64 tier1-macos-amd64 tier1-ubuntu2404-arm64 tier1-ubuntu2604-arm64 tier1-ubuntu2404-amd64 tier1-ubuntu2604-amd64 tier1-ubuntu2404-mingw64-amd64 tier1-ubuntu2604-mingw64-amd64 tier1-ubuntu2404-inteloneapi-amd64 tier1-ubuntu2604-inteloneapi-amd64 tier1-debian12-i386 tier1-debian13-i386
.PHONY: ccache-stats ccache-clear

# === Main Targets ===

all:
	@./build-all.sh

branch:
	@./build-all.sh branch

dist:
	@./build-all.sh dist

tarball:
	@TARBALL= ./build-all.sh tarball

tarball-with:
	@test -n "$(TARBALL)" || (echo "Usage: make tarball-with TARBALL=path/to/mplapack-x.y.z.tar.xz" && exit 1)
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

tier1:
	@$(RUN_REMOTE_TARGETS) $(TIER1_TARGETS)

tier2:
	@./build-all.sh branch

tier1-macos:
	@$(RUN_REMOTE_TARGETS) $(TIER1_MACOS_TARGETS)

tier1-linux:
	@$(RUN_REMOTE_TARGETS) $(TIER1_LINUX_TARGETS)

tier1-macos-arm64:
	@MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) ./run-remote-macos.sh tier1-macos-arm64

tier1-macos-amd64:
	@MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) ./run-remote-macos.sh tier1-macos-amd64

tier1-ubuntu2404-arm64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2404-arm64

tier1-ubuntu2604-arm64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2604-arm64

tier1-ubuntu2404-amd64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2404-amd64

tier1-ubuntu2604-amd64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2604-amd64

tier1-ubuntu2404-mingw64-amd64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2404-mingw64-amd64

tier1-ubuntu2604-mingw64-amd64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2604-mingw64-amd64

tier1-ubuntu2404-inteloneapi-amd64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2404-inteloneapi-amd64

tier1-ubuntu2604-inteloneapi-amd64:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-ubuntu2604-inteloneapi-amd64

tier1-debian12-i386:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-debian12-i386

tier1-debian13-i386:
	@REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier1-debian13-i386

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
	@TARBALL= FILTER_NAME=ubuntu ./build-all.sh tarball

tarball-amd64:
	@TARBALL= FILTER_ARCH=amd64 ./build-all.sh tarball

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
	@printf "%-38s %-45s %-72s %s\n" "NAME" "BASE" "ARCHS" "SOURCE"
	@echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
	@grep -v '^#' $(CONF_FILE) | grep -v '^$$' | \
		awk -F'|' '{src=($$6 ~ /^remote-/ ? $$6 : $$5); printf "%-38s %-45s %-72s %s\n", $$1, $$2, $$3, src}'

clean: clean-remote
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


clean-remote:
	@CONF_FILE=$(CONF_FILE) MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./clean-remote-workdirs.sh

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
	@echo "  make branch      - branch matrix builds only"
	@echo "  make dist        - make dist only"
	@echo "  make tarball     - run make dist, then test that tarball"
	@echo "  make tarball-with TARBALL=path/to/foo.tar.xz"
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
	@echo "  make tier1              - run all Tier 1 remote tests"
	@echo "  make tier1-macos        - run Tier 1 macOS remote tests"
	@echo "  make tier1-linux        - run Tier 1 Linux remote tests"
	@echo "  make tier1-macos-amd64  - run Tier 1 macOS amd64 test via SSH"
	@echo "  make tier1-macos-arm64  - run Tier 1 macOS arm64 test via SSH"
	@echo "  make tier1-ubuntu2404-arm64 - run Tier 1 Ubuntu 24.04 arm64 test via SSH/Docker"
	@echo "  make tier1-ubuntu2404-inteloneapi-amd64 - run Tier 1 Ubuntu 24.04 Intel oneAPI amd64 test"
	@echo "  make tier1-ubuntu2604-inteloneapi-amd64 - run Tier 1 Ubuntu 26.04 Intel oneAPI amd64 test"
	@echo "  make tier1-debian12-i386 - run Tier 1 Debian 12 i386 test"
	@echo "  make tier1-debian13-i386 - run Tier 1 Debian 13 i386 test"
	@echo "  make tier2              - run Tier 2 Docker build matrix"
	@echo "  make arm64       - aarch64"
	@echo "  make ppc64le     - PowerPC 64-bit LE"
	@echo "  make i386        - x86 32-bit"
	@echo "  make s390x       - IBM Z mainframe"
	@echo "  make mips64le    - MIPS 64-bit LE"
	@echo ""
	@echo "Combined filters:"
	@echo "  make branch-ubuntu    - Ubuntu branch builds"
	@echo "  make branch-debian    - Debian branch builds"
	@echo "  make branch-rocky     - Rocky Linux branch builds"
	@echo "  make branch-suse      - openSUSE branch builds"
	@echo "  make branch-cuda      - CUDA branch builds"
	@echo "  make tarball-amd64    - amd64 tarball tests"
	@echo ""
	@echo "ccache:"
	@echo "  make ccache-stats     - show ccache statistics"
	@echo "  make ccache-clear     - clear ccache"
	@echo "  HOST_CCACHE_DIR=/path - override ccache directory (default: /work/ccache)"
	@echo ""
	@echo "Mount options:"
	@echo "  HOST_WORK_DIR=/path   - override source directory (default: project root)"
	@echo "  MACOS_REMOTE_SSH=ssh - override remote macOS SSH command"
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
