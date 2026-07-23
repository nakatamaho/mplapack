# Makefile
# MPLAPACK Release Engineering

SCRIPT_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
PROJECT_ROOT := $(abspath $(SCRIPT_DIR)/..)
CONF_FILE ?= $(SCRIPT_DIR)/build-matrix.conf
ifeq ($(origin LOGDIR), undefined)
LOGDIR := $(SCRIPT_DIR)/logs/$(shell date +%Y%m%d_%H%M%S)
endif
HOST_WORK_DIR ?= /work/mplapack-work
HOST_CCACHE_DIR ?= /work/ccache
SUCCESS_DIR ?= $(SCRIPT_DIR)/success
FAILED_DIR ?= $(SCRIPT_DIR)/failed
MACOS_REMOTE_SSH ?= ssh
REMOTE_LINUX_SSH ?= ssh
REMOTE_TARBALL_SSH ?= ssh
SOURCE_INFO_FILE ?= $(LOGDIR)/source/source.env

export CONF_FILE LOGDIR PROJECT_ROOT HOST_WORK_DIR HOST_CCACHE_DIR SUCCESS_DIR FAILED_DIR REMOTE_TARBALL_SSH MPLAPACK_REF
export MPLAPACK_SOURCE_MODE MPLAPACK_SOURCE_TARBALL MPLAPACK_SOURCE_METADATA MPLAPACK_SOURCE_PATCH MPLAPACK_SOURCE_STATUS MPLAPACK_SOURCE_LABEL MPLAPACK_SOURCE_TARBALL_SHA256 MPLAPACK_SOURCE_SNAPSHOT_REUSE MPLAPACK_DIST_CACHE_DIR MPLAPACK_DIST_CACHE_REUSE MPLAPACK_HOST_LOCK_DIR MPLAPACK_HOST_LOCK_POLL_SECONDS MPLAPACK_LOCAL_HOST_LOCK

TIER1_TARGETS := tier1-macos-arm64 tier1-macos-amd64 tier1-ubuntu2404-arm64 tier1-ubuntu2604-arm64 tier1-ubuntu2404-amd64 tier1-ubuntu2604-amd64 tier1-ubuntu2404-mingw64-amd64 tier1-ubuntu2604-mingw64-amd64 tier1-ubuntu2404-inteloneapi-amd64 tier1-ubuntu2604-inteloneapi-amd64 tier1-debian12-i386 tier1-debian13-i386
TIER1_MACOS_TARGETS := tier1-macos-arm64 tier1-macos-amd64
TIER1_LINUX_TARGETS := tier1-ubuntu2404-arm64 tier1-ubuntu2604-arm64 tier1-ubuntu2404-amd64 tier1-ubuntu2604-amd64 tier1-ubuntu2404-mingw64-amd64 tier1-ubuntu2604-mingw64-amd64 tier1-ubuntu2404-inteloneapi-amd64 tier1-ubuntu2604-inteloneapi-amd64 tier1-debian12-i386 tier1-debian13-i386
TIER2_TARGETS :=
TIER3_TARGETS := tier3-ubuntu2604-cxxstd-arm64 tier3-ubuntu2604-cxxstd-amd64
LEGACY_QA_IMAGE_TAGS := \
    mplapack-tier1-ubuntu-amd64:latest \
    mplapack-tier1-ubuntu-arm64:latest \
    mplapack-tier1-debian-i386:latest \
    mplapack-tier1-mingw64-amd64:latest \
    mplapack-tier1-ubuntu-inteloneapi-amd64:latest
RUN_REMOTE_TARGETS = CONF_FILE=$(CONF_FILE) LOGDIR=$(LOGDIR) PROJECT_ROOT=$(PROJECT_ROOT) SUCCESS_DIR=$(SUCCESS_DIR) FAILED_DIR=$(FAILED_DIR) MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-targets-by-host.sh
WITH_SOURCE_TARBALL = set -a; . "$(SOURCE_INFO_FILE)"; set +a;

.PHONY: all ref tarball dist prepare-dist-source clean clean-remote clean-remote-realclean clean-all realclean clean-success clean-failed clean-test-results list help
.PHONY: intel mingw cuda ubuntu debian fedora rocky suse alpine
.PHONY: amd64 arm64 ppc64le i386 s390x mips64le
.PHONY: ref-ubuntu ref-debian ref-amd64 ref-cuda ref-rocky ref-suse
.PHONY: tarball-ubuntu tarball-amd64 tarball-arm64
.PHONY: tier1 tier1-macos tier1-linux tier2 tier3 tier1-macos-arm64 tier1-macos-amd64 tier1-ubuntu2404-arm64 tier1-ubuntu2604-arm64 tier1-ubuntu2404-amd64 tier1-ubuntu2604-amd64 tier1-ubuntu2404-mingw64-amd64 tier1-ubuntu2604-mingw64-amd64 tier1-ubuntu2404-inteloneapi-amd64 tier1-ubuntu2604-inteloneapi-amd64 tier1-debian12-i386 tier1-debian13-i386 tier3-ubuntu2604-cxxstd-arm64 tier3-ubuntu2604-cxxstd-amd64
.PHONY: ccache-stats ccache-clear

# === Main Targets ===

all:
	@./build-all.sh

ref:
	@./build-all.sh ref

dist:
	@./build-all.sh dist

prepare-dist-source: $(SOURCE_INFO_FILE)

$(SOURCE_INFO_FILE):
	@PROJECT_ROOT=$(PROJECT_ROOT) LOGDIR=$(LOGDIR) MPLAPACK_SOURCE_INFO_FILE=$(SOURCE_INFO_FILE) ./make-source-snapshot.sh

tier1 tier1-macos tier1-linux $(TIER1_TARGETS) tier2 tier3 $(TIER3_TARGETS) tarball tarball-ubuntu tarball-amd64 tarball-arm64: prepare-dist-source

tarball:
	@$(WITH_SOURCE_TARBALL) TARBALL="$$MPLAPACK_SOURCE_TARBALL" ./build-all.sh tarball

tarball-with:
	@test -n "$(TARBALL)" || (echo "Usage: make tarball-with TARBALL=path/to/mplapack-x.y.z.tar.xz" && exit 1)
	@TARBALL=$(TARBALL) ./build-all.sh tarball

# === Distro Filters ===

intel:
	@FILTER_NAME=intel ./build-all.sh ref

mingw:
	@FILTER_NAME=mingw ./build-all.sh ref

cuda:
	@FILTER_NAME=cuda ./build-all.sh ref

ubuntu:
	@FILTER_NAME=ubuntu ./build-all.sh ref

debian:
	@FILTER_NAME=debian ./build-all.sh ref

fedora:
	@FILTER_NAME=fedora ./build-all.sh ref

rocky:
	@FILTER_NAME=rocky ./build-all.sh ref

suse:
	@FILTER_NAME=suse ./build-all.sh ref

alpine:
	@FILTER_NAME=alpine ./build-all.sh ref

# === Arch Filters ===

amd64:
	@FILTER_ARCH=amd64 ./build-all.sh ref

arm64:
	@FILTER_ARCH=arm64 ./build-all.sh ref

tier1:
	@$(RUN_REMOTE_TARGETS) $(TIER1_TARGETS)

tier2:
	@MPLAPACK_LOG_PREFIX=tier2- ./build-all.sh ref

tier3:
	@MPLAPACK_SOURCE_MODE=$${MPLAPACK_SOURCE_MODE:-ref} $(RUN_REMOTE_TARGETS) $(TIER3_TARGETS)

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

tier3-ubuntu2604-cxxstd-arm64:
	@MPLAPACK_SOURCE_MODE=$${MPLAPACK_SOURCE_MODE:-ref} REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier3-ubuntu2604-cxxstd-arm64

tier3-ubuntu2604-cxxstd-amd64:
	@MPLAPACK_SOURCE_MODE=$${MPLAPACK_SOURCE_MODE:-ref} REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./run-remote-linux-docker.sh tier3-ubuntu2604-cxxstd-amd64

i386:
	@FILTER_ARCH=386 ./build-all.sh ref

ppc64le:
	@FILTER_ARCH=ppc64le ./build-all.sh ref

s390x:
	@FILTER_ARCH=s390x ./build-all.sh ref

mips64le:
	@FILTER_ARCH=mips64le ./build-all.sh ref

# === Combined Filters ===

ref-ubuntu:
	@FILTER_NAME=ubuntu ./build-all.sh ref

ref-debian:
	@FILTER_NAME=debian ./build-all.sh ref

ref-amd64:
	@FILTER_ARCH=amd64 ./build-all.sh ref

ref-cuda:
	@FILTER_NAME=cuda ./build-all.sh ref

ref-rocky:
	@FILTER_NAME=rocky ./build-all.sh ref

ref-suse:
	@FILTER_NAME=suse ./build-all.sh ref

tarball-ubuntu:
	@$(WITH_SOURCE_TARBALL) TARBALL="$$MPLAPACK_SOURCE_TARBALL" FILTER_NAME=ubuntu ./build-all.sh tarball

tarball-amd64:
	@$(WITH_SOURCE_TARBALL) TARBALL="$$MPLAPACK_SOURCE_TARBALL" FILTER_ARCH=amd64 ./build-all.sh tarball

tarball-arm64:
	@$(WITH_SOURCE_TARBALL) TARBALL="$$MPLAPACK_SOURCE_TARBALL" FILTER_ARCH=arm64 ./build-all.sh tarball

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
	@mkdir -p logs/ $(SUCCESS_DIR) $(FAILED_DIR)
	@tmp="$$(mktemp)"; \
	trap 'rm -f "$$tmp"' EXIT; \
	{ find "$(SUCCESS_DIR)" -maxdepth 1 -type l -name '*.ok'; find "$(FAILED_DIR)" -maxdepth 1 -type l -name '*.failed'; } | while IFS= read -r link; do \
		target="$$(readlink -f "$$link" 2>/dev/null || true)"; \
		if [ -n "$$target" ] && [ -e "$$target" ]; then \
			dir="$$(dirname "$$target")"; \
			case "$$dir" in \
				"$(abspath $(SCRIPT_DIR)/logs)"/*) printf '%s\n' "$$dir" ;; \
			esac; \
		fi; \
	done | sort -u > "$$tmp"; \
	if [ -s "$$tmp" ]; then \
		echo "Preserving linked success/failed log directories:"; \
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
		echo "No protected success/failed logs; removing logs/"; \
		rm -rf logs/; \
		mkdir -p logs/; \
	fi
	@# Stop and remove any running containers using our labeled images first.
	@# The image-name fallback catches legacy containers created before labels were added.
	@# Images are removed only when they carry the MPLAPACK QA label.
	@# "docker rmi -f" cannot remove images used by *running* containers even with -f.
	-docker ps -q  --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker stop
	-docker ps -aq --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rm -f
	-@for image in $(LEGACY_QA_IMAGE_TAGS); do docker ps -q  --filter ancestor=$$image; done | sort -u | xargs -r docker stop
	-@for image in $(LEGACY_QA_IMAGE_TAGS); do docker ps -aq --filter ancestor=$$image; done | sort -u | xargs -r docker rm -f
	-docker images -q --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rmi -f


clean-remote:
	@CONF_FILE=$(CONF_FILE) MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./clean-remote-workdirs.sh

clean-remote-realclean:
	@CONF_FILE=$(CONF_FILE) MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./clean-remote-workdirs.sh --realclean

clean-success:
	@rm -rf $(SUCCESS_DIR)
	@mkdir -p $(SUCCESS_DIR)

clean-failed:
	@rm -rf $(FAILED_DIR)
	@mkdir -p $(FAILED_DIR)

clean-test-results:
	@if git -C "$(PROJECT_ROOT)" rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
		git -C "$(PROJECT_ROOT)" clean -fd -- mplapack/test/eig/results mplapack/test/lin/results; \
	else \
		version=$$(awk -F'[][]' '/^m4_define\(\[MPLAPACK_VER_MAJOR\]/{maj=$$4} /^m4_define\(\[MPLAPACK_VER_MINOR\]/{min=$$4} /^m4_define\(\[MPLAPACK_VER_PATCH\]/{patch=$$4} END{if (maj != "" && min != "" && patch != "") print maj "." min "." patch}' "$(PROJECT_ROOT)/configure.ac"); \
		if [ -n "$$version" ]; then \
			for d in "$(PROJECT_ROOT)/mplapack/test/eig/results/$$version" "$(PROJECT_ROOT)/mplapack/test/lin/results/$$version"; do \
				[ -e "$$d" ] || continue; \
				echo "RM   $$d"; \
				rm -rf "$$d"; \
			done; \
		fi; \
	fi

clean-all:
	@rm -rf logs/ $(SUCCESS_DIR) $(FAILED_DIR)
	@mkdir -p logs/ $(SUCCESS_DIR) $(FAILED_DIR)
	@# Stop and remove any running containers using our labeled images first.
	@# The image-name fallback catches legacy containers created before labels were added.
	@# Images are removed only when they carry the MPLAPACK QA label.
	@# "docker rmi -f" cannot remove images used by *running* containers even with -f.
	-docker ps -q  --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker stop
	-docker ps -aq --filter label=org.mplapack.project=mplapack --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rm -f
	-@for image in $(LEGACY_QA_IMAGE_TAGS); do docker ps -q  --filter ancestor=$$image; done | sort -u | xargs -r docker stop
	-@for image in $(LEGACY_QA_IMAGE_TAGS); do docker ps -aq --filter ancestor=$$image; done | sort -u | xargs -r docker rm -f
	-docker images -q --filter label=org.mplapack.purpose=mplapack-qa | xargs -r docker rmi -f

realclean:
	@CONF_FILE=$(CONF_FILE) MACOS_REMOTE_SSH=$(MACOS_REMOTE_SSH) REMOTE_LINUX_SSH=$(REMOTE_LINUX_SSH) ./clean-remote-workdirs.sh --realclean
	@$(MAKE) --no-print-directory clean-all
	@$(MAKE) --no-print-directory clean-test-results
	@rm -rf logs/ $(SUCCESS_DIR) $(FAILED_DIR)
	@for d in "$(HOME)"/mplapack_logs* "$(HOME)"/mplapack_build_log* "$(HOME)"/tmp/tier* mplapack_logs* mplapack_build_log*; do \
		[ -e "$$d" ] || continue; \
		echo "RM   $$d"; \
		rm -rf "$$d"; \
	done

help:
	@echo "MPLAPACK Release Engineering"
	@echo ""
	@echo "Full test cycle:"
	@echo "  make             - run all (ref -> dist -> tarball)"
	@echo ""
	@echo "Individual phases:"
	@echo "  make ref         - Git ref matrix builds only"
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
	@echo "  make ref-ubuntu       - Ubuntu ref builds"
	@echo "  make ref-debian       - Debian ref builds"
	@echo "  make ref-rocky        - Rocky Linux ref builds"
	@echo "  make ref-suse         - openSUSE ref builds"
	@echo "  make ref-cuda         - CUDA ref builds"
	@echo "  make tarball-amd64    - amd64 tarball tests via 172.27.109.80 Docker"
	@echo "  make tarball-arm64    - arm64 tarball tests via 172.27.109.40 Docker"
	@echo ""
	@echo "ccache:"
	@echo "  make ccache-stats     - show ccache statistics"
	@echo "  make ccache-clear     - clear ccache"
	@echo "  HOST_CCACHE_DIR=/path - override ccache directory (default: /work/ccache)"
	@echo ""
	@echo "Mount options:"
	@echo "  HOST_WORK_DIR=/path   - override source directory (default: project root)"
	@echo "  MACOS_REMOTE_SSH=ssh - override remote macOS SSH command"
	@echo "  REMOTE_TARBALL_SSH=ssh - override remote tarball SSH command"
	@echo "  MPLAPACK_LOCAL_HOST_LOCK=host - override local Docker host lock key"
	@echo ""
	@echo "Utilities:"
	@echo "  make list        - show build matrix"
	@echo "  make clean       - remove unprotected logs and images"
	@echo "  make clean-success - remove success stamps only"
	@echo "  make clean-failed  - remove failed log links only"
	@echo "  make clean-test-results - remove generated eig/lin test results"
	@echo "  make clean-all   - remove all logs, success stamps, and images"
	@echo "  make realclean   - clean-all plus remote tmp/log/Docker and test results cleanup"
	@echo ""
	@echo "Architecture notes:"
	@echo "  - i386/s390x/mips64le: Debian only (QEMU, very slow)"
	@echo "  - arm64/ppc64le: QEMU on non-native hosts"
	@echo "  - amd64: Native on x86_64 hosts"
