# MPLAPACK Release Engineering

Multi-environment build and test tools for release validation.

## Quick Start

```bash
cd release/

# Full test cycle (branch → make dist → tarball)
make

# Branch tests only
make branch

# Tarball tests only (with existing tarball)
make tarball-with TARBALL=~/mplapack-2.1.0.tar.gz
```

## Build Matrix

Defined in `build-matrix.conf`. Format:
```
name|base_image|archs|dockerfile|test_cmd|source_type
```

### Source Types

- `branch`: Build from git source (requires autoreconf)
- `tarball`: Build from release tarball (no autoreconf)

### Supported Configurations

| Category | Environments |
|----------|-------------|
| Debian-based | Ubuntu 22.04, 24.04, Debian 12 |
| RHEL-based | Fedora 39, 40, Rocky Linux 8, 9 |
| SUSE-based | openSUSE Leap 15.5, Tumbleweed |
| musl libc | Alpine 3.19 |
| Compilers | GCC (default), Intel oneAPI, MinGW (cross) |
| CUDA | 11.8, 12.1, 12.4 (Ubuntu 22.04/24.04) |

### Supported Architectures

| Architecture | Docker Platform | Supported Distros | Notes |
|--------------|-----------------|-------------------|-------|
| x86_64 | linux/amd64 | All | Native on x86_64 hosts |
| aarch64 | linux/arm64 | All | QEMU on non-native hosts |
| ppc64le | linux/ppc64le | Debian | PowerPC 64-bit LE, QEMU |
| i386 | linux/386 | Debian | x86 32-bit, QEMU |
| s390x | linux/s390x | Debian | IBM Z mainframe, QEMU |
| mips64le | linux/mips64le | Debian | MIPS 64-bit LE, QEMU |

## Usage

### Full Release Validation

```bash
make
```

This runs:
1. Branch tests on all configurations
2. `make dist` to generate tarball
3. Tarball tests on all configurations

### Filtered Runs

```bash
# By distro/toolchain
make ubuntu
make debian
make fedora
make rocky
make suse
make alpine
make intel
make mingw
make cuda

# By architecture
make amd64
make arm64

# Tier 1 remote tests
make tier1
make tier1-macos
make tier1-linux
make tier1-macos-amd64   # SSH to macOS amd64 host and run Tier 1 buildtest
make tier1-macos-arm64   # SSH to macOS arm64 host and run Tier 1 buildtest
make tier1-ubuntu-arm64  # SSH to arm64 host and run Ubuntu arm64 Docker distcheck
make ppc64le
make i386
make s390x
make mips64le

# Combined
make branch-ubuntu
make branch-debian
make branch-rocky
make branch-suse
make branch-cuda
make tarball-amd64
make tarball-debian
```

### Testing Existing Tarball

```bash
make tarball-with TARBALL=/path/to/mplapack-2.1.0.tar.gz
```

## Output

Logs are saved to `logs/YYYYMMDD_HHMMSS/`:

```
logs/20250201_143000/
├── dist/
│   ├── autoreconf.log
│   ├── configure.log
│   ├── make_dist.log
│   └── mplapack-2.1.0.tar.gz
├── results.csv
├── results_branch.csv
├── results_tarball.csv
├── ubuntu22_amd64_build.log
├── ubuntu22_amd64_test.log
├── debian12_s390x_build.log
├── debian12_s390x_test.log
├── cuda124-ubuntu22_amd64_build.log
├── cuda124-ubuntu22_amd64_test.log
└── ...
```

### Results CSV Format

```csv
name,arch,base,stage,result,elapsed,source_type
ubuntu22,amd64,ubuntu:22.04,test,OK,342,branch
debian12,s390x,debian:12,test,OK,4521,branch
cuda124-ubuntu22,amd64,nvidia/cuda:12.4.0-devel-ubuntu22.04,test,OK,567,branch
rocky9,arm64,rockylinux:9,test,FAILED,892,branch
```

## Prerequisites

### Docker with buildx

```bash
# Verify buildx is available
docker buildx version
```

### QEMU for Cross-Architecture Builds

Required for arm64, ppc64le, i386, s390x, mips64le on non-native hosts:

```bash
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
```

### CUDA Requirements

- **Build**: No GPU required (uses nvcc for compilation only)
- **Test**: NVIDIA GPU + drivers required on host
- NVIDIA Container Toolkit for `--gpus` flag

```bash
# Install NVIDIA Container Toolkit (Ubuntu)
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
    sudo tee /etc/apt/sources.list.d/nvidia-docker.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

GPU auto-detection:
- Script automatically detects GPU availability
- If GPU found: tests run with `--gpus all`
- If no GPU: tests run without GPU (may fail for GPU-specific tests)

Override with environment variable:
```bash
USE_GPU=yes make cuda   # Force GPU usage
USE_GPU=no make cuda    # Force no GPU (build-only validation)
```

## Adding New Configurations

Edit `build-matrix.conf`:

```conf
# name|base|archs|dockerfile|test_cmd|source_type
rocky9|rockylinux:9|linux/amd64,linux/arm64|Dockerfile.redhat|make check|branch
cuda130-ubuntu24|nvidia/cuda:13.0.0-devel-ubuntu24.04|linux/amd64|Dockerfile.cuda|make check|branch
```

Create corresponding Dockerfile in `docker/` if needed.

### Dockerfile Naming Convention

| Dockerfile | Description |
|------------|-------------|
| Dockerfile.debian | Debian/Ubuntu, branch build |
| Dockerfile.debian-tarball | Debian/Ubuntu, tarball build |
| Dockerfile.redhat | Fedora/Rocky 9+, branch build |
| Dockerfile.redhat-tarball | Fedora/Rocky 9+, tarball build |
| Dockerfile.redhat-el8 | Rocky 8 / RHEL 8, branch build |
| Dockerfile.redhat-el8-tarball | Rocky 8 / RHEL 8, tarball build |
| Dockerfile.suse | openSUSE, branch build |
| Dockerfile.suse-tarball | openSUSE, tarball build |
| Dockerfile.alpine | Alpine, branch build |
| Dockerfile.alpine-tarball | Alpine, tarball build |
| Dockerfile.intel | Intel oneAPI, branch build |
| Dockerfile.mingw | MinGW cross-compile |
| Dockerfile.cuda | CUDA, branch build |
| Dockerfile.cuda-tarball | CUDA, tarball build |

## Troubleshooting

### Build fails on arm64/ppc64le/s390x/mips64le

QEMU emulation is slow and may timeout. These architectures can take 10-50x longer than native builds.

Check logs:
```bash
cat logs/*/debian12_s390x_build.log
```

Increase Docker timeout if needed:
```bash
export DOCKER_CLIENT_TIMEOUT=600
export COMPOSE_HTTP_TIMEOUT=600
```

### i386 build fails

32-bit builds may have issues with large file support or memory addressing. Check for specific compiler flags needed.

### Rocky 8 package installation fails

Rocky 8 requires PowerTools repository for development packages. The `Dockerfile.redhat-el8` handles this automatically.

### Intel compiler license issues

Intel oneAPI base toolkit is free but some features require registration. The base compilers (icx, icpx, ifx) should work without registration.

### MinGW tests fail

Wine configuration may need initialization:
```bash
docker run --rm mplapack-ubuntu22-mingw64-amd64 winecfg
```

### CUDA build succeeds but tests fail

1. Check GPU availability:
```bash
docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
```

2. Verify NVIDIA Container Toolkit:
```bash
nvidia-ctk --version
```

3. Check driver compatibility with CUDA version:
   - CUDA 11.8: Driver >= 450.80.02
   - CUDA 12.1: Driver >= 525.60.13
   - CUDA 12.4: Driver >= 550.54.14

### Out of disk space

Multi-arch builds with many configurations consume significant disk space.

```bash
# Clean up build artifacts
make clean

# Remove all Docker build cache
docker builder prune -af

# Remove unused images
docker image prune -af
```

### QEMU not working

Re-register QEMU handlers:
```bash
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes
```

Verify registration:
```bash
cat /proc/sys/fs/binfmt_misc/qemu-aarch64
```

## Performance Considerations

### Estimated Build Times (on x86_64 host)

| Architecture | Relative Speed | Typical Time |
|--------------|----------------|--------------|
| amd64 (native) | 1x | 5-10 min |
| arm64 (QEMU) | 5-10x | 30-60 min |
| ppc64le (QEMU) | 10-20x | 60-120 min |
| i386 (QEMU) | 3-5x | 15-30 min |
| s390x (QEMU) | 10-30x | 60-180 min |
| mips64le (QEMU) | 20-50x | 120-300 min |

### Recommendations

1. **Quick validation**: Test amd64 only first
   ```bash
   make branch-amd64
   ```

2. **Pre-release**: Add arm64
   ```bash
   make amd64 arm64
   ```

3. **Full release**: Run all architectures overnight
   ```bash
   make
   ```

4. **CI/CD**: Consider running exotic architectures (s390x, mips64le) on schedule rather than every commit

## Directory Structure

```
release/
├── README.md              # This file
├── Makefile               # Build targets
├── build-all.sh           # Main build script
├── build-matrix.conf      # Build configuration matrix
├── .gitignore             # Ignore logs and artifacts
├── docker/
│   ├── Dockerfile.debian
│   ├── Dockerfile.debian-tarball
│   ├── Dockerfile.redhat
│   ├── Dockerfile.redhat-tarball
│   ├── Dockerfile.redhat-el8
│   ├── Dockerfile.redhat-el8-tarball
│   ├── Dockerfile.suse
│   ├── Dockerfile.suse-tarball
│   ├── Dockerfile.alpine
│   ├── Dockerfile.alpine-tarball
│   ├── Dockerfile.intel
│   ├── Dockerfile.mingw
│   ├── Dockerfile.cuda
│   └── Dockerfile.cuda-tarball
└── logs/                  # Build logs (gitignored)
    └── YYYYMMDD_HHMMSS/
        ├── dist/
        ├── results.csv
        └── *.log
```

## License

Same as MPLAPACK main project.
