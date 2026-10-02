# libmplapack_mpfr_opt_cuda

MPFR flavor of MPLAPACK whose `Rgemm` runs on an NVIDIA GPU. Every other
routine, and every `Rgemm` call that is not eligible, uses the same code as
`libmplapack_mpfr_opt`.

The GPU arithmetic comes from [mpc_cuda](https://github.com/tkouya/mpc_cuda)
(bundled as `external/mpc_cuda`):

| Precision | Kernels | Arithmetic | Result vs CPU `Rgemm` |
|---|---|---|---|
| 512, 1024 bits | `Rgemm_device_cuda.cu` | `cu_fp::cu_freal<PB>`: fixed precision, significand in registers | identical, except that a zero is always `+0` |
| any other | `Rgemm_rt_cuda.cu` (port of mpc_cuda `demos/matmul_mpfr.cu`) | `cu_mpfr`: MPFR 4.2.2 compiled for the device | identical, signed zeros included |

Both perform the same operations in the same order as
`openmp/Rgemm_*_omp.cpp`.

## Building

mpc_cuda needs CUDA 13 (nvcc 12.0 rejects its early-clobber `"=&l"` asm
constraints) and is built for a single GPU architecture.

CMake (builds `external/mpc_cuda` unless `MPLAPACK_MPC_CUDA_INCLUDE_DIR` names
an installed mpc_cuda):

    cmake -B build -DMPLAPACK_ENABLE_MPFR=ON -DMPLAPACK_ENABLE_OPT=ON \
          -DMPLAPACK_ENABLE_MPFR_CUDA=ON \
          -DMPLAPACK_CUDA_ARCHITECTURES=121        # one value, e.g. 121 for GB10
    # installed mpc_cuda instead of the bundled one:
    #     -DMPLAPACK_MPC_CUDA_INCLUDE_DIR=/opt/mpc_cuda/include
    #     [-DMPLAPACK_MPC_CUDA_LIBRARY=/opt/mpc_cuda/lib/libmpc_cuda.a]

autotools (builds `external/mpc_cuda` unless `--with-mpc-cuda-includedir`
names an installed mpc_cuda):

    ./configure --enable-mpfr --enable-mpfr-cuda --with-mpfr-cuda-arch=121 \
                --with-cudatoolkithome=/usr/local/cuda
    # installed mpc_cuda instead of the bundled one:
    #     --with-mpc-cuda-includedir=/opt/mpc_cuda/include
    #     [--with-mpc-cuda-library=/opt/mpc_cuda/lib/libmpc_cuda.a]

The CUDA sources are compiled as relocatable device code and device-linked
with `libmpc_cuda.a` inside the library, so programs link it with the C++
compiler: `-lmplapack_mpfr_opt_cuda -lcudart` (pkg-config:
`mplapack_mpfr_opt_cuda`).

## When the GPU is used

`Rgemm` sends a call to the GPU when all of the following hold; otherwise it
runs on the CPU and the result is unchanged:

- the MPFR default precision, `alpha`, `beta` and every referenced element of
  `A`, `B` and `C` have the same precision;
- the MPFR default rounding mode is `MPFR_RNDN`;
- every referenced value is finite (`C` is not read when `beta == 0`);
- exponents are small enough that no intermediate can leave the MPFR exponent
  range (host and device), and every result fits in it;
- `m*n*k >= MPLAPACK_MPFR_CUDA_MIN_MNK` (default 32768);
- a CUDA device is available and `MPLAPACK_MPFR_CUDA` is not `0`.

Environment variables:

| Variable | Meaning |
|---|---|
| `MPLAPACK_MPFR_CUDA=0` | never use the GPU |
| `MPLAPACK_MPFR_CUDA_MIN_MNK=N` | minimum `m*n*k` for the GPU path |
| `MPLAPACK_MPFR_CUDA_RUNTIME=0` | no GPU for precisions other than 512/1024 |
| `MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1` | use the `cu_mpfr` kernels for 512/1024 bits too |
| `MPLAPACK_MPFR_CUDA_RT_BLOCKS`, `MPLAPACK_MPFR_CUDA_RT_THREADS` | thread pool of the `cu_mpfr` kernels (default 256 x 32) |
| `MPLAPACK_MPFR_CUDA_VERBOSE=1` | print CUDA errors that cause a CPU fallback |

The `cu_mpfr` kernels take the limbs MPFR allocates internally from a
per-thread arena of 16 KB per 1024 bits of precision, reset at every step of
the dot product (as in mpc_cuda's matmul demo); they raise the device stack
limit to 64 KB while they run.

## Tests (CMake)

`mpfr_cuda_Rgemm_{host,device}[_runtime]` compare the GPU `Rgemm` bit for bit
with the CPU `Rgemm` at 64, 200, 256, 333, 512, 1024 and 2048 bits, for all
transpose combinations and several `alpha`/`beta`, with zeros and mixed
magnitudes, and check that ineligible calls fall back. `_runtime` runs every
precision through the `cu_mpfr` kernels.

- `host`: the kernels' element code runs on the host, with `cu_mpfr` mapped
  to the system MPFR (libmpc_cuda's own host code needs a GPU). No GPU needed.
- `device`: the kernels run on the GPU; skipped when no device is present.

They run with `OMP_NUM_THREADS=1`: the OpenMP CPU `Rgemm` used as the
reference creates its temporaries in each worker thread at that thread's MPFR
default precision, which is not the caller's precision.

## Limitations

- Only `Rgemm` is accelerated.
- One thread per element of `C`, no shared-memory tiling yet.
- Not yet measured on a GPU.
