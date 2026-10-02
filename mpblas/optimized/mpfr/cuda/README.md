# libmplapack_mpfr_opt_cuda

MPFR flavor of MPLAPACK whose `Rgemm` runs on an NVIDIA GPU for 512-bit and
1024-bit operands. Every other routine, and every `Rgemm` call that is not
eligible, uses the same code as `libmplapack_mpfr_opt`.

The GPU arithmetic is `cu_fp::cu_freal<PB>` from
[mpc_cuda](https://github.com/tkouya/mpc_cuda): a fixed-precision,
register-resident type whose `+ - *` are bit-exact with MPFR round-to-nearest.
The kernels perform the same operations in the same order as
`openmp/Rgemm_*_omp.cpp`, so a GPU result equals the CPU result, except that a
zero result is always `+0` (`cu_freal` has no signed zero).

## Building

mpc_cuda is used header-only (`mpc_cuda/cu_freal.cuh`); pass its include
directory.

CMake:

    cmake -B build -DMPLAPACK_ENABLE_MPFR=ON -DMPLAPACK_ENABLE_OPT=ON \
          -DMPLAPACK_ENABLE_MPFR_CUDA=ON \
          -DMPLAPACK_MPC_CUDA_INCLUDE_DIR=/path/to/mpc_cuda/include \
          -DMPLAPACK_CUDA_ARCHITECTURES=121        # e.g. GB10; default 70;80;90

autotools:

    ./configure --enable-mpfr --enable-mpfr-cuda \
                --with-mpc-cuda-includedir=/path/to/mpc_cuda/include \
                --with-cudatoolkithome=/usr/local/cuda

Link with `-lmplapack_mpfr_opt_cuda -lcudart` (pkg-config:
`mplapack_mpfr_opt_cuda`).

## When the GPU is used

`Rgemm` sends a call to the GPU when all of the following hold; otherwise it
runs on the CPU and the result is unchanged:

- the MPFR default precision, `alpha`, `beta` and every referenced element of
  `A`, `B` and `C` have the same precision, 512 or 1024 bits;
- every referenced value is finite (`C` is not read when `beta == 0`);
- exponents are small enough that no intermediate can leave the MPFR exponent
  range, and every result fits in it;
- `m*n*k >= MPLAPACK_MPFR_CUDA_MIN_MNK` (default 32768);
- a CUDA device is available and `MPLAPACK_MPFR_CUDA` is not `0`.

Environment variables:

| Variable | Meaning |
|---|---|
| `MPLAPACK_MPFR_CUDA=0` | never use the GPU |
| `MPLAPACK_MPFR_CUDA_MIN_MNK=N` | minimum `m*n*k` for the GPU path |
| `MPLAPACK_MPFR_CUDA_VERBOSE=1` | print CUDA errors that cause a CPU fallback |

## Tests

- `mpfr_cuda_Rgemm_host` runs the CUDA kernels' element code on the host and
  compares it bit for bit with the CPU `Rgemm` (512/1024 bits, all transpose
  combinations, several `alpha`/`beta`, zeros and mixed magnitudes), and checks
  that ineligible calls fall back. It needs no GPU and no CUDA toolkit; CMake
  builds it whenever `MPLAPACK_MPC_CUDA_INCLUDE_DIR` is set.
- `mpfr_cuda_Rgemm_device` is the same check with the kernels on the GPU; it
  is skipped when no device is present.

Both run with `OMP_NUM_THREADS=1`: the OpenMP CPU `Rgemm` used as the reference
creates its temporaries in each worker thread at that thread's MPFR default
precision, which is not the caller's precision.

## Limitations

- Only `Rgemm` is accelerated; precisions other than 512 and 1024 bits use the
  CPU.
- One thread per element of `C`, no shared-memory tiling yet.
- mpc_cuda targets CUDA 13. With CUDA 12.0, nvcc rejects the early-clobber
  `"=&l"` asm constraints in `cu_freal.cuh`.
