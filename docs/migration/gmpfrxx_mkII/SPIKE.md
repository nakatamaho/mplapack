# P0 Compile Spike

## Locked inputs

- MPLAPACK source: `/tmp/mplapack-p0-base`, exact base
  `b875e74d4b927282c907c3a29e6cadda48a7d57b`.
- MPLAPACK Autotools install: `/tmp/mplapack-p0-final`.
- gmpfrxx_mkII install: `/tmp/gmpfrxx-p0-install`, exact source
  `2f06785c3f1b62f92e1e2026c2c975df73d1e426`.
- Compiler: GCC C++ 15.2.0, `-std=gnu++17 -fopenmp
  -ffp-contract=off`.
- Correctness macros `GMPFRXX_MKII_ENABLE_FMA`,
  `GMPFRXX_MKII_FAST_FIXED_PREC`, and
  `GMPFRXX_MKII_FAST_STABLE_RND` were not defined.

Command:

```text
python3 docs/migration/gmpfrxx_mkII/tools/spike/run_spikes.py \
  --mplapack-prefix /tmp/mplapack-p0-final \
  --mplapack-source /tmp/mplapack-p0-base \
  --gmpfrxx-prefix /tmp/gmpfrxx-p0-install \
  --work-dir /tmp/mplapack-p0-spikes \
  --results docs/migration/gmpfrxx_mkII/spike/results.json
```

The command passed because the mandatory fresh-process default smoke passed
and every other failure has an explicit later owner. Complete compiler and
runtime output is in `spike/spike-logs/`.

## Wrapper surface

| Case | Result | Finding | Owner |
|---|---|---|---|
| `default-precision-512` | PASS | One fresh process default-constructed `mpfrxx::mpfr_class` and `gmpxx::mpf_class`; both reported exactly 512 bits before any setter. | P0 |
| `same-family-arithmetic` | PASS | Ordinary MPFR/MPC and GMP MPF/MPFC arithmetic compiled and ran. | P2C verification |
| `utility-print-hex` | PASS | Replacement real/complex math, decimal stream output, and hexadecimal stream output compiled and ran. | P2C verification |
| `mplapack-mpfr-surface` | COMPILE_FAIL | Locked overloads `pow2`, `sign`, `nint`, `pi`, `sprintnum`, and `sprinthex_mpfr_fixed` accept legacy `mpfr::mpreal`, not `mpfrxx::mpfr_class`. | P2C/P3 |
| `mplapack-gmp-surface` | COMPILE_FAIL | Locked overloads `pow2`, `sign`, `nint`, `pi`, and `sprintnum` accept global `::mpf_class`, not `gmpxx::mpf_class`. | P2C/P4 |
| `legacy-name-collision` | PASS | Global legacy `::mpc_class`, `mpfrxx::mpc_class`, and `gmpxx::mpfc_class` coexist when fully qualified. | P2C verification |
| `expression-lifetime` | PASS | Nested same-family MPFR/MPC expressions formed from temporary operands materialized correctly. | P2C verification |

The utility gaps are ordinary REAL/COMPLEX wrapper-family work. They do not
justify mixed-backend adapter arithmetic.

## Required real embeddings

| Matrix ID | Result | Exact gap or quality result | Owner |
|---|---|---|---|
| `REAL-DOUBLE` | PASS | `mpfrxx::mpfr_class(double)` compiled and preserved `1.25`. | P2A verification |
| `REAL-DD` | COMPILE_FAIL | No matching `mpfrxx::mpfr_class(dd_real&)` constructor. | P2A |
| `REAL-QD` | COMPILE_FAIL | No matching `mpfrxx::mpfr_class(qd_real&)` constructor. | P2A |
| `REAL-BINARY80` | COMPILE_FAIL | No binary80 constructor; overload resolution reaches the deleted `mpfr_class(bool)` path. | P2A |
| `REAL-BINARY128` | COMPILE_FAIL | No binary128 constructor; overload resolution reaches the deleted `mpfr_class(bool)` path. | P2A |
| `REAL-GMP` | COMPILE_FAIL | No matching `mpfrxx::mpfr_class(gmpxx::mpf_class&)` constructor. | P2A |

The dd and qd sources use `1 + 2^-100`. Binary80 uses `1 + 2^-60`, and
binary128 uses `1 + 2^-100`. A passing implementation must leave a positive
MPFR difference after subtracting one, so an extended source cannot silently
pass through binary64.

## Required complex embeddings

| Matrix ID | Result | Exact gap or quality result | Owner |
|---|---|---|---|
| `COMPLEX-DOUBLE` | PASS | `mpfrxx::mpc_class(std::complex<double>)` compiled and preserved both components. | P2B verification |
| `COMPLEX-DD` | RUN_FAIL | Construction compiled, but the `2^-100` real component was lost; the anti-binary64-fallback check returned 1. | P2B |
| `COMPLEX-QD` | RUN_FAIL | Construction compiled, but the `2^-100` real component was lost; the anti-binary64-fallback check returned 1. | P2B |
| `COMPLEX-BINARY80` | COMPILE_FAIL | No matching `mpfrxx::mpc_class(std::complex<mplapack_binary80_t>&)` constructor. | P2B |
| `COMPLEX-BINARY128` | COMPILE_FAIL | No matching `mpfrxx::mpc_class(std::complex<mplapack_binary128_t>&)` constructor. | P2B |
| `COMPLEX-GMP` | COMPILE_FAIL | No matching `mpfrxx::mpc_class(gmpxx::mpfc_class&)` constructor. | P2B |

For dd/qd, retaining information around `2^-100` is only the practical
comparison-quality smoke required by the migration. It is not an exact
100-bit conversion or exact-rounding contract.

## Scope boundaries

The spike does not test source/destination precision metadata, alternate
default precision, runtime default mutation, reverse conversion, mixed
backend arithmetic, both operand orders, adapter compound assignment,
round-trip identity, exact rounding, ULP guarantees, edd, or td. Those
operations are forbidden, not used, or out of scope in
`interop_requirements.tsv`.
