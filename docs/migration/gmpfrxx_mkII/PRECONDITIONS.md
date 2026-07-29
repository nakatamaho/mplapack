# P0 Preconditions

## Locked repositories

- MPLAPACK repository: `git@github.com:nakatamaho/mplapack.git`.
- Migration branch: `topic/gmpfrxx_mkII_migration`.
- Verified `origin/master` and immutable P0 base:
  `b875e74d4b927282c907c3a29e6cadda48a7d57b`.
- Branch merge base with `origin/master`: the same commit; ahead/behind at
  branch creation was `0/0`.
- gmpfrxx_mkII repository:
  `git@github.com:nakatamaho/gmpfrxx_mkII.git`.
- Locked upstream adapter base:
  `2f06785c3f1b62f92e1e2026c2c975df73d1e426` (`origin/main`).
- Bootstrap release tag `v1.0.1` resolves to
  `7e912f300443dd1f53dd69285e0c8302a7cdb3b3`.

## Prompt and prototype provenance

The canonical prompt pack is
`docs/migration/gmpfrxx_mkII/CODEX_PROMPTS.md`, identifies itself as v5.4,
and has SHA-256
`c4175afa3ec59c385f00a54c7b762b09c0f05b7c8418b33e940ce340d46d2c0e`.
It and the maintainer decisions are authoritative.

The external historical source filename is
`FABLE_PROTOTYPE_PROMPT.v4.md`, at
`/home/docker/mplapack-gmpfrxx-v5.4.1-input/FABLE_PROTOTYPE_PROMPT.v4.md`,
with SHA-256
`d14b03964672da96bb73c86960913d9a566a209bcba86379083bd71c56f262a4`.
It was imported byte-for-byte as
`docs/migration/gmpfrxx_mkII/INTEROP_PROTOTYPE.md`; the imported SHA-256 is
the same, and `cmp -s` returned 0.

The v4 file is provenance, not an executable specification. Its phase-specific
MPLAPACK branches, reverse conversions, mixed adapter arithmetic, both operand
orders, compound adapter assignment, round-trip/exact-rounding requirements,
precision sweeps, and edd/td support are obsolete and must not be executed.

## MPLAPACK build layouts

At the locked SHA, Autotools is rooted in `configure.ac` and the committed
`Makefile.am` tree. A Git checkout is bootstrapped with
`autoreconf --force --install`; `gen_configure.sh` invokes the individual
Autotools programs but assumes libtool support files are already present.
The supported full configuration enables backends with
`--enable-{gmp,mpfr,qd,dd,double,binary80,binary128}=yes` and tests with
`--enable-test=yes`.

Bundled dependencies are committed under `external/*/download/`, verified by
`external/distfiles.sha256`, and use separate
`external/*/work/internal` and `external/*/work/install` source trees.
The internal GMP, MPFR, MPC, and QD installations are staged under
`external/i/<PACKAGE>`; the final tree installs to MPLAPACK's `$(prefix)`.
The locked versions are GMP 6.3.0, MPFR 4.2.2, MPC 1.4.1, and QD 2.3.24.
The gmpfrxx_mkII P0 compile installation was configured with upstream
auto-fetch disabled and with all three dependency paths pointing at this
MPLAPACK-owned internal prefix.

CMake is rooted in `CMakeLists.txt` and `cmake/`. It is a dependency consumer:
backend selection is controlled by `MPLAPACK_ENABLE_*` options, dependencies
are found with `find_package`, and it does not build or download GMP, MPFR,
MPC, QD, or gmpfrxx_mkII. Tests are controlled by
`MPLAPACK_BUILD_TESTS`. Installed package metadata is produced by
`cmake/MplapackInstall.cmake` and `cmake/MplapackPackage.cmake`.

## Generated-file policy

Per root `AGENTS.md`, routine bodies under `mpblas/reference/` and
`mplapack/reference/` are FABLE-generated and are read-only. Numerical changes
must go through the FABLE pipeline. Generated backend declarations under
`include/mplapack_*_<backend>.h` and generated test sources are inventoried as
generated outputs, not hand-edit targets. Source-tree
`include/mplapack_config.h` may be stale; CMake's build-tree generated header
must win by include order. The dd backend requires `-ffp-contract=off` on every
translation unit.

## Locked gmpfrxx_mkII public surface

Header roles at `2f06785...` are:

- `gmpxx_mkII.h`: GMP integer, rational, MPF, and MPFC wrapper family.
- `mpfrxx_mkII.h`: MPFR real wrapper, MPFR environment, and GMP base types.
- `mpcxx_mkII.h`: MPC complex wrapper layered on MPFR and GMP types.
- `gmpfrxx_mkII.h`: aggregate header for all enabled components.
- `gmpfrxx_mkII/detail/*.hpp`: implementation, environment, evaluation,
  math, and generated build/version configuration.

The locked tree contains no `include/gmpfrxx_mkII/adapters/` directory and no
adapter headers. That contradicts the historical v4 document's statement that
prototype adapters already exist.

Installed ordinary CMake package targets are
`gmpxx_mkII::gmpxx_mkII`, `mpfrxx_mkII::mpfrxx_mkII`,
`mpcxx_mkII::mpcxx_mkII`, and
`gmpfrxx_mkII::gmpfrxx_mkII`. The locked release also installs
`gmpxx_mkII::default_context_provider` and its shared library. P0 identified
an existing non-512 runtime-default dependency in MPLAPACK. The maintainer
accepted retirement of that variable-precision behavior for the migration, so
the provider must not be selected or added solely for MPLAPACK.

The correctness configuration requires the following compile-time semantic
macros to remain undefined/OFF in every translation unit:
`GMPFRXX_MKII_ENABLE_FMA`, `GMPFRXX_MKII_FAST_FIXED_PREC`, and
`GMPFRXX_MKII_FAST_STABLE_RND`.

## Default precision

gmpfrxx_mkII establishes the builtin MPFR default in
`include/gmpfrxx_mkII/detail/environment.hpp` and the builtin GMP MPF default
in `include/gmpfrxx_mkII/detail/gmp_default_context.hpp`; both are 512 binary
bits. The ordinary APIs are `mpfrxx::default_precision_bits()`,
`mpfrxx::set_default_precision_bits()`,
`gmpxx::default_mpf_precision_bits()`, and
`gmpxx::set_default_mpf_precision_bits()`.

The legacy MPLAPACK wrappers establish MPFR's default through
`___MPREAL_DEFAULT_PRECISION___` in `mpfrc++/mpreal.h` and GMP MPF defaults at
call sites using `mpf_set_default_prec`; the inventory classifies every setter.

The locked MPLAPACK tree contains the following current variable-precision
behavior:

- `mpblas/reference/mplapackinit.cpp` applies `MPLAPACK_GMP_PRECISION` through
  `mpf_set_default_prec()` and applies `MPLAPACK_MPFR_PRECISION` plus
  `MPLAPACK_MPC_REAL_PRECISION`/`MPLAPACK_MPC_IMAG_PRECISION` to the legacy
  MPFR/MPC wrapper defaults.
- `mplapack/test/lin/mpfr/Makefile.am` runs MPFR profiles at 53 and 113 bits,
  and `mplapack/test/compare/common/Rlamch.test.cpp` mutates GMP default
  precision through 4096, 128, and 64 bits.
- `MIGRATION.md` and `CHANGES.2.1.0.md` describe the MPFR environment variable
  as supported user-facing behavior.

The maintainer explicitly decided during P0 that migration to gmpfrxx_mkII may
retire precision changes performed by `mplapackinit.cpp`. Under the binding
512-bit policy, later MPLAPACK phases must also remove or replace the existing
non-512 default-mutation test profiles rather than carry them into the new
wrapper implementation. MPFR exponent-range configuration is a separate
feature and remains in scope, but later implementation must establish the
normal 512-bit wrapper default before exponent-range changes can prevent
gmpfrxx_mkII first-use initialization. No provider, cross-DSO synchronization,
precision sweep, or runtime default-mutation compatibility layer is approved.

This accepted maintainer decision resolves the Task 0 blocker. The required
single fresh-process smoke will default-construct one `mpfrxx::mpfr_class` and
one `gmpxx::mpf_class`, call no setter, and require both object precisions to
be exactly 512.

## Material discrepancies

The historical v4 document conflicts with the canonical v5.4 decisions in the
following material ways:

- multiple phase-specific MPLAPACK branches versus one migration branch;
- reverse `cast2*` conversion and round-trip goals versus one-way comparison
  embedding only;
- general mixed adapter arithmetic and expression composition versus explicit
  materialization before MPFR/MPC-only comparison;
- exact-rounding and precision-sweep requirements versus practical comparison
  quality at the ordinary 512-bit default;
- edd/td adapters versus explicit exclusion;
- a generic allowed-improvement baseline notion versus immutable,
  metric-specific directional rules;
- claimed existing adapter files versus none at the locked upstream SHA.

The obsolete v4 requirements are recorded and are not executed. The locked
repository also exposed existing variable-precision behavior that conflicted
with the fixed 512-bit migration policy. The maintainer resolved that mismatch
by explicitly approving retirement of `mplapackinit.cpp` precision changes and
non-512 default-mutation profiles during migration. No material contradiction
remains after recording that decision.
