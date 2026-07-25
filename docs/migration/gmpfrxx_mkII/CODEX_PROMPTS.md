# MPLAPACK -> gmpfrxx_mkII migration — Codex prompt pack v5.4 (single MPLAPACK migration branch)

Canonical in-repo location: `docs/migration/gmpfrxx_mkII/CODEX_PROMPTS.md`.
Supersedes v5.3, v4, and the initial v5 reviewed draft.

This revision incorporates the maintainer's binding interoperability and
precision decisions:

- interoperability exists only to materialize backend results in MPFR/MPC
  reference objects for comparison;
- `mpfrxx::mpfr_class` and `gmpxx::mpf_class` use a 512-bit default in both
  gmpfrxx_mkII and MPLAPACK; the corresponding complex wrappers use the normal
  512-bit component defaults supplied by those wrapper families;
- conversion and assignment do not inspect, compare, negotiate, or test source
  and destination precision metadata;
- reverse MPFR/MPC -> GMP/other-backend conversion, mixed adapter arithmetic,
  exact-rounding proofs, and precision sweeps are not migration requirements;
- every MPLAPACK phase runs on the single maintainer-selected branch
  `topic/gmpfrxx_mkII_migration`, created once from a verified `origin/master`
  commit; phase boundaries are immutable accepted commit SHAs on that branch,
  not separate phase branches.

This prompt pack is executable only against the repository revisions recorded
by P0 in `docs/migration/gmpfrxx_mkII/LOCK.json`. The branch name alone is
never a sufficient phase base: every phase after P0 must verify that HEAD is
exactly the accepted SHA of its required preceding MPLAPACK phase.

## Binding architectural decisions

1. Target release: MPLAPACK 3.0.0.
2. A C++ ABI break is allowed for the MPFR and GMP backend libraries only.
   The public ABI and exported-symbol set of dd, qd, double, binary80, and
   binary128 must not change.
3. No stubs, placeholder implementations, disabled assertions, or temporary
   binary64 fallbacks. If a phase is blocked, stop and report the blocker.
4. gmpfrxx_mkII comparison-embedding code lives upstream under opt-in adapter
   headers. MPLAPACK vendors an official gmpfrxx_mkII release archive; the
   normal `external/gmpfrxx_mkII/patches/` state is empty.
5. edd and td are outside the MPLAPACK migration scope. Existing upstream
   support may remain untouched, but no phase extends, wires, or tests it for
   MPLAPACK.
6. Autotools and CMake remain supported, but their dependency models are not
   the same:
   - Autotools uses the committed archive under `external/gmpfrxx_mkII/` and
     the established internal/final two-prefix external-package flow.
   - MPLAPACK CMake remains `find_package`-based and does not download or build
     GMP, MPFR, MPC, QD, or gmpfrxx_mkII. CMake tests may point
     `CMAKE_PREFIX_PATH` either at the Autotools-staged internal prefix or at an
     independently installed gmpfrxx_mkII prefix.
7. MPLAPACK-specific names such as `sign`, `nint`, `iceil`, `cabs1`, `pow2`,
   `pow4`, `castREAL_*`, and `castINTEGER_*` remain in MPLAPACK utility
   headers. Upstream gmpfrxx_mkII provides generally reusable wrapper APIs and
   the narrow comparison embeddings; it must not include MPLAPACK headers or
   depend on `mplapackint`.
8. Interoperability is one-way and explicit. MPLAPACK comparison code
   materializes finite results from double, dd, qd, binary80, binary128, and
   GMP/MPFC into `mpfrxx::mpfr_class` or `mpfrxx::mpc_class`, then performs
   subtraction, norms, tolerances, and pass/fail decisions entirely in
   MPFR/MPC. Reverse conversion and implicit mixed arithmetic are out of scope.
   In particular, `mpfrxx::mpfr_class -> gmpxx::mpf_class` and
   `mpfrxx::mpc_class -> gmpxx::mpfc_class` are not required.
9. The canonical default precision is exactly 512 binary bits for
   `mpfrxx::mpfr_class` and `gmpxx::mpf_class` in both projects. Default-
   constructed comparison targets therefore use 512 bits. The migration does
   not query source precision, compare precision metadata, negotiate a common
   precision, or add per-assignment precision tests. Ordinary wrapper
   construction/assignment semantics are accepted.
10. Do not introduce a process-wide default-context provider, cross-DSO
    precision synchronization, precision mutation tests, or worker-thread
    precision setup solely for this migration. P0 inventories existing
    runtime precision setters. If current MPLAPACK genuinely depends on a
    non-512 process-wide runtime default mutation across library images, stop
    and report that concrete dependency for a maintainer decision.
11. Correctness builds use the same gmpfrxx_mkII compile-time semantic options
    in every translation unit and final linked image. Unless a later approved
    performance phase says otherwise:
    - `GMPFRXX_MKII_ENABLE_FMA` is OFF;
    - `GMPFRXX_MKII_FAST_FIXED_PREC` is OFF;
    - `GMPFRXX_MKII_FAST_STABLE_RND` is OFF.

## Branch and integration model

### MPLAPACK repository

- The only migration branch is `topic/gmpfrxx_mkII_migration`.
- P0 creates that branch once from the exact verified `origin/master` SHA and
  records the SHA in `LOCK.json`. Do not use a local moving `master` as the
  base without first fetching and resolving `origin/master` to a full SHA.
- P1, P1R, P2E, and P3-P8 continue on the same branch. Do not create, merge,
  or switch to phase-specific `topic/gmpfrxx-*` branches.
- Every phase ends at a committed candidate SHA. The maintainer accepts that
  SHA before the next phase begins. The next phase verifies that branch HEAD is
  exactly the accepted SHA; branch continuity never replaces SHA verification.
- Do not rebase, force-reset, merge `master`, or rewrite accepted phase commits
  unless the maintainer explicitly instructs it. New upstream changes are not
  incorporated opportunistically during the migration.
- The intended integration unit is one MPLAPACK migration branch and one final
  MPLAPACK pull request, with per-phase reports and accepted SHAs inside it.
- If a local or remote branch with this name already exists, never overwrite
  it. Fetch it, inspect its ancestry and worktree state, and stop on any
  mismatch with the expected locked base or accepted phase SHA.

### gmpfrxx_mkII repository

P2A-P2D are upstream work in a different repository. Their branch names remain
local to that repository and are listed in the P2 prompts. The maintainer's
`topic/gmpfrxx_mkII_migration` single-branch rule applies to MPLAPACK phases.
P2D produces the official upstream release consumed by MPLAPACK P1R.

## Repository ownership

- `nakatamaho/mplapack`: P0, P1, P1R, P2E, P3-P8, all on
  `topic/gmpfrxx_mkII_migration`.
- `nakatamaho/gmpfrxx_mkII`: P2A-P2D on their upstream branches.

## Phase graph

```text
P0 ──► P1 ───────────────────────────────┐
 │                                       │
 └────► P2A ► P2B ► P2C ► P2D ─────────┤  (gmpfrxx_mkII)
                                         ▼
                                        P1R ► P2E ► P3 ► P4 ► P5 ► P6 ► P7 ► P8
```

P1 may run in parallel with P2A-P2D. P1R must vendor the P2D release. P2E must
finish before either backend switch. P3 and P4 must use generator-driven
changes for generated public headers; hand-editing those headers is forbidden.

## Required phase-report format

Phase reports use repository-local paths:

- MPLAPACK phases: `docs/migration/gmpfrxx_mkII/REPORT-<phase>.md`;
- gmpfrxx_mkII phases: `docs/mplapack_migration/REPORT-<phase>.md`.

Each report contains, in this order:

1. repository, branch, approved base SHA, final SHA;
2. relevant entries and SHA-256 digest of `LOCK.json` or the imported
   cross-repository requirements manifest;
3. files changed, grouped by purpose;
4. exact gate command and complete output;
5. test/build environment and dependency versions;
6. numerical or ABI deviations, including zero-deviation statements;
7. remaining blockers or `none`;
8. worktree status and remote branch SHA.

Each phase commits an executable `tools/gate-<phase>.sh` (or a Python gate plus
thin shell entry point). The single acceptance command is that gate. A gate
must use `set -euo pipefail`, fail on missing tools or skipped tests, and write
its full log into the report directory. Merely showing a grep count or a
`git diff --stat` is not a gate.

## COMMON CONTEXT (prepend verbatim to every Codex session)

```text
You are executing one phase of
`docs/migration/gmpfrxx_mkII/CODEX_PROMPTS.md` v5.4.

Read `AGENTS.md` at the repository root before making changes and obey it.
Read the phase prompt, the local `LOCK.json` or imported lock excerpt, the
immediately preceding accepted phase report, and every requirements manifest
named by the phase. Verify that HEAD is exactly the approved base SHA before
editing. If it is not, stop and report.

Migration-wide invariants:

I1. No stubs, TODO implementations, disabled tests, unconditional skips,
    placeholder return values, or silent binary64 fallback. If blocked, write
    the blocker into the report and stop with a clean worktree.
I2. Never hand-edit Fable-generated routine bodies or generated public routine
    declarations. Change generators/templates and regenerate through the
    repository pipeline. P2E defines the generated-file manifest.
I3. gmpfrxx_mkII arithmetic operators may return expression nodes. Do not store
    an expression in `auto` when it can outlive or observe later mutation of a
    borrowed operand. Materialize explicitly into the owning wrapper type.
I4. Raw `long double` and raw `__float128` remain rejected scalar leaves. Use
    the explicit binary80 and binary128 wrapper types at adapter boundaries.
I5. No multiprecision conversion or math path may pass through binary64, except
    when the source itself is `double`.
I6. edd and td are out of scope for MPLAPACK.
I7. Interop is explicit, one-way comparison embedding into MPFR/MPC. Do not add
    reverse conversion, cross-family mixed arithmetic, both-operand-order
    overloads, compound assignment, or expression-template composition solely
    for interop.
I8. The default precision of `mpfrxx::mpfr_class` and `gmpxx::mpf_class` is
    512 binary bits in gmpfrxx_mkII and MPLAPACK. Conversion/assignment code
    must not query, compare, negotiate, preserve, or test source/destination
    precision metadata. Do not add precision sweeps or alternate-precision
    interop tests.
I9. A single fresh-process default-construction smoke test is sufficient to
    guard the 512-bit invariant. Do not add runtime default-mutation,
    cross-DSO, or worker-thread synchronization machinery unless P0 proves an
    existing MPLAPACK dependency and the maintainer approves it.
I10. Keep `GMPFRXX_MKII_ENABLE_FMA`, `GMPFRXX_MKII_FAST_FIXED_PREC`, and
     `GMPFRXX_MKII_FAST_STABLE_RND` identical across all TUs and exported build
     targets. Treat a mismatch as an ODR/configuration bug.
I11. Code, comments, tests, scripts, reports, and commit messages are in English.
     All MPLAPACK phases use `topic/gmpfrxx_mkII_migration`; one phase ends at
     one immutable candidate SHA and one report, while the entire MPLAPACK
     migration is one branch and one final PR. P2 upstream branches are governed
     by their own prompts. Do not reformat or clean unrelated code.
I12. Every acceptance statement is implemented by the committed phase gate and
     must fail closed. Run the gate and commit its complete output.
I13. Do not edit the frozen P0 baseline or ABI manifests in P3/P4. Unexpected
     changes produce a candidate-delta report and stop for maintainer review.
I14. Autotools may build committed externals. MPLAPACK CMake must not download
     or build precision dependencies; it consumes installed package targets.
```

---

## P0 — Revision lock, inventory, baseline, ABI freeze, compile spike

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: exact verified `origin/master` SHA selected at P0 start

Create the branch once from that full SHA. If the branch already exists locally
or remotely, do not reset or overwrite it; verify that it is an acceptable P0
starting point or stop and report.

This phase may change only `docs/migration/gmpfrxx_mkII/`. Production source,
public headers, generators, build files, and CI are read-only.

### Task 0 — Freeze repository inputs

Create `LOCK.json` with at least:

```json
{
  "schema": 1,
  "mplapack": {
    "repository": "https://github.com/nakatamaho/mplapack",
    "base_sha": "<40-hex SHA>"
  },
  "gmpfrxx_mkII": {
    "repository": "https://github.com/nakatamaho/gmpfrxx_mkII",
    "adapter_base_sha": "<40-hex SHA>",
    "bootstrap_release": "1.0.1",
    "bootstrap_archive": "gmpfrxx_mkII.1.0.1.tar.xz",
    "bootstrap_sha256": "c0816b3538b6b77009f714bb391cebe11abb2fdb69e07aa3bb305ff822764afb"
  },
  "precision_policy": {
    "default_bits": 512,
    "assignment_precision_policy": "unchecked",
    "interop_precision_parameters": "none"
  }
}
```

Replace the two SHA placeholders with verified full SHAs. Record why those
revisions were selected. Record the SHA-256 of the prototype interop document
and import it verbatim as `INTEROP_PROTOTYPE.md`; do not silently reconcile it
with repository reality.

Create `PRECONDITIONS.md` that records:

- the exact Autotools and CMake layouts at the locked MPLAPACK SHA;
- the generated-file policy from `AGENTS.md`;
- the existing gmpfrxx_mkII header roles, package targets, compile-time semantic
  options, current adapter file names, and ordinary default-precision APIs;
- where the 512-bit defaults are established in gmpfrxx_mkII and MPLAPACK;
- one minimal fresh-process probe showing that default-constructed
  `mpfrxx::mpfr_class` and `gmpxx::mpf_class` report 512 binary bits before any
  setter. Do not add alternate-precision, mutation, thread, or cross-DSO probes;
- any actual current MPLAPACK call that sets a non-512 process-wide default and
  requires that mutation to propagate across library images. Record such a
  dependency as a blocker; do not pre-implement provider machinery;
- any mismatch between this prompt and either locked repository. A material
  mismatch is a phase blocker, not an invitation to improvise.

### Task 1 — Reproducible inventory

Create `tools/inventory.py`. It emits both `INVENTORY.md` and
`inventory.json` from the source tree. Every hit records path, line, category,
generated/non-generated status, and owning build targets.

Inventory all of the following:

- `mpreal`, `mpcomplex`, `mpreal.h`, `mpcomplex.h`, `using namespace mpfr`;
- public declarations that hard-code legacy types, including MPBLAS and
  MPLAPACK backend headers; do not assume the switch is typedef-only;
- `gmpxx.h`, global `::mpf_class`, legacy global `::mpc_class`, and
  `include/mpc_class.h`;
- `-lgmpxx`, `libgmpxx`, GMP `--enable-cxx`, pkg-config metadata, installed
  CMake targets, packaging, release scripts, and CI;
- existing precision/default calls, including setters. Classify them as
  `sets-512`, `sets-non-512`, `query-only`, or `per-object`. This inventory is
  for compatibility review, not permission to add precision negotiation to
  interop;
- all comparison-side cast/conversion helpers involving MPFR/MPC and double,
  dd, qd, binary80, binary128, GMP MPF/MPFC, and complex counterparts;
- all reverse `cast2*` uses. Record whether each is genuinely used by MPLAPACK;
- all utility/math names needed by generated sources, tests, examples,
  benchmarks, and print paths;
- Fable generators/templates and generated output files that emit backend type
  names or signatures;
- all test executables and numerical metrics for MPFR and GMP backends;
- exported symbols, SONAME/install-name data, installed headers, pkg-config
  fields, CMake package targets/properties, and consumer examples for every
  backend.

Inventory commands must be reproducible and fail on unreadable files. Commit a
small fixture test proving that generated files, build metadata, and history
whitelist categories are distinguished correctly.

### Task 2 — Normative interoperability matrix

Create `interop_requirements.tsv` with columns:

```text
id  source_type  target_type  form  requirement  evidence  test_id
```

Allowed `requirement` values are `REQUIRED`, `FORBIDDEN`, `OUT_OF_SCOPE`, and
`NOT_USED`. Every REQUIRED row cites at least one current MPLAPACK call site or
comparison pattern and names a test.

The default REQUIRED direction is only:

- `double` -> `mpfrxx::mpfr_class`;
- dd real -> `mpfrxx::mpfr_class`;
- qd real -> `mpfrxx::mpfr_class`;
- binary80 real -> `mpfrxx::mpfr_class`;
- binary128 real -> `mpfrxx::mpfr_class`;
- `gmpxx::mpf_class` -> `mpfrxx::mpfr_class`;
- corresponding finite complex source types -> `mpfrxx::mpc_class`.

Use the exact concrete type names discovered by Task 1. The permitted form is
explicit construction, explicit assignment, or one named comparison helper,
whichever matches the smallest existing call-site change. The matrix does not
require all three forms.

Mark the following FORBIDDEN unless the maintainer explicitly changes the
matrix:

- MPFR/MPC -> GMP/MPFC or another backend;
- `mpfrxx::mpfr_class -> gmpxx::mpf_class`;
- `mpfrxx::mpc_class -> gmpxx::mpfc_class`;
- cross-family mixed arithmetic;
- adapter arithmetic in both operand orders;
- compound assignment with adapter operands;
- interop-driven expression-template composition.

Mark edd/td rows OUT_OF_SCOPE. Mark unused prototype operations NOT_USED; do not
implement them because they happen to exist in the prototype.

All REQUIRED rows share this binding policy:

- default-construct the MPFR/MPC comparison target and use ordinary
  construction/assignment semantics; the normal target has 512-bit precision;
- do not query or compare source/destination precision metadata;
- do not pass a precision argument solely for interop;
- do not negotiate, preserve, or sweep precision;
- no correctly-rounded conversion, exact round-trip, ULP, or adversarial
  component-expansion guarantee;
- no binary64 fallback for extended-precision sources.

### Task 3 — Numerical baseline freeze

Build the locked tree through the repository-supported Autotools path and run
the complete discovered MPFR and GMP test sets. Configure the ordinary 512-bit
default exactly as the locked project does. Do not introduce alternate-
precision runs solely for this migration.

Store raw logs and a normalized `baseline.json`. Each result records test ID,
metric, unit, direction, value, status, and source line.

Create `baseline_rules.tsv` with only metric-specific modes such as:

```text
exact
status-equal
upper-bound
lower-bound
nonincreasing-error
```

There is no generic `allowed-improvement`. Freeze units and directions in P0.
Create `tools/compare_baseline.py` and unit tests for identical input, worsening,
missing/duplicate/unknown tests, malformed logs, and an allowed improvement
under a predeclared directional rule. Later phases may not edit these files.

Record CMake's actual test level separately. Do not represent a smoke subset as
a full numerical baseline.

### Task 4 — ABI and install-manifest freeze

For every backend, capture normalized:

- exported mangled and demangled symbols;
- SONAME or platform install-name/version data;
- installed public headers and checksums;
- pkg-config fields;
- exported CMake target names, link interfaces, compile definitions, include
  directories, and package dependencies;
- a downstream compile/link/run consumer for Autotools/pkg-config and CMake.

Create `tools/compare_abi.py`. It permits a declared change only for the MPFR
backend in P3 and GMP backend in P4. Every non-migrated backend must remain
identical.

### Task 5 — Compile spike

Under `tools/spike/`, compile minimal translation units against the locked
MPLAPACK headers and locked gmpfrxx_mkII installation. Do not change production
code.

Exercise:

- ordinary same-family REAL/COMPLEX arithmetic used by generated MPLAPACK code;
- MPLAPACK utility, math, print, and hexadecimal-print surfaces;
- every REQUIRED one-way comparison embedding from Task 2;
- default construction at the ordinary 512-bit default;
- name collisions while legacy and replacement complex types coexist;
- expression-template lifetime hazards in non-interop arithmetic.

Do not add source/destination precision comparisons to the spike. A conversion
compiles and preserves a representative value sufficiently for the current
comparison code, or it does not.

Write exact missing signatures and semantic gaps into `SPIKE.md`. Distinguish
wrapper-family functionality needed for MPLAPACK REAL/COMPLEX from the narrow
comparison embeddings; the former does not justify broad adapter arithmetic.

### Acceptance

`tools/gate-P0.sh` must prove:

- only the migration documentation subtree changed;
- repository SHAs and archive hash are complete and valid;
- inventory generation and its fixture test pass;
- all REQUIRED matrix rows have evidence and tests;
- forbidden reverse/mixed operations have no REQUIRED row;
- the one fresh-process 512-bit default smoke passes for MPFR and MPF;
- no assignment/conversion precision-metadata test was added;
- baseline comparator unit tests and self-comparison pass;
- ABI/install manifests are complete for all backends;
- compile spikes run and `SPIKE.md` maps every failure to a later owner.

Deliverables: `LOCK.json`, `PRECONDITIONS.md`, inventory outputs,
`interop_requirements.tsv`, baseline data/rules/comparator, ABI manifests,
spike sources/results, gate, and `REPORT-P0.md`.

---

## P1 — Vendor the bootstrap gmpfrxx_mkII release

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: approved P0 SHA on the same migration branch

Verify that branch HEAD equals the approved P0 SHA; do not create another branch.

Vendor the exact bootstrap release from `LOCK.json` only to prove package
plumbing. No MPLAPACK source or public header may include gmpfrxx_mkII yet.

Create:

```text
external/gmpfrxx_mkII/
  Makefile.am
  download/gmpfrxx_mkII.1.0.1.tar.xz
  patches/
```

The patches directory must be empty and checked by the gate.

### Autotools external-package requirements

Follow the locked repository's current external-package conventions, including
`external/distfiles.sha256`, archive verification, extraction stamps, separate
`work/internal` and `work/install` trees, clean/distclean behavior, and
`make dist` inclusion.

Configure both gmpfrxx_mkII trees with CMake, but force:

```text
GMPFRXX_MKII_DEPS_AUTO_FETCH=OFF
GMPFRXX_MKII_BUILD_EXAMPLES=OFF
GMPFRXX_MKII_BUILD_BENCHMARKS=OFF
GMPFRXX_MKII_COMPONENTS=GMP,MPFR,MPC
```

Use the exact upstream option spelling discovered at the locked SHA. Point GMP,
MPFR, and MPC discovery at `external/i/{GMP,MPFR,MPC}` for bundled builds. A
missing dependency must fail; it must never trigger an upstream download.

Install `work/internal` into `external/i/GMPFRXX_MKII/` for the MPLAPACK build.
Install `work/install` into `$(prefix)` from MPLAPACK `make install`. This
intentional two-prefix installation applies to the Autotools bundled flow.

Install the headers and all ordinary exported CMake package targets required by
the locked release. Do not require or add a default-context provider artifact
solely for this migration.

Add bundled/system selection using the closest existing MPLAPACK option style.
A system prefix must be validated by compiling and linking a consumer against
all required installed targets. Do not accept a header-only path probe.

### MPLAPACK CMake requirements

Preserve the existing CMake architecture: it consumes installed dependencies
and does not build externals. Add package discovery and imported-target plumbing
for future phases, but do not link a backend target yet.

Test discovery in two modes:

1. staged-prefix mode: `CMAKE_PREFIX_PATH` includes the already populated
   `external/i/GMPFRXX_MKII` and matching GMP/MPFR/MPC prefixes;
2. independent-prefix mode: a clean installation of the bootstrap release.

There is no CMake "bundled build" mode in this phase.

### Acceptance

`tools/gate-P1.sh` must prove:

- archive SHA-256 matches `LOCK.json` and `external/distfiles.sha256`;
- corruption makes archive verification fail, followed by restoration and a
  clean successful verification;
- no network access is required after the committed archives exist;
- internal and final-prefix installations contain expected headers, package
  config, and exported targets;
- the two install manifests contain no build-tree absolute paths;
- CMake discovers the package in both staged and independent prefix modes;
- all P0 MPFR/GMP baseline builds remain green;
- no MPLAPACK production source includes a new wrapper header;
- `patches/` is empty;
- `make dist` contains the exact archive.

Deliverables: external package, build-system plumbing, gate, and
`REPORT-P1.md`.

---

## P2 — Complete the narrow comparison embeddings and wrapper surface

Repository: gmpfrxx_mkII
Branches: one branch/PR per P2 subphase

At the start of P2A, copy `interop_requirements.tsv`, `SPIKE.md`, and the
relevant `LOCK.json` entries into
`docs/mplapack_migration/requirements/`. Record SHA-256 values and never edit
those imported copies. A requirements change needs a new MPLAPACK P0 review.

P2 has two distinct responsibilities:

1. the narrow one-way embeddings used only by MPLAPACK comparison code;
2. ordinary same-wrapper-family functionality needed for
   `mpfrxx::mpfr_class`/`mpfrxx::mpc_class` and
   `gmpxx::mpf_class`/`gmpxx::mpfc_class` to serve as MPLAPACK REAL/COMPLEX.

Do not use responsibility 2 to justify general cross-backend adapter arithmetic.

### P2A — Real comparison embeddings

Branch: `topic/mplapack-compare-embed-real`

Implement only REQUIRED real rows from the imported matrix. Sources are the
exact discovered double, dd, qd, binary80, binary128, and
`gmpxx::mpf_class` types; the target is `mpfrxx::mpfr_class`.

Use the smallest explicit API compatible with current call sites: an existing
constructor/assignment, or one named helper such as
`to_mpfr_for_compare(value)`. Do not add several equivalent APIs.

Binding contract:

- the target is ordinarily default-constructed and therefore uses the project
  default of 512 binary bits;
- conversion and assignment do not inspect, compare, negotiate, or preserve
  source/destination precision metadata;
- no precision argument is required solely for interop;
- no exact-rounding, ULP, exact round-trip, or arbitrary non-canonical
  expansion guarantee is made;
- ordinary finite normalized values produced by MPLAPACK backends are in scope;
- no extended source may be routed through binary64.

Implementation guidance:

- reconstruct dd and qd from their stored components directly in MPFR using a
  straightforward component sum; exact treatment of adversarial exponent gaps
  is not required;
- use the native MPFR import route for binary80 and binary128 when that source
  wrapper is enabled; absence of the native route must be a configure/build
  failure for that adapter, not a fallback through `double`;
- import `gmpxx::mpf_class` directly from its GMP MPF representation;
- materialize the owning `mpfrxx::mpfr_class` before subtraction or other
  comparison arithmetic.

Do not add reverse conversion, mixed `+ - * /`, both operand orders, compound
assignment, or adapter expression-template composition.

Tests:

- each REQUIRED source compiles through the selected explicit API;
- one ordinary representative value per source compares successfully through
  existing MPFR comparison arithmetic;
- anti-binary64-fallback sentinels retain information below 53 bits where the
  source supports it. For dd, a low component around `2^-100` is sufficient;
  this is only a loose smoke test, not a 100-bit exactness proof;
- do not query source or destination precision in conversion tests;
- do not sweep precision or rounding modes.

### P2B — Complex comparison embeddings

Branch: `topic/mplapack-compare-embed-complex`

Implement only REQUIRED finite complex-source -> `mpfrxx::mpc_class` rows.
Use the selected real embedding independently for real and imaginary
components, then materialize an owning MPC wrapper before comparison
arithmetic.

The P2A precision and non-goals apply verbatim. Do not add mixed complex
arithmetic, reverse conversion, compound assignment, operand-order overloads,
round-trip tests, asymmetric precision tests, or per-component precision
inspection solely for interop.

Tests compile every REQUIRED complex source, preserve both components on an
ordinary representative value, and include one anti-binary64 sentinel. Existing
MPLAPACK complex comparison tests remain the main acceptance criterion.

### P2C — Ordinary wrapper compatibility and 512-bit default

Branch: `topic/mplapack-wrapper-compat`

This phase is not adapter interop. Implement only the ordinary same-family
operators, constructors, math overloads, and wrapper APIs required by
`SPIKE.md` for the replacement types to serve as MPLAPACK REAL/COMPLEX.

Requirements:

1. Keep MPLAPACK-specific utilities downstream. Upstream may add generic math
   overloads and wrapper facilities, but not `mplapackint`, `castREAL_*`,
   `nint`, `iceil`, `cabs1`, or other project-specific names.
2. Preserve expression-template lifetime rules and raw scalar-leaf rejection.
3. High-precision GMP real math required by the SPIKE must not use `get_d()` or
   another binary64 detour.
4. The ordinary default for `mpfrxx::mpfr_class` and `gmpxx::mpf_class` is
   exactly 512 binary bits. Preserve the existing gmpfrxx_mkII mechanism that
   establishes that default; do not invent a second precision state.
5. Add one fresh-process default-construction smoke test for the two real
   wrapper families. Do not add default mutation, non-512 sweeps, worker-thread
   initialization, cross-DSO synchronization, or existing-object precision
   tests solely for this migration.
6. Do not add a default-context provider unless the imported P0 evidence proves
   an existing MPLAPACK requirement and records explicit maintainer approval.

Add one 512-bit no-binary64-fallback regression based on an actual MPLAPACK
GMP-side operation identified by the SPIKE. It proves that high-precision math
does not plateau near binary64; it is not an interop or correct-rounding test.

### P2D — Hardening, documentation, and official release

Branch: `topic/mplapack-adapter-release`

1. Run the full upstream test suite in Debug and Release with dependency
   auto-fetch OFF.
2. Keep compile-fail tests for rejected raw scalar leaves and forbidden
   cross-family mixed arithmetic.
3. Add a lightweight expression-lifetime scanner for MPLAPACK to use in P3/P4;
   it is advisory but true positives must be fixed.
4. Document the exact one-way comparison API, its finite-value scope, the
   absence of precision metadata handling, the fixed 512-bit default
   assumption, the prohibition on binary64 fallback, and all explicit
   non-goals.
5. Document the ordinary wrapper API surface separately from interop.
6. Produce an official release archive using the upstream release machinery;
   record repository SHA, version, archive name, SHA-256, dependency versions,
   CMake options, exported targets, and license.
7. In a clean directory, install the archive with auto-fetch OFF and compile a
   consumer that exercises all REQUIRED embeddings, ordinary wrapper SPIKE
   calls, and the one fresh-process 512-bit smoke.

Do not add round-trip suites, source/destination precision checks, exact-
rounding matrices, or general mixed-adapter arithmetic.

The single gate is `tools/gate-P2D.sh`. Its report contains the final API/test
manifest consumed by P1R, P3, and P4.

---

## P1R — Re-vendor the P2D release

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P1 SHA on the same migration branch plus accepted P0 documentation

Verify that branch HEAD equals the accepted P1 SHA; do not create another branch.

Replace the bootstrap archive with the exact P2D official archive. Update
`LOCK.json`, `external/distfiles.sha256`, and the external package version in a
single commit. Record both upstream source SHA and archive SHA-256.

No source, generator, public header, or backend target change is allowed.

`tools/gate-P1R.sh` reruns P1 gates verbatim, verifies the new archive against
the P2D report, and proves that no bootstrap-version artifact remains outside
history documents.

---

## P2E — Generator readiness before type switching

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P1R SHA on the same migration branch

Verify that branch HEAD equals the accepted P1R SHA; do not create another branch.

This phase changes Fable/header-generation templates and their tests, not
backend semantics.

1. Generate `generated_files.json`, the authoritative manifest of generated
   routine bodies, public declarations, signature tables, and name maps.
2. Parameterize backend public-type emission so the generator can emit:
   - legacy MPFR: `mpreal` / `mpcomplex`;
   - new MPFR: `mpfrxx::mpfr_class` / `mpfrxx::mpc_class`;
   - legacy GMP: global `::mpf_class` / legacy `::mpc_class`;
   - new GMP: `gmpxx::mpf_class` / `gmpxx::mpfc_class`.
3. Preserve routine names and parameter passing exactly except for the approved
   backend type substitutions.
4. Add generator golden tests for all four modes.
5. Run the complete legacy regeneration pipeline. Production output in legacy
   mode must be byte-for-byte unchanged. A changed production generated file
   is a blocker in P2E.
6. Document the exact commands P3 and P4 must use to regenerate only the
   affected backend outputs.

Do not hand-edit generated files. Do not switch a production backend in this
phase.

`tools/gate-P2E.sh` verifies the generated-file manifest, golden tests, full
legacy regeneration byte identity, no production backend switch, and both
build systems at P0 gate levels.

---

## P3 — Switch the MPFR/MPC backend

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P2E SHA on the same migration branch

Verify that branch HEAD equals the accepted P2E SHA; do not create another branch.

Switch the MPFR backend from mpfrc++ to gmpfrxx_mkII. The legacy wrapper may
remain in-tree only for later deletion; no MPFR target may include or link it.

### Required work

1. Use the P2E generator mode to regenerate MPFR public declarations with
   `mpfrxx::mpfr_class` and `mpfrxx::mpc_class`. Never hand-edit those
   declarations.
2. Change the non-generated MPFR backend binding/include layer to the strict
   upstream header layering and only the adapter headers actually required by
   the P0 matrix.
3. Implement MPLAPACK-owned utility shims over generic gmpfrxx_mkII APIs,
   preserving signatures found in P0. Keep `sign`, `nint`, `iceil`, `cabs1`,
   `pow2`, `pow4`, `castREAL_*`, and `castINTEGER_*` downstream.
4. Port print and hexadecimal-print paths without binary64 detours.
5. Rewrite comparison/test call sites so dd/qd/double/binary80/binary128 results
   are explicitly materialized into MPFR/MPC before subtraction, norm, or
   tolerance evaluation. Do not preserve old syntax by adding mixed adapter
   arithmetic.
6. Rely on the established 512-bit wrapper default. Mechanically map an existing
   required 512-bit initialization call if the locked tree has one, but do not
   add precision arguments, source/destination precision queries, precision
   negotiation, precision sweeps, or per-assignment precision assertions.
7. Qualify every `mpfrxx::mpc_class` use while legacy global `::mpc_class`
   still exists.
8. Run the P2D expression-lifetime scanner over all non-generated MPFR code and
   fix every true positive by explicit materialization.
9. Autotools: use `external/i/GMPFRXX_MKII` for in-tree builds and install
   public dependency metadata correctly.
10. CMake: consume installed MPFR/MPC package targets with `find_package`; do
    not build the external archive. Export dependencies transitively to
    installed MPLAPACK MPFR consumers.
11. Propagate locked semantic compile definitions consistently.

### Gates

`tools/gate-P3.sh` must prove:

- generated changes came from the P2E command and match its expected diff
  class;
- full Autotools MPFR tests pass;
- CMake MPFR build/tests at the P0 CMake level pass in staged-prefix and
  independent-prefix modes;
- P0 numerical comparison has no unclassified delta and no baseline file was
  edited;
- MPFR ABI change is recorded, while every non-migrated backend ABI/install
  manifest is identical to P0;
- installed Autotools and CMake consumer programs compile and run;
- no MPFR target/source/public header uses `mpreal`, `mpcomplex`, `mpreal.h`,
  `mpcomplex.h`, `gmpxx.h`, `-lgmpxx`, or an unqualified legacy `mpc_class`;
- the one fresh-process default-construction smoke reports 512 bits for
  `mpfrxx::mpfr_class`;
- all REQUIRED one-way embeddings compile and existing MPLAPACK comparison
  tests pass without reverse conversion, mixed adapter arithmetic, or
  assignment precision checks;
- no hidden network fetch occurs.

An unexpected baseline delta writes candidate files and stops for maintainer
approval; the agent must not relax a rule.

---

## P4 — Switch the GMP/MPFC backend

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P3 SHA on the same migration branch

Verify that branch HEAD equals the accepted P3 SHA; do not create another branch.

Switch the GMP backend from libgmpxx and legacy `include/mpc_class.h` to
`gmpxx::mpf_class` and `gmpxx::mpfc_class`.

### Required work

1. Use the P2E generator mode for GMP public declarations. Do not hand-edit
   generated headers.
2. Replace `gmpxx.h` and legacy complex includes with `gmpxx_mkII.h` and the
   exact required headers.
3. In GMP comparison/test drivers, explicitly materialize
   `gmpxx::mpf_class` results into `mpfrxx::mpfr_class` and
   `gmpxx::mpfc_class` results into `mpfrxx::mpc_class` before error
   arithmetic. Do not add MPFR -> GMP conversion or GMP/MPFR mixed arithmetic.
4. Implement MPLAPACK-owned utility shims over generic upstream APIs and remove
   legacy GMP transcendental workarounds that round through binary64.
5. Preserve documented MPF/MPFC semantics; do not advertise MPFR/MPC correct
   rounding.
6. Rely on the established 512-bit default for `gmpxx::mpf_class`. Do not add a
   default-context provider, cross-image synchronization, runtime default
   mutation test, non-512 sweep, or per-assignment precision check solely for
   this migration. If P0 recorded an actual required cross-image mutation,
   stop until the maintainer supplies a separate approved design.
7. Run one 512-bit no-binary64-fallback regression. It is a backend math test,
   not an interop exactness test.
8. Remove all MPLAPACK links and package requirements for libgmpxx. If P0 finds
   no remaining third-party dependency on GMP C++, remove `--enable-cxx` from
   the bundled GMP configuration. If a real dependency exists, stop for a
   maintainer decision rather than preserving it speculatively.
9. Run the expression-lifetime scanner and fix all true positives.

### Baseline policy

The historical GMP implementation may have binary64-limited results. P4 may
produce genuine numerical improvements, but it may not edit P0 rules. The gate
accepts only improvements already expressed by a frozen metric rule. Any other
changed result goes to `baseline_delta_candidates.tsv` with old/new values and
an independent 512-bit MPFR reference, then the phase stops for maintainer
review. Accuracy regressions are always failures.

### Gates

`tools/gate-P4.sh` must prove:

- full Autotools GMP tests pass;
- CMake GMP build/tests pass in staged-prefix and independent-prefix modes;
- the one fresh-process default-construction smoke reports 512 bits for
  `gmpxx::mpf_class`;
- no precision metadata is queried by the comparison embeddings and no
  provider/mutation machinery was introduced without an approved P0 blocker;
- no unclassified numerical delta and no P0 baseline edit;
- GMP ABI change is recorded, while every non-migrated backend ABI/install
  manifest remains identical;
- no MPLAPACK source, public header, build file, installed metadata, or target
  includes/links/requires `gmpxx.h` or libgmpxx;
- no legacy unqualified `mpc_class` remains in the GMP path;
- no binary64 precision plateau;
- required GMP-result -> MPFR/MPC comparison embeddings pass existing tests,
  with no MPFR/MPC -> GMP/MPFC path;
- installed consumer programs compile and run through both build systems.

---

## P5 — Remove mpfrc++ and legacy wrapper infrastructure

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P4 SHA on the same migration branch

Verify that branch HEAD equals the accepted P4 SHA; do not create another branch.

1. Delete `mpfrc++/` and `include/mpc_class.h`.
2. Remove every include path, install rule, distribution entry, package record,
   license aggregation entry, build option, and CI reference identified by the
   current inventory.
3. Re-run the inventory against the whole source tree, generated outputs,
   packaging, release scripts, and installed trees.
4. Inspect the exact bundled mpfrc++ license files and header notices before
   editing license documentation. Do not assume LGPL: record the exact license
   text/version of the removed copy and remove only its corresponding notices.
   Record the exact gmpfrxx_mkII license from the vendored release.
5. Add `tools/purge_gate.py` with a committed structured whitelist for history
   documents. It rejects, outside the whitelist:
   - `gmpxx.h` includes or libgmpxx linkage;
   - `mpreal`, `mpcomplex`, mpfrc++, and legacy `mpc_class.h`;
   - stale installed target properties, pkg-config flags, archive contents, and
     generated outputs.
6. Finalize the downstream migration table: type names, headers, one-way
   comparison embeddings, fixed 512-bit default assumption, unchecked
   assignment precision policy, expression-template lifetime example, and
   MPF/MPFC limitations. Do not document a provider requirement.

`tools/gate-P5.sh` runs full MPFR/GMP Autotools suites, the CMake gate level,
source-tree purge, installed-prefix purge for both build systems, archive
content purge, ABI checks for non-migrated backends, and consumer tests.

---

## P6 — Full code-generation and residual-reference verification

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P5 SHA on the same migration branch

Verify that branch HEAD equals the accepted P5 SHA; do not create another branch.

P2E made generators switchable; P3/P4 used the new modes. P6 is now a complete
regeneration and residual audit, not the first time codegen learns the new
types.

1. Run the full supported Fable regeneration pipeline for library routines,
   public headers, signatures, name maps, and generated tests.
2. The resulting diff must be empty. Any non-empty diff means P3/P4 failed to
   commit a generator-derived output or another generator remains stale.
3. Prove that no generator/template can emit legacy wrapper names in an active
   backend mode.
4. Extend purge coverage to all Fable scripts and archived generated fixtures,
   with only explicit historical-document exceptions.
5. Re-run both backend suites and non-migrated ABI checks.

`tools/gate-P6.sh` fails on any regeneration diff or active legacy token.

---

## P7 — CI and supported-platform matrix

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P6 SHA on the same migration branch

Verify that branch HEAD equals the accepted P6 SHA; do not create another branch.

Update committed CI and release buildtest scripts without weakening existing
non-migrated backend coverage.

Required coverage, pruned with comments rather than a Cartesian explosion:

- Linux, macOS, and MinGW/Wine paths already supported by the locked tree;
- Autotools bundled gmpfrxx_mkII archive flow;
- CMake `find_package` flow against a clean installed P2D release;
- at least one CMake staged-prefix job using the Autotools internal prefix;
- MPFR and GMP backends in both build systems;
- one fresh-process 512-bit default smoke for MPFR and MPF in the primary
  Linux job; do not multiply precision checks across the matrix;
- binary80 only on supported x86/x86_64 runners;
- a compile-only job that deliberately changes one semantic definition in one
  TU and verifies the ODR/configuration guard fails;
- purge gate, baseline comparison, ABI comparison, archive verification,
  installed consumer tests, and no-network build checks;
- all pre-existing dd/qd/double/binary80/binary128 coverage retained.

Do not call a CMake job "bundled" when CMake merely consumes a staged prefix.
Do not let upstream gmpfrxx_mkII auto-fetch private GMP/MPFR/MPC copies. Do not
add provider-linkage or precision-mutation matrix axes.

`tools/gate-P7.sh` validates workflow structure locally where possible and the
report records every remote workflow result and any intentionally pruned
combination.

---

## P8 — MPLAPACK 3.0.0 release engineering

Repository: MPLAPACK
Branch: `topic/gmpfrxx_mkII_migration`
Base: accepted P7 SHA on the same migration branch

Verify that branch HEAD equals the accepted P7 SHA; do not create another branch.

1. Write `CHANGES.3.0.0.md`: MPFR/GMP ABI break, replacement types, removal of
   mpfrc++/libgmpxx use, exact vendored release/hash, one-way comparison
   embedding policy, 512-bit default assumption, unchecked assignment
   precision policy, no-binary64-fallback fixes, MPF/MPFC limitations, and
   edd/td non-scope.
2. Update package versions in both `configure.ac` and `CMakeLists.txt`.
3. Compute libtool `current:revision:age` independently for each affected
   installed library based on the actual locked layout. Show symbol and
   interface reasoning. Do not equate package version with libtool version.
4. Finalize `MIGRATION.md` and installed documentation.
5. Run `make distcheck` with the release configuration and committed external
   archive.
6. Clean-room Autotools release-tarball tests must build and run both MPFR and
   GMP smoke/full release-required subsets using bundled externals and no
   network.
7. Clean-room CMake source-tarball tests must use installed dependencies, build
   MPFR and GMP, install MPLAPACK, and compile/run downstream `find_package`
   consumers. CMake must not build the committed external archive.
8. Verify release archive contents, licenses, install manifests, package
   metadata, the one fresh-process 512-bit smoke, and exact gmpfrxx_mkII
   SHA-256. No provider artifact is required by this migration.
9. Produce a maintainer release checklist. The agent does not tag or publish a
   GitHub Release.

`tools/gate-P8.sh` is the final release gate and includes all P5 purge checks,
P0 baseline/ABI checks, `make distcheck`, both clean-room paths, and version
consistency checks.

---

## Appendix A — Normative public type and comparison mapping

| Removed API/use | Replacement |
|---|---|
| `mpfr::mpreal` | `mpfrxx::mpfr_class` |
| `mpfr::mpcomplex` | `mpfrxx::mpc_class` |
| global `::mpf_class` from libgmpxx | `gmpxx::mpf_class` |
| legacy global `::mpc_class` | `gmpxx::mpfc_class` |
| backend real result compared with MPFR reference | explicit construction/assignment or the single approved `to_mpfr_for_compare`-style helper, then MPFR arithmetic |
| backend complex result compared with MPC reference | explicit construction/assignment or the single approved `to_mpc_for_compare`-style helper, then MPC/MPFR arithmetic |
| legacy MPFR/MPC -> dd/qd/binary/GMP `cast2*` use | removed; no replacement in this migration |
| raw `long double` adapter leaf | MPLAPACK binary80 wrapper |
| raw `__float128` adapter leaf | MPLAPACK binary128 wrapper |

The exact accepted helper names come from the P2D release report. The names
above describe roles, not permission to create parallel APIs.

## Appendix B — Narrow comparison-embedding contract

1. Interop is one-way: double/dd/qd/binary80/binary128/GMP results are
   explicitly materialized into MPFR/MPC reference objects.
2. The ordinary target is default-constructed at 512 binary bits. Conversion
   and assignment do not query, compare, negotiate, preserve, or test source
   and destination precision metadata.
3. dd/qd import sums stored components directly in the MPFR target. Ordinary
   finite normalized backend results are in scope. Exact adversarial expansion
   handling is not required.
4. binary80, binary128, and GMP MPF use native direct MPFR import paths.
   Extended sources must never be routed through binary64.
5. Conversion need not be correctly rounded and is not tested by ULP sweeps,
   precision sweeps, precision-metadata comparison, retained-bit accounting, or
   round-trip identity. A loose anti-binary64 sentinel plus existing MPLAPACK
   comparison tests is sufficient; for dd, retaining a term around `2^-100`
   satisfies the intended approximately-100-bit comparison use.
6. After materialization, subtraction, norms, tolerances, and decisions are
   performed entirely in MPFR/MPC.
7. MPFR/MPC -> source conversion, MPFR -> MPF, MPC -> MPFC, source-type/MPFR
   mixed arithmetic, both-operand-order overloads, and compound assignment are
   outside scope.
8. MPF/MPFC ordinary backend math remains arbitrary-precision and free of
   binary64 fallback but does not claim MPFR/MPC correct rounding.
9. Baseline rules are frozen in P0. Later agents report unexpected changes and
   stop; they do not edit the rules.

## Appendix C — 512-bit default and unchecked-assignment policy

1. `512` means binary precision bits, not decimal digits.
2. gmpfrxx_mkII and MPLAPACK both use 512 bits as the normal default for
   `mpfrxx::mpfr_class` and `gmpxx::mpf_class`. Their complex counterparts use
   the ordinary wrapper-family component defaults.
3. A single fresh-process smoke test guards this invariant. The migration does
   not add alternate-precision, mutation, worker-thread, cross-DSO, or
   existing-object precision tests.
4. Interop uses ordinary default construction and assignment. It does not
   inspect precision metadata or pass precision parameters solely for
   conversion.
5. Existing explicit per-object precision APIs are not removed, but their
   non-512 behavior is outside this migration's acceptance contract.
6. A process-wide provider is not introduced unless P0 proves an existing
   MPLAPACK dependency on cross-image runtime default mutation and the
   maintainer approves a separate design.

## Appendix D — Build-model distinction

| Path | Dependency behavior |
|---|---|
| MPLAPACK Autotools | Builds committed external archives; stages internal prefix; installs final prefix |
| MPLAPACK CMake | Uses `find_package`; never downloads/builds precision dependencies |
| gmpfrxx_mkII external build | Uses its CMake project with dependency auto-fetch forced OFF |
