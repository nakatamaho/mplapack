# P0 Report

Status: ACCEPTED P0 candidate; gate-P0 PASS.

## Repository and lock verification

- Repository root: `/home/docker/mplapack`
- `origin`: `git@github.com:nakatamaho/mplapack.git` (`nakatamaho/mplapack`)
- Branch: `topic/gmpfrxx_mkII_migration`; existing local branch reused, not recreated, reset, rebased, renamed, or force-updated.
- Initial/current pre-commit HEAD: `b875e74d4b927282c907c3a29e6cadda48a7d57b`
- Verified `origin/master`: `b875e74d4b927282c907c3a29e6cadda48a7d57b`
- Initial merge base: the same SHA; initial ahead/behind: `0/0`
- Remote migration branch: absent at verification time.
- The only pre-existing supplied worktree change was the canonical v5.4 prompt pack; unrelated `external/gmp/work/` build output was removed after the maintainer explicitly allowed cleanup.

The P0 base is immutable. No newer `master` commit was merged.

## Prompt and historical prototype provenance

- Canonical prompt: `docs/migration/gmpfrxx_mkII/CODEX_PROMPTS.md`
- Prompt identity: v5.4, SHA-256 `c4175afa3ec59c385f00a54c7b762b09c0f05b7c8418b33e940ce340d46d2c0e`
- Historical source filename: `FABLE_PROTOTYPE_PROMPT.v4.md`
- Historical source path: `/home/docker/mplapack-gmpfrxx-v5.4.1-input/FABLE_PROTOTYPE_PROMPT.v4.md`
- Source SHA-256: `d14b03964672da96bb73c86960913d9a566a209bcba86379083bd71c56f262a4`
- Imported path: `docs/migration/gmpfrxx_mkII/INTEROP_PROTOTYPE.md`
- Imported SHA-256: `d14b03964672da96bb73c86960913d9a566a209bcba86379083bd71c56f262a4`
- Byte-for-byte identity: PASS, `cmp -s` exit 0.

The v4 document is historical provenance only. Canonical v5.4 and the maintainer decisions control this phase. Reverse conversion, mixed adapter arithmetic, precision sweeps, exact-rounding proofs, phase-specific branches, and edd/td support were not executed.

## Environment

The complete environment command output is in `logs/environment.log`:

- Linux x86_64, GCC/G++/GFortran 15.2.0, CMake 4.2.3
- Autoconf 2.72, Automake 1.18.1, Libtool 2.5.4
- `nproc` 64; requested compilation/test scheduling used `-j32`
- dd translation units retained `-ffp-contract=off`.

The locked dependency versions were GMP 6.3.0, MPFR 4.2.2, MPC 1.4.1, and QD 2.3.24. The locked gmpfrxx_mkII source is `2f06785c3f1b62f92e1e2026c2c975df73d1e426`.

## Build and install commands

Autotools was configured with:

```sh
env CC='ccache gcc' CXX='ccache g++' FC='ccache gfortran' ./configure \
  --prefix=/tmp/mplapack-p0-final --disable-dependency-tracking \
  --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes \
  --enable-qd=yes --enable-dd=yes --enable-double=yes \
  --enable-binary80=yes --enable-test=yes --enable-benchmark=no
```

The effective arguments are also recorded by `/tmp/mplapack-p0-base/config.status --config` and `logs/autotools-config.log`. The clean archive lacked `ltmain.sh`, so `gen_configure.sh` could not complete; the documented Git-checkout bootstrap `autoreconf --force --install` succeeded. This is recorded as a bootstrap deviation, not a source change.

```sh
make -j32
make -j32 install
```

Both commands exited 0. Complete output is in `logs/autotools-build.log` and `logs/autotools-install.log`. The earlier interrupted `-j16` attempt is retained as `logs/autotools-build-j16-partial.log` and is not the accepted build.

The locked gmpfrxx_mkII tree was configured with auto-fetch OFF, components GMP/MPFR/MPC, and include/library paths under MPLAPACK's internal dependency prefix. Its complete cache is `logs/gmpfrxx-cmake-cache.txt`; it was installed to `/tmp/gmpfrxx-p0-install`.

CMake was configured with all seven backends ON, OPT ON, OpenCL ON, CUDA OFF, and tests/examples/benchmarks OFF. It built and installed with `-j32`; complete output is in `logs/cmake-configure.log`, `logs/cmake-build.log`, and `logs/cmake-install.log`. The actual CMake test level is recorded separately in `CMAKE-TEST-LEVEL.md`; it is not represented as the numerical baseline.

## Task 1: inventory

`tools/inventory.py` deterministically generated `INVENTORY.md` and `inventory.json`. The inventory contains 42,831 hits, including generated/build/history classifications and owning targets. The fixture test passed. Migration output files are excluded from source scanning except the historical prototype, preventing self-referential reports and logs.

## Task 2: interoperability matrix

`interop_requirements.tsv` contains 12 REQUIRED one-way materialization rows: real and complex double, dd, qd, binary80, binary128, and GMP MPF/MPFC into MPFR/MPC. Every REQUIRED row has a current call-site evidence citation and spike test ID. Reverse conversions, mixed arithmetic, both operand orders, compound adapter assignment, and expression-template adapter composition are FORBIDDEN. edd and td are OUT_OF_SCOPE; unused historical prototype operations are NOT_USED.

## Task 3: numerical baseline

The accepted baseline is produced only after all three `-j32` test sessions finish:

```sh
python3 docs/migration/gmpfrxx_mkII/tools/capture_baseline.py \
  --build-root /tmp/mplapack-p0-base \
  --output docs/migration/gmpfrxx_mkII/baseline.json \
  --raw-output docs/migration/gmpfrxx_mkII/logs/baseline-raw.log.gz
```

The scope is the complete discovered Autotools MPBLAS, LIN, and EIG sets for GMP and the ordinary default MPFR profile at exactly 512 bits. Capture produced 874 normalized results from 696 raw files; all normalized statuses are PASS. Raw output is archived in `logs/baseline-raw.log.gz`; complete command logs are `logs/baseline-mpblas-gmp-check-all.log`, `logs/baseline-mpblas-mpfr-check-all.log`, `logs/baseline-lin-gmp.log`, `logs/baseline-lin-mpfr-default.log`, `logs/baseline-eig-gmp.log`, and `logs/baseline-eig-mpfr-default.log`.

`baseline_rules.tsv` uses exact/status-equal/upper-bound/nonincreasing-error metric-specific rules only. The parser records both normal threshold summaries and the locked EIG largest-error/info/total-example and DMD summary formats. It preserves arbitrarily large decimal exponents without rounding. Comparator and capture fixture tests cover identical, worsening, missing, duplicate, unknown, malformed, directional improvement, special EIG, DMD, and huge-exponent cases.

## Task 4: ABI and install manifest

`abi/manifest.json` covers all seven backends, 15 installed shared libraries, 64 public-header sets, mangled/demangled symbol files, SONAMEs, pkg-config fields, CMake targets/properties/dependencies, and 14 downstream compile/link/run consumers. The ABI comparator permits future changes only for MPFR in P3 and GMP in P4; P0 self-comparison passes.

An existing package metadata gap is recorded: Autotools `.pc` files do not expose all transitive GMP/MPFR/MPC/QD link dependencies. Consumers therefore use the dependency package closure explicitly in the P0 capture tool. This does not modify production metadata or hide the defect.

## Task 5: compile spike

The mandatory single fresh-process smoke default-constructs `mpfrxx::mpfr_class` and `gmpxx::mpf_class`; both report exactly 512 bits before any setter. Same-family arithmetic, utility/print/hex, name collision, expression lifetime, and double embeddings pass. Missing real/complex extended/GMP embeddings and the dd/qd complex binary64-fallback result are mapped explicitly to P2A/P2B/P2C/P3/P4 in `SPIKE.md`; no production source was changed.

## Gate

The committed gate is `tools/gate-P0.sh`. It was run exactly as:

```sh
bash docs/migration/gmpfrxx_mkII/tools/gate-P0.sh
```

The complete gate output is recorded in `logs/gate-P0.log`; result: PASS. It verifies the immutable base, docs-only changes, hashes and prototype identity, inventory determinism, interop matrix, the one 512-bit smoke, baseline and ABI comparator tests/self-comparisons, all spike results, and the separate CMake test-level record.

## Contradictions, deviations, and blockers

- Historical v4 and canonical v5.4 differ materially; all differences are recorded in `PRECONDITIONS.md` and were resolved in favor of v5.4 and maintainer decisions.
- The locked MPLAPACK tree contains legacy runtime precision setters and non-512 profiles. The maintainer explicitly accepted retiring those `mplapackinit` precision changes and profiles during migration; no provider or cross-DSO synchronization machinery was added.
- The clean archive bootstrap required `autoreconf --force --install` because `ltmain.sh` was absent.
- CMake numerical tests were OFF by explicit P0 test-level decision; CMake package consumers and all CMake build/install checks passed.
- Existing pkg-config transitive dependency metadata is incomplete; the gap is captured in the ABI manifest and report.
- CUDA was OFF and OpenCL was ON for the CMake package/build manifest; accelerator migration is not a P0 numerical requirement.
- Blockers after the maintainer precision decision: none.

## Push

Nothing was pushed. P0 stops after its commit; P1 and later phases were not started.
