# MPLAPACK -> gmpfrxx_mkII migration — Codex prompt pack v4

Canonical in-repo location: `docs/migration/gmpfrxx_mkII/CODEX_PROMPTS.md`.
Supersedes v3. Carry-over decisions from v3 remain binding:

- Target version: MPLAPACK 3.0.0. C++ ABI break is ALLOWED for the MPFR and
  GMP backends only. All other backends (dd, qd, double, binary80, binary128)
  keep their ABI.
- No stubs, no placeholder implementations, ever. If blocked, stop and report.
- Backend switch happens only after the interop surface (P2A–P2D) is complete
  and released upstream.
- Fable/codegen alignment is an independent phase (P6), never mixed into the
  type-switch phases.
- P0 produces a compile spike and a machine-readable baseline comparator
  (`baseline_rules.tsv` + `tools/compare_baseline.py`).

## What changed since v3

1. The interop adapters now live UPSTREAM in gmpfrxx_mkII as opt-in headers
   under `include/gmpfrxx_mkII/adapters/` (see
   `docs/migration/gmpfrxx_mkII/INTEROP_STATUS.md`, imported from the
   prototype write-up). Implemented today: construction, assignment,
   `operator-` in both operand orders, `cast2*` reverse conversion, and
   expression-operand use, for: dd/td/qd/edd real+complex, binary80/binary128
   real+complex, and GMP `mpf_class`/`mpfc_class` as targets. NOT yet
   implemented: `operator+`, `operator*`, `operator/`, compound assignment,
   and the math-function surface.
2. The "interop as patches under external/gmpfrxx_mkII/patches/" plan is
   dropped. Interop is upstream work in the gmpfrxx_mkII repository; MPLAPACK
   vendors release tarballs. `external/gmpfrxx_mkII/patches/` exists only as
   an emergency escape hatch and must normally be empty.
3. Scope restriction: edd and td adapters exist upstream but are OUT OF SCOPE
   for MPLAPACK. MPLAPACK has no edd/td backend. Do not wire them in, do not
   extend them, do not test them from the MPLAPACK side.
4. Vendoring model is settled: gmpfrxx_mkII release tarball under
   `external/gmpfrxx_mkII/`, following the existing external package pattern
   exactly (`distfiles.sha256`, `verify-distfile`, stamp files, `work/internal`
   installed into `external/i/GMPFRXX_MKII/` for use during the MPLAPACK
   build, `work/install` installed into `$(prefix)` at `make install`). Yes,
   this installs the headers twice; that is the established convention for
   every external package and is intentional.
5. MPLAPACK now has two build systems (autotools and CMake). Both must be
   wired and both must stay green in every phase.

## Repositories

- `nakatamaho/mplapack` — main repository. Phases P0, P1, P1R, P3–P8.
- `nakatamaho/gmpfrxx_mkII` — upstream wrapper library. Phases P2A–P2D only.

## Phase graph

```
P0 ──► P1 ──────────────┐
        │               ▼
        │   P2A ► P2B ► P2C ► P2D   (gmpfrxx_mkII repo)
        │               │
        │               ▼
        └────────────► P1R ► P3 ► P4 ► P5 ► P6 ► P7 ► P8
```

P1 may run in parallel with P2 (vendor the current release to prove the
plumbing; P1R re-vendors the interop release). P3 must not start before P1R.

## COMMON CONTEXT (prepend verbatim to every Codex session)

```text
You are working under docs/migration/gmpfrxx_mkII/CODEX_PROMPTS.md (v4).
Read AGENTS.md at the repository root first and obey its hard rules. In
particular for mplapack: never edit routine bodies under mpblas/reference/
or mplapack/reference/ (FABLE-generated); dd requires -ffp-contract=off on
every TU; version numbers live in BOTH configure.ac and CMakeLists.txt;
autotools AND CMake builds must both stay green.

Migration-wide invariants:

I1. No stubs, no TODO placeholders, no "temporarily return 0". If you cannot
    complete something, stop, write what blocks you into the phase report,
    and end the session.
I2. Expression-template contract: public arithmetic operators involving
    gmpfrxx_mkII types return expression nodes that reference their
    operands. Never bind such an expression to `auto` and let it outlive an
    operand; materialize into mpfr_class / mpc_class / mpf_class /
    mpfc_class explicitly.
I3. Raw `long double` and raw `__float128` are rejected as scalar leaves by
    design. The compile-fail tests enforcing this must keep passing. Use the
    explicit wrapper types binary80 / binary128.
I4. No silent double fallback: no conversion or math path for a
    multiprecision type may round through double or long double unless the
    SOURCE type is double or binary80 itself. This is a known historical bug
    class (GMP transcendentals falling back via get_d(), capping accuracy at
    2^-53); reintroducing it is a hard failure.
I5. edd and td are out of scope on the MPLAPACK side.
I6. Code, comments, and commit messages in English. One phase = one branch
    (topic/gmpfrxx-<phase>) = one PR.
I7. Every acceptance gate is a command. Run it, paste the output into
    docs/migration/gmpfrxx_mkII/REPORT-<phase>.md, commit the report.
I8. Do not reformat, rename, or "clean up" code outside the stated scope of
    the phase.
```

---

## P0 — Inventory, baseline freeze, compile spike (mplapack repo)

```text
[COMMON CONTEXT]

Phase P0. Branch topic/gmpfrxx-p0. This phase changes NOTHING outside
docs/migration/gmpfrxx_mkII/. Production sources, headers, and build files
are read-only.

Task 1 — INVENTORY.md
Produce docs/migration/gmpfrxx_mkII/INVENTORY.md: an exhaustive, per-file
inventory of every dependency on mpfrc++ and gmpxx.h:
  - the token `mpreal` and `mpcomplex` (declaration sites, using-directives
    such as `using namespace mpfr;` in include/mplapack_mpfr.h, direct uses
    in benchmarks/examples/test drivers/utils/print headers);
  - `#include <gmpxx.h>` / `#include "gmpxx.h"` and uses of global
    ::mpf_class;
  - include/mpc_class.h (MPLAPACK's own GMP complex class) and its users;
  - -lgmpxx and gmpxx references in configure.ac, Makefile.am tree, CMake
    tree, pkg-config templates, packaging/;
  - cast helpers between mpreal/mpcomplex and dd/qd/binary80/binary128
    types (these define the interop surface P2 must cover);
  - the exact typedef layer that binds REAL/COMPLEX/INTEGER to backend types
    for the mpfr and gmp backends (name the files and lines).
Categorize each hit: typedef-layer / utils / print / benchmark / example /
test-driver / build-system / codegen-template / doc. The grep commands must
be committed as docs/migration/gmpfrxx_mkII/tools/inventory.sh so the
inventory is reproducible.

Task 2 — Baseline freeze
Build current master with the autotools build for the mpfr and gmp backends
(bundled externals), run the full test suites for those two backends, and
store the logs under docs/migration/gmpfrxx_mkII/baseline/. Then:
  - commit tools/compare_baseline.py: reads a stored log and a new log,
    applies rules from baseline_rules.tsv, exits nonzero on any
    unclassified deviation;
  - commit baseline_rules.tsv with columns:
    test_id <TAB> metric <TAB> mode <TAB> parameter
    where mode is one of exact | ratio<= | abs<= | allowed-improvement.
    Initially every test is `exact` unless the log format forces a looser
    rule; document each looser rule inline with a # comment.
  - CMake build: configure + build the mpfr and gmp backend libraries and
    run a smoke subset; record which subset in BASELINE.md. Full-suite
    numerical baseline is taken from the autotools build only.
Write BASELINE.md: exact commit SHA, configure flags, compiler versions,
how to re-run everything.

Task 3 — Compile spike
Under docs/migration/gmpfrxx_mkII/tools/spike/, add small translation units
(NOT wired into any build system; a standalone Makefile in the same
directory is fine) that include mpfrxx_mkII.h / mpcxx_mkII.h /
gmpxx_mkII.h from a locally installed gmpfrxx_mkII together with the
MPLAPACK headers the mpfr and gmp typedef layers would need, and attempt:
REAL arithmetic as generated code performs it, utils calls (sign, nint,
iceil, castREAL_*, castINTEGER_*, pow2, pow4, min/max, mod, frac, pi,
cabs1), and printnum paths. Do NOT fix anything. Enumerate every compile
failure and semantic mismatch in SPIKE.md, grouped into:
  (a) missing operators (expected: +, *, /, compound assignment vs adapter
      types),
  (b) missing math/utils functions,
  (c) precision/rounding-control API differences (mpreal global default
      precision vs mpfrxx TLS default precision),
  (d) name collisions — in particular mpfrxx::mpc_class vs MPLAPACK's
      legacy global ::mpc_class (GMP complex), which coexist until P5,
  (e) anything else.
SPIKE.md is the requirements input for P2; be precise about signatures.

Acceptance:
  git diff --stat master..HEAD touches only docs/migration/gmpfrxx_mkII/;
  inventory.sh runs and reproduces INVENTORY.md counts;
  compare_baseline.py baseline/<log> baseline/<log> exits 0 (self-compare);
  SPIKE.md lists at least the operator and utils gaps with signatures.
Deliverables: INVENTORY.md, BASELINE.md, SPIKE.md, baseline logs,
baseline_rules.tsv, tools/, REPORT-P0.md.
```

---

## P1 — Vendor gmpfrxx_mkII under external/ (mplapack repo)

```text
[COMMON CONTEXT]

Phase P1. Branch topic/gmpfrxx-p1. Vendor gmpfrxx_mkII (current upstream
release, even though interop is incomplete — P1R re-vendors later) as a new
external package. NOTHING in mplapack sources may include it yet; this
phase is build plumbing only.

Model it on external/gmp/Makefile.am and external/mpfr/Makefile.am:

  external/gmpfrxx_mkII/
    Makefile.am
    download/gmpfrxx_mkII.<version>.tar.xz   (official dist tarball)
    patches/                                  (empty; escape hatch only)

Requirements:
  - SHA256 entry in external/distfiles.sha256; reuse the verify-distfile
    recipe pattern.
  - Stamp-file structure identical in spirit to the other externals
    (extract/patch/configure/build/check/install_internal, plus the final
    install tree).
  - Two extracted trees like the other packages: work/internal configured
    with CMAKE_INSTALL_PREFIX=$(abs_builddir)/../i/GMPFRXX_MKII/ and
    installed during the normal build; work/install configured with
    CMAKE_INSTALL_PREFIX=$(prefix) and installed at `make install`. The
    package is header-only (CMake INTERFACE targets), so "build" is
    configure + install of headers and CMake package files; keep the stamp
    sequence anyway for uniformity.
  - CMake configure arguments for both trees: point GMP/MPFR/MPC discovery
    at external/i/{GMP,MPFR,MPC} when MPLAPACK builds bundled externals,
    honor system libraries otherwise; disable the package's tests,
    examples, and benchmarks by default; add a single knob that controls
    MPFR_WANT_FLOAT128 consistency and document its default in the
    Makefile.am header comment. Optionally run the upstream ctest suite in
    the check stage behind a maintainer flag.
  - Ordering: gmpfrxx_mkII depends on gmp, mpfr, mpc. Wire it into the
    external build order after mpc in external/Makefile.am and into
    configure.ac (AC_CONFIG_FILES etc.) following how the other externals
    are registered. Regenerate autotools via ./gen_configure.sh per repo
    convention.
  - Add --with-system-gmpfrxx-mkII=<prefix> (naming consistent with
    existing --with-system-* options if any; otherwise follow the closest
    existing convention) to skip the bundled build.
  - CMake side of MPLAPACK: teach the CMake build to locate gmpfrxx_mkII
    from external/i/GMPFRXX_MKII (bundled) or a user-provided prefix
    (system), exporting an imported/INTERFACE target for later phases. No
    target links against it yet.

Acceptance:
  - autotools: full configure && make with bundled externals stages headers
    into external/i/GMPFRXX_MKII/include/; make install places them under
    $(prefix)/include/ plus the CMake package files; make dist includes the
    distfile; verify-distfile fails on a corrupted tarball (demonstrate in
    the report, then restore).
  - CMake: configure succeeds in both bundled and system modes; the
    imported target resolves.
  - Everything that was green in P0 baseline still builds.
Deliverables: external/gmpfrxx_mkII/, build-system wiring, REPORT-P1.md.
```

---

## P2 — Complete the interop surface (gmpfrxx_mkII repo)

P2 runs in `nakatamaho/gmpfrxx_mkII`, on top of the branch containing the
current adapters. SPIKE.md from P0 is the requirements document; keep it in
the session context. MPLAPACK-relevant adapter types are exactly:
`dd_real/dd_complex`, `qd_real/qd_complex`, `binary80_real/binary80_complex`,
`binary128_real/binary128_complex`, and the GMP targets
`gmpxx::mpf_class`/`gmpxx::mpfc_class`. Native `double` is already a scalar
leaf and needs tests only. edd/td: leave existing code untouched; do not
extend (extension may follow the same traits mechanically later, outside
this migration).

### P2A — Real arithmetic completion

```text
[COMMON CONTEXT]

Phase P2A. gmpfrxx_mkII repo, branch topic/adapter-arith-real.

Extend the existing external_mpfr_real_traits / external_mpf_real_traits
hook so that for each MPLAPACK-relevant real adapter type T
(dd_real, qd_real, binary80_real, binary128_real):

  1. operator+, operator*, operator/ are provided against mpfr_class and
     against mpf_class, BOTH operand orders, returning expression nodes
     under the existing expression-template contract (same shape as the
     already-implemented operator-).
  2. Compound assignment +=, -=, *=, /= on mpfr_class and mpf_class with a
     T right-hand side.
  3. Mixed expressions compose: (x - dd) * dd + dd style chains must
     evaluate with the same single-rounding discipline as below.

Numerical contract (binding):
  - Ingesting T into MPFR must not perform any intermediate rounding. For
    component types (dd, qd) use mpfr_sum over the components, or, if you
    use another scheme, include a written proof in the header comment that
    no intermediate rounding can occur for arbitrary component exponent
    gaps (double-double components may span the full double exponent
    range). binary80 via mpfr_set_ld, binary128 via mpfr_set_float128, as
    already established.
  - Each arithmetic result is the correctly rounded (current rounding mode)
    value of the exact mathematical result of <materialized operand values>
    — i.e. exactly one rounding at the destination precision per operation
    node, matching the library's existing semantics for its own types.
  - MPF targets route through an MPFR accumulator and a single final
    mpfr_get_f, as the existing adapters do.
  - No path through double/long double except when T is binary80 (long
    double is its storage) — see I4.

Do not change core public headers' responsibilities: everything stays in
opt-in adapter headers and detail/ helpers. Raw long double / __float128
scalar-leaf rejection stays intact.

Acceptance: new unit tests per type per operator per operand order,
including correct-rounding checks against an mpfr_sum-based reference at
mismatched precisions (e.g. destination 24, 53, 113, 256 bits), and at
least one directed test where a dd operand has a large hi/lo exponent gap.
Full ctest green, including all pre-existing 162+ tests and compile-fail
tests. REPORT in docs/ per upstream convention.
```

### P2B — Complex arithmetic completion

```text
[COMMON CONTEXT]

Phase P2B. gmpfrxx_mkII repo, branch topic/adapter-arith-complex.

Mirror P2A for the complex adapters (dd_complex, qd_complex,
binary80_complex, binary128_complex) against mpc_class and mpfc_class:
operator+, *, / in both operand orders, compound assignment, expression
composition. Rounding contract: per-part single rounding consistent with
the library's existing mpc/mpfc semantics; document explicitly how
real*complex and complex/complex ingest the adapter operand (no
intermediate rounding on ingestion, per P2A rules). Same acceptance
structure as P2A, plus mixed real-adapter-with-complex-class cases if the
existing traits already admit them (if they do not, state so in the report
instead of adding new categories).
```

### P2C — MPLAPACK compatibility surface

```text
[COMMON CONTEXT]

Phase P2C. gmpfrxx_mkII repo, branch topic/adapter-mplapack-compat.

Add opt-in headers (NOT in core headers):
  <gmpfrxx_mkII/adapters/mplapack_compat_mpfr.hpp>
  <gmpfrxx_mkII/adapters/mplapack_compat_gmp.hpp>

Contents are dictated by SPIKE.md section (b): the function surface
MPLAPACK's utils layers need on top of mpfr_class/mpc_class and
mpf_class/mpfc_class respectively. From the current inventory that is at
least: sign, castREAL_*, castINTEGER_*, nint, iceil, cabs1, pow2, pow4,
variadic min/max, mod, frac, pi(REAL&-style or value-returning per
SPIKE.md signatures), plus whatever SPIKE.md adds. Match the SPIKE.md
signatures exactly so MPLAPACK's utils headers become thin includes.

GMP side hard requirement (I4): every transcendental/elementary function
exposed here (exp, log, sin, cos, pow, atan2, hypot, log10, log2, ...)
must be computed at working precision via the library's own GMP math
implementations — never via get_d()/double round-trips. Add a directed
test: reconstruct a Givens rotation (c, s = ±sqrt(1-c*c)) at 512-bit
working precision and verify the residual scales with working precision
rather than saturating near 2^-53.

Also verify (and fix if missing) parity of abs/floor/ceil/trunc/fmod/
copysign overload sets between the mpfr and gmp class surfaces, so ADL can
never silently resolve a gmp-side call into std:: via an implicit
conversion.

Acceptance: compat headers compile standalone; new tests green; full ctest
green; a table in the report mapping every SPIKE.md requirement to the
providing header and test.
```

### P2D — Test hardening and release

```text
[COMMON CONTEXT]

Phase P2D. gmpfrxx_mkII repo, branch topic/adapter-release.

1. Round-trip and exactness suite across the MPLAPACK-relevant matrix:
   T -> mpfr_class(prec) -> T identity for prec large enough per type
   (binary80: 64; binary128: 113; dd/qd: verify against an mpfr_sum
   reference rather than assuming a fixed embedding width), plus
   double <-> mpfr/mpf sanity tests.
2. Compile-fail suite: re-verify raw long double and raw __float128
   rejection still holds INCLUDING in the presence of every new operator
   overload from P2A/P2B (overload resolution changes are the risk here).
3. Advisory lint: scripts/check_auto_expr.sh — a grep-level heuristic that
   flags `auto <name> = <expr involving adapter ops>;` patterns in a given
   source tree, for MPLAPACK to run in P3/P4. Advisory, not a build gate.
4. Update SPECIFICATIONS.md/STATUS.md/CHANGES for the new surface; bump
   version; produce the official dist tarball via the existing dist
   machinery; record its SHA256 in the report.

Acceptance: full ctest green on Debug and Release; compile-fail suite
green; tarball builds standalone in a clean directory (configure + install
+ compile a sample TU using the adapters against installed headers).
```

### P1R — Re-vendor the interop release (mplapack repo)

```text
[COMMON CONTEXT]

Phase P1R. Branch topic/gmpfrxx-p1r. Replace the vendored tarball in
external/gmpfrxx_mkII/download/ with the P2D release, update
external/distfiles.sha256, adjust the version variable in
external/gmpfrxx_mkII/Makefile.am. Re-run the P1 acceptance gates
verbatim. No other changes. Deliverable: REPORT-P1R.md with gate outputs.
```

---

## P3 — MPFR backend switch (mplapack repo)

```text
[COMMON CONTEXT]

Phase P3. Branch topic/gmpfrxx-p3. Switch the mpfr backend from mpfrc++ to
gmpfrxx_mkII. FABLE-generated routine bodies are untouched; the switch
happens in the typedef layer, utils, print, benchmarks, examples, and test
drivers, exactly as enumerated in INVENTORY.md.

1. include/mplapack_mpfr.h (and the mpblas counterpart / wherever
   INVENTORY.md locates the REAL/COMPLEX binding): drop mpreal.h,
   mpcomplex.h, gmpxx.h, and `using namespace mpfr;`. Bind REAL to
   mpfrxx::mpfr_class and COMPLEX to mpfrxx::mpc_class (include
   mpfrxx_mkII.h and the MPC complex header per upstream layering), plus
   the adapter headers required by the cast paths INVENTORY.md found
   (dd/qd/binary80/binary128) and the mplapack_compat_mpfr.hpp surface.
   The mpfr backend must no longer include gmpxx.h, directly or
   transitively through MPLAPACK headers.
2. Rewrite include/mplapack_utils_mpfr.h as a thin layer over
   mplapack_compat_mpfr.hpp, preserving every public signature the rest of
   the tree uses (P0 SPIKE.md signatures are normative). Same for the mpfr
   paths in mplapack_print.h and the hex printers.
3. Precision and rounding control: map every site that used
   mpreal::set_default_prec / per-object precision to the mpfrxx TLS
   default-precision API. Document the semantic delta (global vs
   thread-local) in MIGRATION.md draft form, and make the test drivers set
   precision such that numerical results are comparable to the P0
   baseline. If any test spawns threads and relied on the global default,
   set the precision explicitly per thread and note it in the report.
4. Name collision audit: the legacy global ::mpc_class (GMP complex) still
   exists until P5. Every use of the token mpc_class in the mpfr backend
   path must be namespace-qualified or eliminated; add a temporary grep
   gate for unqualified `[^:]mpc_class` in mpfr-backend TUs to the report.
5. Expression-template hygiene: run scripts/check_auto_expr.sh (from P2D)
   over the mpfr backend's non-generated sources, benchmarks, examples,
   and test drivers; fix every true positive by explicit materialization;
   list false positives in the report.
6. Both build systems: link lines, include paths (external/i/GMPFRXX_MKII
   during build; installed prefix for installed-tree examples), removal of
   any mpfrc++ include path from the mpfr backend targets.

Acceptance:
  - full mpfr backend test suite (autotools) passes;
  - tools/compare_baseline.py against the P0 baseline exits 0 under
    baseline_rules.tsv; any rule you must loosen requires a written
    numerical justification in the report and review by the maintainer —
    do not loosen silently;
  - CMake mpfr backend builds and passes the P0 smoke subset;
  - grep gates: no gmpxx.h include and no unqualified mpc_class in the
    mpfr backend path; `mpreal` / `mpcomplex` tokens appear nowhere in the
    mpfr backend path (mpfrc++ itself still exists in-tree for the gmp
    backend's sake only if INVENTORY.md showed such coupling; otherwise it
    is simply unused until P5).
Deliverables: the switch, REPORT-P3.md with all gate outputs.
```

---

## P4 — GMP backend switch (mplapack repo)

```text
[COMMON CONTEXT]

Phase P4. Branch topic/gmpfrxx-p4. Switch the gmp backend off gmpxx.h and
include/mpc_class.h.

1. include/mplapack_gmp.h and the gmp typedef layer: replace
   #include <gmpxx.h> with gmpxx_mkII.h; REAL binds to gmpxx::mpf_class,
   COMPLEX binds to gmpxx::mpfc_class. include/mpc_class.h is no longer
   included by the gmp backend (file deletion happens in P5).
2. Rewrite include/mplapack_utils_gmp.h over mplapack_compat_gmp.hpp,
   preserving public signatures. Replace
   include/mplapack_gmp_transcendents.h usage with the mkII math surface;
   whatever workarounds it contained for the double-fallback era must be
   either obsoleted (preferred) or ported with a comment explaining why
   they are still needed.
3. Directed regression for the historical bug (I4): at 512-bit working
   precision, the Givens/rotation-style constructions used by test matrix
   generators must show error scaling with working precision, not a
   ~2^-53 cap. Add this as a permanent test if the suite lacks one.
4. Baseline: this phase MAY legitimately change numbers where the old code
   silently detoured through double. Every deviation from the P0 baseline
   must be classified in baseline_rules.tsv as allowed-improvement with a
   one-line justification per rule; compare_baseline.py must exit 0 with
   zero unclassified deviations. Accuracy regressions are not acceptable
   in any case.
5. Remove -lgmpxx and libgmpxx references from configure.ac, Makefile.am
   tree, CMake tree, pkg-config templates, packaging/ — everywhere
   INVENTORY.md found them. GMP itself may still be built with
   --enable-cxx if other externals expect it; changing external/gmp
   configure flags is allowed only if nothing else in the tree needs
   libgmpxx, and must be justified in the report.
6. Run scripts/check_auto_expr.sh over the gmp backend path; fix true
   positives.

Acceptance: full gmp backend suite (autotools) green;
compare_baseline.py exit 0; CMake gmp backend builds + smoke subset; grep
gates: no gmpxx.h anywhere in the gmp backend path, no -lgmpxx anywhere in
the tree, no unqualified legacy mpc_class in the gmp backend path.
Deliverables: the switch, REPORT-P4.md.
```

---

## P5 — Excision of mpfrc++ and legacy headers (mplapack repo)

```text
[COMMON CONTEXT]

Phase P5. Branch topic/gmpfrxx-p5. Point of no return; both backends are
already green on gmpfrxx_mkII.

1. Delete the mpfrc++/ tree and include/mpc_class.h. Remove every install
   rule, EXTRA_DIST entry, include-path reference, and packaging reference
   to them (INVENTORY.md is the checklist; re-run inventory.sh to catch
   drift since P0).
2. Licensing: mpfrc++ is LGPL and shipped its own copying files. Update
   AUTHORS/README/COPYING structure and any license aggregation in
   packaging/ to drop mpfrc++ and to state the gmpfrxx_mkII license for
   the vendored external. Do not touch the licenses of other externals.
3. Commit docs/migration/gmpfrxx_mkII/tools/purge_gate.sh: exits nonzero
   if any of the following match outside an explicit whitelist
   (docs/migration/, CHANGES.*, MIGRATION.md, ChangeLog, doc/ history):
     #include <gmpxx.h>     #include "gmpxx.h"
     \bmpreal\b             \bmpcomplex\b
     mpfrc\+\+              include/mpc_class\.h      -lgmpxx
   The whitelist is a committed file, not inline in the script.
4. MIGRATION.md: downstream-facing section with the mapping table
   (mpreal -> mpfrxx::mpfr_class, mpcomplex -> mpfrxx::mpc_class,
   ::mpf_class -> gmpxx::mpf_class, ::mpc_class -> gmpxx::mpfc_class,
   cast2* naming, precision API mapping, the `auto` pitfall with a
   before/after example, header names to include).

Acceptance: purge_gate.sh exit 0; full test suites for mpfr and gmp
backends green under both build systems' gate levels (full autotools,
CMake smoke); make dist / CMake package produce trees free of mpfrc++.
Deliverables: deletions, purge_gate.sh + whitelist, MIGRATION.md section,
REPORT-P5.md.
```

---

## P6 — Fable / codegen alignment (mplapack repo, independent)

```text
[COMMON CONTEXT]

Phase P6. Branch topic/gmpfrxx-p6. Scope: the fable/ pipeline and any
code-generation templates (cout.py and friends) only.

1. Audit fable templates, emitted-header lists, and generator scripts for
   references to mpreal, mpcomplex, gmpxx.h, mpfrc++, mpc_class.h. Fix
   the templates, not the outputs.
2. If any template changed, regenerate the affected outputs through the
   pipeline and verify the regeneration is diff-clean against the
   committed tree except for the intended reference changes. Per AGENTS.md
   this is the ONLY sanctioned way generated sources may change; no hand
   edits.
3. Extend purge_gate.sh's search space to cover fable/ and the generator
   scripts if P5's version excluded them.

Acceptance: regeneration diff report in REPORT-P6.md (empty diff, or
intended-only diff with the template commit that explains it);
purge_gate.sh exit 0; mpfr and gmp suites still green.
```

---

## P7 — CI matrix (mplapack repo)

```text
[COMMON CONTEXT]

Phase P7. Branch topic/gmpfrxx-p7. Update CI (GitHub Actions and any other
committed CI) to cover the migrated tree:

  axes: {Linux, macOS, MinGW cross} x {autotools, CMake}
        x {bundled gmpfrxx_mkII, system gmpfrxx_mkII}
        (system-mode job may install the P2D tarball as its "system" copy)
  plus: binary80 jobs restricted to x86-64 runners;
        one job with the MPFR_WANT_FLOAT128 knob in the non-default
        position;
        purge_gate.sh and compare_baseline.py (mpfr+gmp, against the
        committed baseline) run as mandatory steps on the Linux/autotools
        job.

Prune redundant combinations rather than exploding the matrix; justify
each pruning in the workflow file comments. Do not weaken any existing
job's coverage of the non-migrated backends (dd/qd/double/binary80/
binary128).

Acceptance: all workflows green on the PR; matrix documented in
REPORT-P7.md.
```

---

## P8 — Release engineering: MPLAPACK 3.0.0 (mplapack repo)

```text
[COMMON CONTEXT]

Phase P8. Branch topic/gmpfrxx-p8.

1. CHANGES.3.0.0.md: the migration, the ABI break for mpfr/gmp backends,
   removal of mpfrc++/gmpxx.h, the vendored gmpfrxx_mkII version, the
   double-fallback fix class, edd/td non-scope note.
2. Version bump in BOTH configure.ac and CMakeLists.txt; libtool
   -version-info bump reflecting the ABI break for the affected libraries
   only (per AGENTS.md, -version-info is independent of package version —
   compute current:revision:age correctly and show the computation in the
   report).
3. MIGRATION.md finalized (P5 draft + precision-semantics note from P3).
4. make dist / CMake packaging: tarball self-contained including the
   vendored gmpfrxx_mkII distfile; clean-room test: unpack the dist
   tarball in an empty directory, configure with bundled externals, build
   mpfr backend, run its smoke subset.
5. Release checklist file for the maintainer (tagging and GitHub Release
   are performed by the maintainer, not by the agent).

Acceptance: clean-room dist build passes; version/ABI numbers consistent
across build systems; checklist complete. Deliverables: release artifacts
prep, REPORT-P8.md.
```

---

## Appendix A — Type mapping (normative)

| legacy (removed in 3.0.0)         | replacement                       |
|-----------------------------------|-----------------------------------|
| `mpfr::mpreal`                    | `mpfrxx::mpfr_class`              |
| `mpfr::mpcomplex`                 | `mpfrxx::mpc_class`               |
| `::mpf_class` (libgmpxx)          | `gmpxx::mpf_class`                |
| `::mpc_class` (include/mpc_class.h GMP complex) | `gmpxx::mpfc_class` |
| `cast2dd_real` etc. (mpfrc++ flavor) | `mpfrxx::cast2dd_real` / `gmpxx::cast2*` per adapter headers |
| raw `long double` scalar          | `binary80` wrapper (adapter)      |
| raw `__float128` scalar           | `binary128` wrapper (adapter)     |

## Appendix B — Numerical contract summary

1. Adapter-operand ingestion into MPFR: no intermediate rounding
   (mpfr_sum for component types, or a written no-double-rounding proof).
2. One correct rounding per operation node at destination precision.
3. MPF/MPFC targets: MPFR accumulator + single final transfer.
4. No double/long-double detours except for binary80's native storage.
5. Baseline deviations only via classified rules in baseline_rules.tsv;
   allowed-improvement requires per-rule justification; regressions never.