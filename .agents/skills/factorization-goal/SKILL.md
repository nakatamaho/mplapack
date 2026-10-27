---
name: factorization-goal
description: >
  Execute ONE numbered goal (01-09, 07a) of the MPLAPACK dedicated LU /
  Cholesky work: CPU and GPU implementations of Rgetrf/Cgetrf, Rpotrf/Cpotrf and the
  routines around them (panels, laswp, solvers, inverses). Use ONLY when the
  user explicitly asks to run a factorization goal (e.g. "run factorization
  goal 03", "$factorization-goal 05"). Not for other build, release or
  numerical work in this repository.
---

# factorization-goal — run one goal of the LU / Cholesky series

The user's message specifies the goal number NN (01-09 or 07a). If no number is
given, ask for it and stop.

Specification files live next to this skill:

    .agents/skills/factorization-goal/references/00-overview.md  (normative spec)
    .agents/skills/factorization-goal/references/NN-*.md         (this session's goal)

Procedure — follow strictly, in order:

1. Read AGENTS.md at the repository root (hard rules, the accelerator entry
   point naming rule, the build quickstarts).
2. Read references/00-overview.md in full: routine table, result modes,
   naming, GPU design, verification gates, environment notes.
3. Read references/NN-*.md in full (07a: references/07a-*.md). Its
   "Scope", "Tasks" and "Acceptance criteria" define the session's work.
   Change nothing outside that scope, except link lines you broke.
4. Check the goal's "Prerequisites" (earlier goals' commits, files they
   create). If one is missing, STOP and report; do not improvise.
5. Before editing, verify the goal file's factual claims against the tree
   (existing files, symbols via `nm`, Makefile.am / CMakeLists.txt
   structure). If reality disagrees, reality wins: note the discrepancy in
   the final report and proceed on facts.
6. Implement for BOTH build systems (autotools and CMake) unless the goal
   says otherwise.
7. Run every acceptance criterion, including misc/check_unique_symbols.sh on
   every archive you changed. Without a GPU or nvcc, mark the GPU checks
   UNVERIFIED (CUDA device tests exit 77 without a device) — never
   fabricate results or timings.
8. Commit once, with the message the goal mandates; append to the body the
   acceptance results (pass / fail / UNVERIFIED), the measured timings the
   goal asks for, and any discrepancies found in step 5.
9. Final answer: a terse report — what changed, verification matrix with
   numbers, UNVERIFIED items, whether the next goal is unblocked.

Constraints:
- Never edit FABLE-generated routine bodies in mpblas/reference or
  mplapack/reference; reuse them by including them under another name
  (`#define Rgetrf Rgetrf_reference_cpu` / `#include
  "../../reference/Rgetrf.cpp"`, as mplapack/optimized/common/Rgetrf.cpp).
- The default results of every public routine stay bit-identical to the
  reference unless 00-overview.md explicitly allows otherwise for that path.
- One goal per session.
