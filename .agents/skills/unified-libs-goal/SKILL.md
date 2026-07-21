---
name: unified-libs-goal
description: >
  Execute ONE numbered goal (01-05) of the MPLAPACK unified-library
  migration on branch topic/unified-libs. Use ONLY when the user
  explicitly asks to run a unified-libs goal (e.g. "run goal 02",
  "$unified-libs-goal 03"). Not for general build, release, or
  numerical work in this repository.
---

# unified-libs-goal — run one migration goal

The user's message specifies the goal number NN (01|02|03|04|05).
If no number is given, ask for it and stop.

All specification files live next to this skill:

    .agents/skills/unified-libs-goal/references/00-overview.md   (normative spec)
    .agents/skills/unified-libs-goal/references/NN-*.md          (this session's goal)

Procedure — follow strictly, in order:

1. Read AGENTS.md at the repository root (hard rules).
2. Read references/00-overview.md in full. It is the normative
   specification: target library taxonomy, the basename-shadowing rule
   (accelerator > optimized > reference), the unique-symbol verification
   gate, and the invariants that bind every goal.
3. Read references/NN-*.md in full. Its "Scope", "Tasks", and
   "Acceptance criteria" sections define this session's work. Nothing
   outside that scope may be changed, except fixing link lines you broke.
4. Confirm you are on branch `topic/unified-libs` and that the previous
   goal's commit is present (goal 01 needs none). If either check fails,
   STOP and report; do not improvise a base.
5. Before editing, verify the goal file's factual claims against the tree
   (file counts, symbol definitions via `nm`, current Makefile.am /
   cmake/ structure). If reality disagrees with the goal file, reality
   wins: note the discrepancy in the final report and proceed on facts.
6. Implement the tasks for BOTH build systems (autotools and CMake)
   unless the goal file says otherwise.
7. Run every item in the goal's "Acceptance criteria". Do not skip the
   unique-symbol check (misc/check_unique_symbols.sh). If a toolchain is
   unavailable (nvcc, OpenCL), follow the goal file's environment note
   and mark the affected checks UNVERIFIED — never fabricate results.
8. Commit exactly once, with the commit message mandated by the goal
   file, and append a body listing: acceptance criteria results
   (pass/fail/UNVERIFIED), and any spec discrepancies found in step 5.
9. Final answer: a terse report — what changed, verification matrix,
   UNVERIFIED items, and whether the next goal is unblocked.

Constraints:
- Never edit FABLE-generated routine bodies in mpblas/reference or
  mplapack/reference.
- No backward-compatibility shims; the ABI break is intentional.
- One goal per session. Do not start goal NN+1 even if trivial.
