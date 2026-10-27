# GOAL 09 (P2): triangular inverse and inverses from the factors

Prerequisites: goals 02, 03 (CPU part), 05-06 (GPU part). Read
00-overview.md first.

## Scope

CPU: dedicated {R,C}trtri and {R,C}lauum in mplapack/optimized/common
(recursive, parallel; default bit-identical, fast mode opt-in), potri /
getri keep the reference code on top of them. GPU: Rtrtri, Rlauum, Rpotri,
Rgetri `_<b>_cuda` for dd, td, qd.

## Acceptance criteria

- Rfactor_test extended to trtri / lauum / potri / getri (bit identity in
  default mode, ‖A·inv(A) − I‖ check in fast mode).
- mplapack/test/lin Rchktr / Rchkpo / Rchkge inverse parts pass.
- Timings in the commit body.
- One commit per side is allowed: "... (factorization goal 09, CPU)" and
  "... (factorization goal 09, GPU)".
