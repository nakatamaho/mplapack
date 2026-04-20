# MPLAPACK 2.2.1 — Changes since 2.2.0

Release date: 2026-04 (tag: v2.2.1, branch: `release/2.2`)

This is a focused maintenance release for the 2.2 series.

---

## MPFR exponent range handling

- Removed MPLAPACK's automatic MPFR exponent range truncation.
  MPLAPACK no longer changes the MPFR global `emin`/`emax` range to
  `±(precision * 64)` during initialization.
- `MPLAPACK_MPFR_PRECISION` now changes only the MPFR working precision.
  It does not imply any exponent range adjustment.
- Explicit exponent range configuration remains supported through
  `MPLAPACK_MPFR_EMIN` and `MPLAPACK_MPFR_EMAX`.  The MPFR binary64 and
  binary128 emulation test profiles continue to use these variables to
  request IEEE-like exponent ranges.
- `ArithmeticParams<mpfr::mpreal>` continues to derive `rmin`, `rmax`,
  `sfmin`, `safmin`, `safmax`, and Blue-scaling inputs from the current
  MPFR global exponent range at runtime.

This keeps the default MPFR backend faithful to the user's MPFR runtime
configuration.  Tests that need an IEEE-like exponent range should set it
explicitly rather than relying on MPLAPACK initialization to clamp it.
