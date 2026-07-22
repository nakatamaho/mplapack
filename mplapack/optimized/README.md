# Optimized MPLAPACK sources

Each backend directory contains backend-specific optimized LAPACK drivers for
the corresponding `libmplapack_<backend>_opt` library.

A file named `X.cpp` replaces `mplapack/reference/X.cpp` in that backend's
optimized library. If the optimized implementation keeps the reference
implementation as a fallback, name the fallback helper `X_ref.cpp`; helper
files do not shadow public routines. Every replacement must preserve the
public signature of the reference routine.
