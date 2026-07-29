# MPLAPACK LAPACK Examples

| Directory | LAPACK problem class | Routines demonstrated |
|---|---|---|
| `00_GeneralLinearEquations` | General linear systems | `Rgesv`, `Cgesv`, `Rgetrf/Rgetrs`, `Cgetrf/Cgetrs`, `Rgetri`, `Cgetri`, `Rgbsv`, `Rgtsv`, `Rgesvx`, `Rgecon` |
| `01_PositiveDefiniteLinearEquations` | Positive definite linear systems | `Rposv`, `Cposv`, `Rpotrf/Rpotrs`, `Cpotrf/Cpotrs`, `Rptsv`, `Rpocon` |
| `02_LeastSquares` | Least squares and constrained least squares | `Rgels`, `Cgels`, `Rgelsy`, `Rgelss`, `Cgelss`, `Rgglse`, `Rggglm` |
| `03_SymmetricEigenproblems` | Symmetric/Hermitian eigenproblems | Existing `*_test` examples |
| `04_NonsymmetricEigenproblems` | Nonsymmetric eigenproblems | Existing `*_test` and matrix-file examples |
| `05_SingularValueDecomposition` | Singular value decomposition | Existing `*_test` and matrix-file examples |
| `06_SymmetricIndefiniteLinearEquations` | Symmetric/Hermitian indefinite systems | `Rsysv`, `Csysv`, `Chesv` |
| `07_GeneralizedSymmetricDefiniteEigenproblems` | Generalized symmetric/Hermitian definite eigenproblems | `Rsygv`, `Chegv` |
| `08_GeneralizedNonsymmetricEigenproblems` | Generalized nonsymmetric eigenproblems | `Rggev`, `Rgges`, `Cggev` |
| `09_GeneralizedSingularValueDecomposition` | Generalized singular value decomposition | `Rggsvd3` |
| `90_PrecisionComparison` | Precision comparison studies | `Rlamch`, Hilbert, Wilkinson, Vandermonde, Pascal, Kahan, Frank, `Rgejsv` vs `Rgesvd` |

## Building

In an in-tree autotools build, configure MPLAPACK with examples enabled and run:

```sh
make -C examples
```

For installed MPLAPACK examples, copy or enter an installed category directory under `share/examples/mplapack/<Category>` and use the platform Makefile, for example:

```sh
make -f Makefile.linux
make -f Makefile.macos
```

## Adding An Example

Only files under `generic/` are hand-maintained templates. Write a new `*_generic.cpp` in the target category, using the pseudo-types `REAL`, `COMPLEX`, and `INTEGER`. The generator substitutes those tokens for each backend and rewrites `Rlamch` to the backend-specific symbol. Use `%%MPLIB%%` when the template needs to print the backend name.

Template filenames beginning with `C` receive the complex precision header. Every other template receives the real precision header, so do not give a real template a leading `C` name.

Regenerate a category with:

```sh
cd examples/mplapack/<Category>/generic
bash ../../generic/generate.sh
```

Or regenerate all examples from `examples/` with:

```sh
bash gen_all.sh
```

Commit both the template and the generated outputs. Outside `generic/`, category generated `.cpp` files and per-category Makefiles are generated and should not be hand-edited. Matrix data files named `Matrix_*.txt` are preserved by the generator cleanup.
