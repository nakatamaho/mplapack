# MPLAPACK 3.0.0 gmpfrxx_mkII Migration: P1R Report

## 1. Repository state

- Repository: `/home/docker/mplapack`
- Remote: `git@github.com:nakatamaho/mplapack.git`
- Branch: `topic/gmpfrxx_mkII_migration`
- Approved P1 starting SHA: `2ebf798fd0081ecdc5c1b53fc117431c406bf884`
- Starting worktree: clean

The local and remote migration branch both resolved to the approved P1 SHA before P1R. The accepted P0 and P1 commits were not amended, rebased, squashed, or rewritten.

## 2. Original migration base

The locked original MPLAPACK base remains `b875e74d4b927282c907c3a29e6cadda48a7d57b`. The accepted P0 `LOCK.json` was not changed and continues to describe the bootstrap input used by P0 and P1.

## 3. Published upstream release

- Version and tag: `1.1.0`, `v1.1.0`
- Release title and ID: `gmpfrxx_mkII 1.1.0`, `359956005`
- Draft / prerelease: `false` / `false`
- Published evidence SHA: `24849d87a222c096cbbf3cdaa9a295c828be0aa3`

GitHub API authentication and release metadata validation passed. Exactly one uploaded asset was present and had the required name and size. Its release notes contain the authoritative SHA-256.

## 4. Upstream tag verification

- Annotated tag object: `880748940902be6b2f5fbfce60efa6fd8cc9fa40`
- Peeled commit: `429fd1b35e1927ebaccc9fda5aa2801300b45bf5`

The tag object and peeled target were resolved independently with `git ls-remote --tags` and matched the approved publication state.

## 5. Published asset metadata

- Asset ID: `490164327`
- Filename: `gmpfrxx_mkII.1.1.0.tar.xz`
- Size: `15169540` bytes
- SHA-256: `e0f3b813463b7a45dd493a818c60a17530075e0e647ea02227b75501c1984c73`
- Forbidden archive paths: none

## 6. Downloaded asset

The initial independently downloaded asset was `/tmp/mplapack-gmpfrxx-p1r-download.YZXTsk/gmpfrxx_mkII.1.1.0.tar.xz`. Its filename, size, and SHA-256 matched the published metadata.

## 7. Byte identity

`cmp --silent` confirmed that the downloaded published asset and the committed vendored archive are byte-for-byte identical. The gate repeated a fresh authenticated download and the same comparison successfully.

## 8. Vendored archive replacement

- Removed: `external/gmpfrxx_mkII/download/gmpfrxx_mkII.1.0.1.tar.xz`
- Added: `external/gmpfrxx_mkII/download/gmpfrxx_mkII.1.1.0.tar.xz`

The old bootstrap archive is absent from active vendoring. The new archive is the manually uploaded GitHub Release asset; it was not regenerated, recompressed, or patched.

## 9. Active version-bearing files

- `external/gmpfrxx_mkII/Makefile.am`: package version is `1.1.0`.
- `external/distfiles.sha256`: filename and authoritative digest are `1.1.0`.
- `Makefile.am`: distribution input names the `1.1.0` archive.

The upstream package is configured with `BUILD_TESTING=OFF` in addition to the accepted P1 options. This is required because upstream 1.1.0 enables its own adapter tests through CTest by default, including QD tests, while MPLAPACK P1R intentionally stages only its binding GMP, MPFR, and MPC components. P1R validates the installed package with independent consumers and does not build upstream's optional adapter test matrix. No component was added and no MPLAPACK acceptance condition was disabled.

## 10. Historical references retained

Bootstrap references remain unchanged in the accepted P0/P1 prompt, lock, preconditions, report, cache evidence, P1 gate, and P0 validator. The exact historical paths are classified in `P1R_HISTORICAL_REFERENCES.tsv`. The P1R gate uses that committed path whitelist and rejects old release references in all other active paths.

## 11. P1R release record

`P1R_RELEASE.json` records version `1.1.0`, tag `v1.1.0`, peeled target `429fd1b35e1927ebaccc9fda5aa2801300b45bf5`, evidence SHA `24849d87a222c096cbbf3cdaa9a295c828be0aa3`, archive size `15169540`, the authoritative SHA-256, published non-draft/non-prerelease state, and successful asset verification.

## 12. Internal installation

PASS. Autotools built the committed package in its internal tree and installed it under `external/i/GMPFRXX_MKII` in disposable storage. The installation contained `gmpfrxx_mkII.h`, package configuration, version metadata, and the `gmpfrxx_mkII::gmpfrxx_mkII` exported target. Installed version metadata was exactly `1.1.0`.

## 13. Final-prefix installation

PASS. The second Autotools package tree installed to a separate final prefix. Its public headers, package configuration, exported targets, and version metadata were present. Neither installation contained source/build absolute paths or an unexpected `libgmpxx` dependency.

## 14. Staged-prefix CMake discovery

PASS. MPLAPACK CMake discovered the staged internal package with `MPLAPACK_ENABLE_GMPFRXX_MKII=ON`. CMake only used package discovery; it did not download or build gmpfrxx_mkII.

## 15. Independent-prefix CMake discovery

PASS. A clean installation made directly from the published archive was found through `CMAKE_PREFIX_PATH`. Configure, build, and install completed with dependency auto-fetch OFF and components `GMP,MPFR,MPC`.

## 16. Installed and relocated consumers

PASS. A consumer using `find_package(gmpfrxx_mkII 1.1.0 EXACT CONFIG REQUIRED)` compiled, linked, and ran against the final installation. The prefix was moved, and the same consumer compiled, linked, and ran against the relocated prefix. All component targets and `gmpfrxx_mkII::gmpfrxx_mkII` were present.

## 17. Corruption and restoration

PASS. The gate copied the archive to disposable storage, appended one byte, and verified that SHA-256 validation rejected it. It restored that copy from a fresh published download and verified digest and byte identity. The primary archive was never damaged or rewritten.

## 18. Offline and auto-fetch evidence

PASS. After committed dependency archives were available, the build ran with HTTP, HTTPS, and ALL proxy variables directed to unreachable `127.0.0.1:9`. Internal, final, and independent CMake caches recorded `GMPFRXX_MKII_DEPS_AUTO_FETCH=OFF`. No network-fetched dependency appeared.

## 19. Distribution result

PASS. `make dist` included the exact byte-identical `gmpfrxx_mkII.1.1.0.tar.xz` archive and excluded the old 1.0.1 archive. The vendored archive extracted from the MPLAPACK distribution matched the published asset's SHA-256 and bytes.

## 20. Accepted P0 baseline

PASS. All nine comparator unit tests passed. The immutable `baseline.json` was compared with itself under the accepted rules and reported `baseline comparison: PASS` for all 874 recorded results. No baseline file, rule, numerical output, or threshold changed.

## 21. Active stale-reference scan

PASS. Every old archive/version/digest occurrence was either absent from active configuration or present at an exact path in the committed historical whitelist. No unexpected active bootstrap reference remained.

## 22. P1R gate

Exact command:

```sh
P1R_JOBS=32 bash docs/migration/gmpfrxx_mkII/tools/gate-P1R.sh
```

Result: `gate-P1R: PASS`.

The gate validated release provenance, archive identity, corruption rejection, offline builds, both Autotools installations, both MPLAPACK CMake discovery paths, installed and relocated consumers, distribution inclusion, the accepted P0 baseline, active-reference classification, and primary worktree idempotence.

## 23. Scope confirmation

No backend type or typedef changed. No numerical source, generated routine, Fable generator, public backend header, interop matrix, compile spike, P0/P1 report, P0/P1 gate, baseline rule, or precision contract changed. No MPLAPACK production source/header includes gmpfrxx_mkII. No upstream patch was added.

## 24. Deviations and blockers

Before the complete pass, validation exposed two build-plumbing issues: upstream CTest attempted optional QD adapter tests outside the P1 component set, and out-of-tree `make dist` needed benchmark/example Makefiles configured even though no benchmark/example was built. The final configuration uses `BUILD_TESTING=OFF` and leaves those distribution directories configured. A Python 3.14 absolute-path unittest invocation was changed to the accepted P0 relative-path form. These are gate/build-plumbing corrections only. Blockers: none.

## 25. Worktree status

The complete gate left HEAD, the committed tree, the index, the worktree, and `git diff --check` unchanged and clean.

## 26. Push status

No P1R commit was pushed. The remote migration branch remains at the approved P1 SHA pending human review.
