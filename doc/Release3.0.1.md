## MPLAPACK 3.0.1 Release process

This release uses the clean source archive generated from release baseline
`953d7a4916554546937a753a30b0619691072841`.  The result snapshot and this
document are committed with the final release tag.

| Action | Date | Status | Description |
| --- | --- | --- | --- |
| Freeze source snapshot | 2026-09-10 | done | Clean source at `953d7a4916554546937a753a30b0619691072841`. |
| Verify release archive | 2026-09-10 | done | `mplapack-3.0.1.tar.xz` verified with the SHA256 below. |
| Tier1 build and QA | 2026-09-10 onward | accepted | 12 platform/compiler jobs completed with `OK`. |
| Tier1 result review | 2026-09-15 | done | 126,515,712 `lin` tests passed; known `eig` threshold-edge results documented. |
| Tier3 C++ standard builds | 2026-09-10 onward | done | Ubuntu 26.04 amd64 and arm64 with GCC 15.2.0. |
| Tarball smoke test | 2026-09-10 onward | done | Ubuntu 24.04 amd64 and arm64 builds completed. |
| Finalize release notes | 2026-09-15 | done | See [CHANGES.3.0.1.md](../CHANGES.3.0.1.md). |
| Create Git tag `v3.0.1` | 2026-09-15 | done | Annotated tag on the final result/documentation commit. |
| Create GitHub Release | — | pending | Upload the archive and checksum files separately. |
| Post-release verification | after release | pending | Download assets, verify checksums, and repeat a tarball smoke build. |

### Release artifact

- Source archive: `release/logs/20260910_175329/source/mplapack-3.0.1.tar.xz`
- SHA256 file: `release/logs/20260910_175329/source/mplapack-3.0.1.tar.xz.sha256sum`
- MD5 file: `release/logs/20260910_175329/source/mplapack-3.0.1.tar.xz.md5sum`
- SHA256: `47ebb653b21f0c62e8144c1e515e94d76965d9d0d4b7ba216b034d778570cbaa`
- MD5: `8963437912a012d9f00e972c4d2dfde4`
- Source metadata: `release/logs/20260910_175329/source/source-metadata.txt`

### QA records

The complete release QA records are retained under:

```text
release/logs/20260910_175329/
```

The Tier1 CSV files contain 12 `OK` results, and the two Tier3 CSV files also
contain `OK`.  The two Ubuntu 24.04 tarball smoke logs completed successfully.
The aggregate Tier2 `results.csv` records 15 skipped entries and two failed
legacy MinGW reference-build entries; Tier2 is not release-gating.

The committed numerical results are:

```text
mplapack/test/eig/results/3.0.1/
mplapack/test/lin/results/3.0.1/
```

They contain 116,811,306 recognized `eig` tests with 30 known threshold-edge
exceedances and 126,515,712 recognized `lin` tests with no failures.
