## MPLAPACK 3.0.0 Release process

This release uses a clean source archive generated from release baseline
`eec51cd962d452563f6bc8703a441d0c3d410c88`.  The archive is accompanied by
SHA256 and MD5 checksum files.

| Action | Date | Status | Description |
| --- | --- | --- | --- |
| Freeze source snapshot | 2026-08-29 | done | Clean source at `eec51cd962d452563f6bc8703a441d0c3d410c88`. |
| Verify release archive | 2026-08-29 | pending | Generate `mplapack-3.0.0.tar.xz` and its checksum files from the release commit. |
| Tier1 build and QA | 2026-08-24 onward | accepted | 11 platform/compiler triplets were covered by the preceding QA snapshot; the release-only QD fix was not rerun through the full matrix. |
| Tier1 `lin` result review | 2026-08-29 | done | 126,515,712 recognized tests; no failures. |
| Tier1 `eig` result review | 2026-08-29 | accepted | 30 known threshold-edge exceedances in 28 records; documented as non-harmful and non-blocking. |
| Tier3 C++ standard builds | 2026-08-24 onward | done | Ubuntu 26.04 amd64 and arm64 with GCC 15.2.0. |
| Tarball smoke test | 2026-08-24–25 | done | Ubuntu 24.04 amd64 and arm64 remote Docker builds completed successfully. |
| Finalize release notes | 2026-08-29 | done | See [CHANGES.3.0.0.md](../CHANGES.3.0.0.md). |
| Create Git tag `v3.0.0` | 2026-08-29 | pending | Tag the exact source commit selected for the release. |
| Create GitHub Release | 2026-08-29 | pending | Upload the regenerated tarball, `.sha256sum`, and `.md5sum` files. |
| Post-release verification | after release | pending | Download the published assets and verify checksums and a tarball smoke build. |

### Release artifact locations

- Source archive: `release/logs/20260824_231945/source/mplapack-3.0.0.tar.xz`
- SHA256 file: `release/logs/20260824_231945/source/mplapack-3.0.0.tar.xz.sha256sum`
- MD5 file: `release/logs/20260824_231945/source/mplapack-3.0.0.tar.xz.md5sum`
- Reusable dist cache: `release/dist-cache/git-eec51cd962d452563f6bc8703a441d0c3d410c88/`

The source archive was supplied to Docker jobs from the host release context;
it is not embedded in the Docker images.
