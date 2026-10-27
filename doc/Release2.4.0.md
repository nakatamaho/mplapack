# MPLAPACK 2.4.0 maintenance release preparation

## Exact release source

- Series: 2.x maintenance; changes since the latest published 2.3.x tag, `v2.3.0`.
- QA commit: `00e0b659d01d073a36a58fa5851fae067fdd3a10`.
- QA source label: `00e0b659d01d-bfdd8c64498a`.
- Source metadata reports a clean tree, with no tracked patch.
- Tested artifact: `release/logs/20261020_143347/source/mplapack-2.4.0.tar.xz`.
- SHA-256: `bfdd8c64498a1b9472b98d3f40bcf94a51fb3e54f2921996cf1eb2c85c6093df`.
- Package versions in Autotools/CMake: 2.4.0. Libtool: 2:0:0. SONAME: 2.

Use the **tested tarball**, not a newly generated `make dist` archive and not
GitHub's automatically generated source archives. These post-QA documentation
files are separate release attachments; they are not in the tested tarball.
Do not silently substitute a later documentation commit as the QA source.

## Status checked during preparation (2026-10-27)

The first QA run (`6de413e4a75f`, `release/logs/20261013_103040/`) found
22 MinGW/GMP EIG failures caused by GMP tests running at 64 bits. The fix
(`00e0b659d`) retains the backend initializers in static test links; that
run's tarball is superseded and must not be published.

- All 12 Tier 1 CSV records: OK.
- All 4 Tier 2 CSV records: OK.
- Tarball SHA-256 checked against QA metadata.
- Published tags/releases checked: latest 2.3.x release is v2.3.0; no
  v2.4.0 tag existed locally or on origin.
- Numerical result aggregation: LIN 0/126,515,712; EIG 8/116,310,045 failed
  (6 records), all matching the exceptions published for 2.3.0.
- GMP runs at 512 bits on all 11 triplets, including both MinGW hosts.
  See [the investigation report](QA.2.4.0-mingw-gmp.md).
- Two MinGW jobs share a result-tree triplet name; both per-host result
  directories were checked separately and agree.
- Review decision (2026-10-27): numerical exceptions accepted as the known
  2.3.0 exceptions; release approved by the maintainer.

## Prepared material

- `CHANGES.2.4.0.md`: maintenance changes and numerical QA exceptions.
- `doc/ANNOUNCE.2.4.0.md`: English announcement draft.
- `release/dist-cache/release-2.4.0/`: local publication bundle containing the
  tested tarball, checksum files, metadata, notes, and numerical summary JSON.

## Verification

From the repository root:

```sh
bundle="$PWD/release/dist-cache/release-2.4.0"
(cd "$bundle" && sha256sum -c mplapack-2.4.0.tar.xz.sha256sum)
cmp release/logs/20261020_143347/source/mplapack-2.4.0.tar.xz \
    "$bundle/mplapack-2.4.0.tar.xz"
python3 misc/summarize_mplapack_tests.py 2.4.0 mplapack/test/eig --only-fail
python3 misc/summarize_mplapack_tests.py 2.4.0 mplapack/test/lin --only-fail
```

## Publication after numerical review and explicit approval

The commands below are a checklist, not steps already executed. Record the
review decision and final publication date here. Update the draft text and
refresh its copies in the local bundle before publishing.

1. Commit the release documentation separately, preserving the QA commit as
   the source target of the release. Preserve the numerical evidence.
2. If `v2.4.0` already exists locally or remotely, inspect it rather than
   replacing or force-pushing it. Otherwise create the tag at the QA commit:

   ```sh
   git tag -a v2.4.0 00e0b659d01d073a36a58fa5851fae067fdd3a10 \
       -m 'MPLAPACK 2.4.0 maintenance release'
   git push origin refs/tags/v2.4.0
   ```

3. Create a **draft** GitHub release. Keep the newer 3.x release as the
   repository's latest release (`--latest=false`).

   ```sh
   bundle="$PWD/release/dist-cache/release-2.4.0"
   gh release create v2.4.0 --verify-tag --draft --latest=false \
       --title 'MPLAPACK 2.4.0 (2.x maintenance release)' \
       --notes-file CHANGES.2.4.0.md \
       "$bundle/mplapack-2.4.0.tar.xz" \
       "$bundle/mplapack-2.4.0.tar.xz.sha256sum" \
       "$bundle/mplapack-2.4.0.tar.xz.md5sum" \
       "$bundle/source-metadata.txt" \
       "$bundle/CHANGES.2.4.0.md"
   ```

4. Check the draft's tag target, attachments, checksum, maintenance scope,
   and numerical caveats. After approval, publish without marking it latest:

   ```sh
   gh release edit v2.4.0 --draft=false --latest=false
   ```

5. Update the README News entry with the actual publication date and release
   URL. Publish the reviewed announcement. Do not claim all numerical tests
   passed or send the announcement while the release remains a draft.
