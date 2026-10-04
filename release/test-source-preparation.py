#!/usr/bin/env python3
"""Check release source preparation without running remote builds."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

MAKEFILE = Path(os.environ.get(
    "QA_RELEASE_MAKEFILE", str(Path(__file__).resolve().with_name("Makefile"))
))
SNAPSHOT = """#!/bin/sh
set -eu
echo snapshot >> "$QA_EVENTS"
mkdir -p "$(dirname "$MPLAPACK_SOURCE_INFO_FILE")"
sleep 0.1
printf 'MPLAPACK_SOURCE_MODE=dist\\nMPLAPACK_SOURCE_TARBALL=%s\\n' "$QA_EXPECTED" > "$MPLAPACK_SOURCE_INFO_FILE"
"""
RUNNER = """#!/bin/sh
set -eu
test "${MPLAPACK_SOURCE_MODE:-dist}" = dist
test "${MPLAPACK_SOURCE_TARBALL:-}" = "$QA_EXPECTED"
echo "$0 $*" >> "$QA_EVENTS"
"""
BUILD = """#!/bin/sh
set -eu
test "${TARBALL:-}" = "$QA_EXPECTED"
echo tarball >> "$QA_EVENTS"
exit "${QA_TARBALL_RC:-0}"
"""


def check(targets, *, explicit=False, tarball_failure=False):
    with tempfile.TemporaryDirectory(prefix="mplapack-release-source-") as directory:
        root = Path(directory)
        shutil.copyfile(MAKEFILE, root / "Makefile")
        for name, content in {
            "make-source-snapshot.sh": SNAPSHOT,
            "run-remote-targets-by-host.sh": RUNNER,
            "run-remote-macos.sh": RUNNER,
            "run-remote-linux-docker.sh": RUNNER,
            "build-all.sh": BUILD,
        }.items():
            path = root / name
            path.write_text(content)
            path.chmod(0o755)
        events = root / "events"
        expected = str(root / "input.tar.xz")
        env = os.environ.copy()
        for name in list(env):
            if name.startswith("MPLAPACK_") or name in (
                "MAKEFLAGS", "MFLAGS", "MAKELEVEL", "TARBALL", "SOURCE_INFO_FILE"
            ):
                env.pop(name)
        env.update(
            LOGDIR=str(root / "logs"),
            MPLAPACK_SOURCE_SNAPSHOT_REUSE="no",
            MPLAPACK_DIST_CACHE_REUSE="no",
            QA_EVENTS=str(events),
            QA_EXPECTED=expected,
            QA_TARBALL_RC="7" if tarball_failure else "0",
        )
        command = ["make", "-j", *targets]
        if explicit:
            command.append("TARBALL=" + expected)
        result = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True)
        lines = events.read_text().splitlines() if events.exists() else []
        assert (result.returncode != 0) == tarball_failure, result.stdout + result.stderr
        assert lines.count("snapshot") == (0 if explicit else 1), lines
        if tarball_failure:
            assert not any("run-remote" in line for line in lines), lines
        else:
            for target in targets:
                if target.startswith("tier1"):
                    expected_target = {
                        "tier1": "tier1-macos-amd64",
                        "tier1-macos": "tier1-macos-amd64",
                        "tier1-linux": "tier1-ubuntu2404-amd64",
                    }.get(target, target)
                    assert any(expected_target in line for line in lines), lines
            if "tier2" in targets:
                assert lines.count("tarball") == 1, lines
                assert any("tier2-ubuntu2604-cxxstd" in line for line in lines), lines
        print("PASS:", " ".join(command))


check(["tier1", "tier2"])
check(["tier1-macos", "tier1-linux"])
check([
    "tier1-macos-arm64", "tier1-macos-amd64",
    "tier1-ubuntu2404-arm64", "tier1-ubuntu2604-arm64",
    "tier1-ubuntu2404-amd64", "tier1-ubuntu2604-amd64",
    "tier1-ubuntu2404-mingw64-amd64", "tier1-ubuntu2604-mingw64-amd64",
    "tier1-ubuntu2404-inteloneapi-amd64", "tier1-ubuntu2604-inteloneapi-amd64",
    "tier1-debian12-i386", "tier1-debian13-i386",
    "tier2-ubuntu2604-cxxstd-arm64", "tier2-ubuntu2604-cxxstd-amd64",
])
check(["tier2"], explicit=True)
check(["tier2"], tarball_failure=True)
check(["tier2"], explicit=True, tarball_failure=True)
