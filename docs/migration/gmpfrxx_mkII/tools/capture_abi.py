#!/usr/bin/env python3
"""Capture the P0 ABI, package metadata, and downstream consumer results."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path


BASE_SHA = "b875e74d4b927282c907c3a29e6cadda48a7d57b"
BACKENDS = ("gmp", "mpfr", "qd", "dd", "double", "binary80", "binary128")
TARGET_PROPERTIES = (
    "INTERFACE_COMPILE_DEFINITIONS",
    "INTERFACE_INCLUDE_DIRECTORIES",
    "INTERFACE_LINK_LIBRARIES",
    "IMPORTED_LOCATION_RELEASE",
    "IMPORTED_SONAME_RELEASE",
)


class CaptureError(RuntimeError):
    pass


def run(
    command: list[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    input_text: str | None = None,
) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        rendered = shlex.join(command)
        raise CaptureError(
            f"command failed ({completed.returncode}): {rendered}\n{completed.stdout}"
        )
    return completed.stdout


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise CaptureError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def normalize(text: str, prefixes: tuple[Path, ...]) -> str:
    normalized = text
    for prefix in prefixes:
        normalized = normalized.replace(str(prefix), "${prefix}")
    return normalized


def exported_symbols(library: Path) -> tuple[list[str], list[str]]:
    output = run(["nm", "-D", "--defined-only", "--format=posix", str(library)])
    mangled = sorted(
        {
            line.split()[0]
            for line in output.splitlines()
            if line.strip() and len(line.split()) >= 2
        }
    )
    if not mangled:
        raise CaptureError(f"no exported symbols found in {library}")
    demangled_output = run(["c++filt"], input_text="\n".join(mangled) + "\n")
    demangled = sorted(set(demangled_output.splitlines()))
    return mangled, demangled


def soname(library: Path) -> str:
    output = run(["readelf", "-d", str(library)])
    match = re.search(r"\(SONAME\).*?\[([^\]]+)\]", output)
    if not match:
        raise CaptureError(f"missing SONAME in {library}")
    return match.group(1)


def backend_headers(prefix: Path, backend: str) -> dict[str, str]:
    include = prefix / "include" / "mplapack"
    common = {"blas.h", "lapack.h", "blaswrapper.h", "mplapack_config.h"}
    selected = set(common)
    selected.update(path.name for path in include.glob(f"*_{backend}.h"))
    if backend == "mpfr":
        selected.update({"mpreal.h", "mpcomplex.h"})
    elif backend == "gmp":
        selected.update({"mpc_class.h", "mplapack_gmp_transcendents.h"})
    elif backend in {"qd", "dd"}:
        selected.update({"dd_complex.h", "qd_complex.h"})
    missing = sorted(name for name in selected if not (include / name).is_file())
    if missing:
        raise CaptureError(f"{backend}: missing installed public headers: {missing}")
    return {name: sha256(include / name) for name in sorted(selected)}


def pkg_config(prefix: Path, package: str) -> dict[str, object]:
    pc_path = prefix / "lib" / "pkgconfig" / f"{package}.pc"
    if not pc_path.is_file():
        raise CaptureError(f"missing pkg-config file: {pc_path}")
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(prefix / "lib" / "pkgconfig")
    fields: dict[str, object] = {"file_sha256": sha256(pc_path)}
    commands = {
        "version": ["pkg-config", "--modversion", package],
        "cflags": ["pkg-config", "--cflags", package],
        "libs": ["pkg-config", "--libs", package],
        "requires": ["pkg-config", "--print-requires", package],
        "requires_private": ["pkg-config", "--print-requires-private", package],
    }
    for key, command in commands.items():
        fields[key] = normalize(run(command, env=env).strip(), (prefix,))
    return fields


def capture_library(
    prefix: Path, output_dir: Path, library_name: str
) -> dict[str, object]:
    symlink = prefix / "lib" / f"{library_name}.so"
    if not symlink.exists():
        raise CaptureError(f"missing library: {symlink}")
    resolved = symlink.resolve()
    mangled, demangled = exported_symbols(resolved)
    symbol_dir = output_dir / "symbols"
    symbol_dir.mkdir(parents=True, exist_ok=True)
    mangled_path = symbol_dir / f"{library_name}.mangled.txt"
    demangled_path = symbol_dir / f"{library_name}.demangled.txt"
    mangled_path.write_text("\n".join(mangled) + "\n", encoding="utf-8")
    demangled_path.write_text("\n".join(demangled) + "\n", encoding="utf-8")
    return {
        "file": f"${{prefix}}/lib/{symlink.name}",
        "resolved_name": resolved.name,
        "sha256": sha256(resolved),
        "soname": soname(resolved),
        "mangled_symbols": {
            "file": str(mangled_path.relative_to(output_dir)),
            "count": len(mangled),
            "sha256": sha256(mangled_path),
        },
        "demangled_symbols": {
            "file": str(demangled_path.relative_to(output_dir)),
            "count": len(demangled),
            "sha256": sha256(demangled_path),
        },
    }


def cmake_metadata(prefix: Path) -> tuple[dict[str, dict[str, object]], list[str]]:
    matches = sorted(prefix.glob("**/cmake/mplapack/mplapackTargets.cmake"))
    if len(matches) != 1:
        raise CaptureError(
            f"expected one installed mplapackTargets.cmake under {prefix}, got {matches}"
        )
    targets_file = matches[0]
    files = [targets_file]
    files.extend(sorted(targets_file.parent.glob("mplapackTargets-*.cmake")))
    combined = "\n".join(path.read_text(encoding="utf-8") for path in files)
    combined = normalize(combined, (prefix,))

    target_types = {
        target: target_type
        for target, target_type in re.findall(
            r"add_library\((mplapack::[^\s]+)\s+([A-Z]+)\s+IMPORTED\)", combined
        )
    }
    if not target_types:
        raise CaptureError(f"no imported MPLAPACK targets in {targets_file}")

    metadata: dict[str, dict[str, object]] = {
        target: {"type": target_type, "properties": {}}
        for target, target_type in sorted(target_types.items())
    }
    blocks = re.findall(
        r"set_target_properties\((mplapack::[^\s]+)\s+PROPERTIES(.*?)\n\)",
        combined,
        flags=re.DOTALL,
    )
    for target, block in blocks:
        if target not in metadata:
            continue
        properties = metadata[target]["properties"]
        assert isinstance(properties, dict)
        for name in TARGET_PROPERTIES:
            match = re.search(rf"\b{name}\s+\"([^\"]*)\"", block)
            if match:
                properties[name] = match.group(1)

    config_matches = sorted(prefix.glob("**/cmake/mplapack/mplapackConfig.cmake"))
    if len(config_matches) != 1:
        raise CaptureError(
            f"expected one installed mplapackConfig.cmake under {prefix}"
        )
    config_text = config_matches[0].read_text(encoding="utf-8")
    dependencies = sorted(
        set(re.findall(r"find_dependency\(\s*([A-Za-z0-9_]+)", config_text))
    )
    return metadata, dependencies


def target_for_backend(targets: dict[str, dict[str, object]], backend: str) -> dict:
    selected = {
        name: data
        for name, data in targets.items()
        if name in {
            f"mplapack::mplapack_{backend}",
            f"mplapack::mplapack_{backend}_opt",
        }
        or (backend == "binary128" and name.endswith("binary128_opt_opencl"))
    }
    required = f"mplapack::mplapack_{backend}"
    if required not in selected:
        raise CaptureError(f"missing CMake target {required}")
    return selected


def consumer_source(backend: str) -> str:
    return (
        f'#include "mpblas_{backend}.h"\n'
        "int main() {\n"
        f'    return Mlsame_{backend}("N", "n") ? 0 : 1;\n'
        "}\n"
    )


def consumer_env(autotools_prefix: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PKG_CONFIG_PATH"] = str(autotools_prefix / "lib" / "pkgconfig")
    prior = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = (
        str(autotools_prefix / "lib") + (f":{prior}" if prior else "")
    )
    return env


def run_pkg_consumer(prefix: Path, backend: str, root: Path) -> dict[str, str]:
    package = f"mplapack_{backend}"
    dependency_packages = {
        "gmp": ["gmpxx", "gmp"],
        "mpfr": ["mpc", "mpfr", "gmp"],
        "qd": ["qd"],
        "dd": ["qd"],
        "double": [],
        "binary80": [],
        "binary128": [],
    }[backend]
    extra_libraries = {
        "binary80": ["-lquadmath"],
        "binary128": ["-lquadmath"],
    }.get(backend, [])
    directory = root / f"pkg-{backend}"
    directory.mkdir(parents=True)
    source = directory / "main.cpp"
    executable = directory / "consumer"
    source.write_text(consumer_source(backend), encoding="utf-8")
    env = consumer_env(prefix)
    flags = shlex.split(
        run(
            [
                "pkg-config",
                "--cflags",
                "--libs",
                package,
                *dependency_packages,
            ],
            env=env,
        ).strip()
    )
    command = [
        "c++",
        "-std=gnu++17",
        str(source),
        f"-I{prefix / 'include'}",
        "-o",
        str(executable),
        *flags,
        *extra_libraries,
    ]
    run(command, env=env)
    run([str(executable)], env=env)
    return {
        "status": "PASS",
        "package": package,
        "dependency_packages": " ".join(dependency_packages),
        "extra_libraries": " ".join(extra_libraries),
        "compile": normalize(shlex.join(command), (root, prefix)),
        "run": "${work}/consumer",
    }


def run_cmake_consumer(
    autotools_prefix: Path, cmake_prefix: Path, backend: str, root: Path
) -> dict[str, str]:
    directory = root / f"cmake-{backend}"
    source_dir = directory / "source"
    build_dir = directory / "build"
    source_dir.mkdir(parents=True)
    (source_dir / "main.cpp").write_text(consumer_source(backend), encoding="utf-8")
    target = f"mplapack::mplapack_{backend}"
    (source_dir / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.20)\n"
        f"project(p0_consumer_{backend} LANGUAGES CXX)\n"
        "find_package(mplapack CONFIG REQUIRED)\n"
        "add_executable(consumer main.cpp)\n"
        f"target_link_libraries(consumer PRIVATE {target})\n",
        encoding="utf-8",
    )
    env = consumer_env(autotools_prefix)
    prefix_path = f"{cmake_prefix};{autotools_prefix}"
    configure = [
        "cmake",
        "-S",
        str(source_dir),
        "-B",
        str(build_dir),
        f"-DCMAKE_PREFIX_PATH={prefix_path}",
    ]
    run(configure, env=env)
    build = ["cmake", "--build", str(build_dir), "-j32"]
    run(build, env=env)
    executable = build_dir / "consumer"
    run([str(executable)], env=env)
    return {
        "status": "PASS",
        "target": target,
        "configure": normalize(
            shlex.join(configure), (root, autotools_prefix, cmake_prefix)
        ),
        "build": normalize(shlex.join(build), (root,)),
        "run": "${work}/consumer",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--autotools-prefix", required=True, type=Path)
    parser.add_argument("--cmake-prefix", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    autotools_prefix = args.autotools_prefix.resolve()
    cmake_prefix = args.cmake_prefix.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    targets, dependencies = cmake_metadata(cmake_prefix)

    backends: dict[str, dict[str, object]] = {}
    with tempfile.TemporaryDirectory(prefix="mplapack-p0-consumers-") as temp:
        consumer_root = Path(temp)
        for backend in BACKENDS:
            library_names = [
                f"libmplapack_{backend}",
                f"libmplapack_{backend}_opt",
            ]
            if backend == "binary128" and (
                autotools_prefix
                / "lib"
                / "libmplapack_binary128_opt_opencl.so"
            ).exists():
                library_names.append("libmplapack_binary128_opt_opencl")
            libraries = {
                name.removeprefix("lib"): capture_library(
                    autotools_prefix, output_dir, name
                )
                for name in library_names
            }
            backends[backend] = {
                "libraries": libraries,
                "public_headers": backend_headers(autotools_prefix, backend),
                "pkg_config": {
                    f"mplapack_{backend}": pkg_config(
                        autotools_prefix, f"mplapack_{backend}"
                    ),
                    f"mplapack_{backend}_opt": pkg_config(
                        autotools_prefix, f"mplapack_{backend}_opt"
                    ),
                },
                "cmake_targets": target_for_backend(targets, backend),
                "consumers": {
                    "pkg_config": run_pkg_consumer(
                        autotools_prefix, backend, consumer_root
                    ),
                    "cmake": run_cmake_consumer(
                        autotools_prefix, cmake_prefix, backend, consumer_root
                    ),
                },
            }

    manifest = {
        "schema": 1,
        "mplapack_base_sha": BASE_SHA,
        "platform": {
            "system": run(["uname", "-s"]).strip(),
            "machine": run(["uname", "-m"]).strip(),
            "shared_library_format": "ELF",
        },
        "cmake_package_dependencies": dependencies,
        "backends": backends,
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"captured ABI/install manifest for {len(backends)} backends")
    print(f"manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
