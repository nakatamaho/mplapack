#!/usr/bin/env python3
"""Generate the per-backend Makefile.am files from backends.txt.

For every backend in backends.txt this writes

    mpblas/optimized/<name>/Makefile.am
    mpblas/test/<name>/Makefile.am
    mplapack/test/compare/<name>/Makefile.am
    benchmark/Makefile.<name>.am

from the templates in misc/backend_makefiles/.  Flags follow the backend's
build traits (backends.txt); settings that do not follow from the traits are
in misc/backend_makefiles/backends.conf.

    misc/gen_backend_makefiles.py          write the files
    misc/gen_backend_makefiles.py --check  exit 1 if any file is out of date

Template syntax: %%b%% and %%B%% are the backend name in lower and upper
case, %%key%% is a value, and a line holding only "%%block key%%" is replaced
by zero or more lines.
"""

import argparse
import configparser
import pathlib
import re
import sys

TOP = pathlib.Path(__file__).resolve().parent.parent
TEMPLATE_DIR = TOP / "misc" / "backend_makefiles"

OUTPUTS = [
    ("mpblas_optimized", "mpblas/optimized/{b}/Makefile.am"),
    ("mpblas_test", "mpblas/test/{b}/Makefile.am"),
    ("compare", "mplapack/test/compare/{b}/Makefile.am"),
    ("benchmark", "benchmark/Makefile.{b}.am"),
]

CONF_KEYS = {
    "accelerator", "smoke_check_only", "smoke_extra", "compare_tolerance",
    "compare_rlamch_per_abi", "compare_test_env", "bench_quadmath_conditional",
    "bench_note",
}

# The MPFR backend is the test oracle; every other backend's tests link it.
ORACLE = "mpfr"


def read_backends():
    backends = []
    for line in (TOP / "backends.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, real, complex_, default, traits = line.split()
        traits = [] if traits == "-" else traits.split(",")
        backends.append((name, traits))
    return backends


def read_conf(names):
    conf = configparser.ConfigParser(interpolation=None)
    conf.read(TEMPLATE_DIR / "backends.conf")
    for section in conf.sections():
        if section not in names:
            sys.exit(f"backends.conf: [{section}] is not in backends.txt")
        unknown = set(conf[section]) - CONF_KEYS
        if unknown:
            sys.exit(f"backends.conf: [{section}]: unknown keys {sorted(unknown)}")
    return conf


def values(name, traits, conf):
    """Template values for one backend."""
    c = conf[name] if conf.has_section(name) else {}
    qd = "qd" in traits
    nofma = "nofma" in traits
    gmp = "gmp" in traits
    mpfr = "mpfr" in traits
    v = {}

    # mpblas/optimized
    if gmp:
        v["opt_cppflags"] = ("$(GMPFRXX_MKII_CPPFLAGS) -I$(top_srcdir)/include $(OPENMP_CXXFLAGS) "
                             "-I$(GMP_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_%%B%%  -finline-functions")
        v["opt_libadd"] = "$(GMP_LIBS) $(OPENMP_LIBS)"
    elif mpfr:
        v["opt_cppflags"] = ("$(GMPFRXX_MKII_CPPFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) "
                             "-I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) $(OPENMP_CXXFLAGS) "
                             "-DMPLAPACK_BUILD_WITH_%%B%%")
        v["opt_libadd"] = "$(MPLAPACK_MPFR_LIBS) $(OPENMP_LIBS)"
    else:
        v["opt_cppflags"] = ("-I. -I$(top_srcdir)/include -I$(QD_INCLUDEDIR) $(OPENMP_CXXFLAGS) "
                             "-DMPLAPACK_BUILD_WITH_%%B%%")
        v["opt_libadd"] = ("$(QD_LIBADD) " if qd else "") + "$(OPENMP_LIBS)"
    v["opt_cxxflags"] = ["libmplapack_%%b%%_opt_la_CXXFLAGS = $(DD_CXXFLAGS)"] if nofma else []
    v["opt_mingw_ldflags"] = " -lquadmath" if "binary128libs" in traits else ""
    accelerator = c.get("accelerator", "")
    if accelerator:
        subdir, conditional = accelerator.split()
        v["opt_accelerator"] = ["SUBDIRS =", f"if {conditional}", f"SUBDIRS += {subdir}", "endif", ""]
    else:
        v["opt_accelerator"] = []

    # mpblas/test and mplapack/test/compare link the MPFR oracle as well.
    if name == ORACLE:
        v["oracle_archive"] = v["oracle_lib"] = v["oracle_opt_lib"] = ""
    else:
        v["oracle_archive"] = ",$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a"
        v["oracle_lib"] = " -lmplapack_mpfr"
        v["oracle_opt_lib"] = " -L$(top_builddir)/mplapack/reference -lmplapack_mpfr"

    # mpblas/test
    smoke = "arithmetic.test_%%b%% complex.test_%%b%%"
    if c.get("smoke_extra", ""):
        smoke += " " + c["smoke_extra"]
    if c.get("smoke_check_only", "no") == "yes":
        v["test_check_programs"] = [
            f"mpblas_%%b%%_smoke_programs = {smoke}",
            "check_PROGRAMS = $(mpblas_%%b%%_smoke_programs)",
            "mpblas_%%b%%_smoke_TESTS = $(mpblas_%%b%%_smoke_programs)",
        ]
    else:
        v["test_check_programs"] = [
            "check_PROGRAMS = $(mpblas_%%b%%_test_PROGRAMS)",
            f"mpblas_%%b%%_smoke_TESTS = {smoke}",
        ]
    if qd:
        v["test_mplibs"] = ("-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmp "
                            "-L$(QD_LIBDIR) -lqd $(QD_RPATH_LDFLAGS)")
    else:
        v["test_mplibs"] = "-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp"
    v["test_cxxflags"] = (
        "$(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include -I$(GMP_INCLUDEDIR) "
        "-I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR)"
        + (" -I$(QD_INCLUDEDIR)" if qd else "")
        + " -DMPLAPACK_BUILD_WITH_%%B%% -DMPLAPACK_INTERNAL"
        + (" $(DD_CXXFLAGS)" if nofma else ""))

    # mplapack/test/compare
    env = c.get("compare_test_env", "").split()
    if env:
        v["compare_env"] = (
            [f"%%B%%_COMPARE_TEST_ENV = {' '.join(env)}", "backend_tests_environment = \\"]
            + [f"    {a}; export {a.split('=')[0]};" + (" \\" if i < len(env) - 1 else "")
               for i, a in enumerate(env)]
            + [""])
        env_prefix = "$(%%B%%_COMPARE_TEST_ENV) "
    else:
        v["compare_env"] = []
        env_prefix = ""
    tolerance = c.get("compare_tolerance", "")
    per_abi = c.get("compare_rlamch_per_abi", "no") == "yes"
    scripts = []
    for routine in ("Rlaruv", "Rlamch"):
        reference = f"{routine}_reference.txt"
        label = routine
        if routine == "Rlamch" and per_abi:
            reference = "Rlamch_reference@ABI_BITS@.txt"
            label = "Rlamch (@ABI_BITS@-bit)"
        scripts += [f"run_{routine}_test.sh: Makefile", "\t@echo '#!/bin/sh' > $@"]
        if tolerance:
            scripts.append("\t@echo 'NUM_DIFF_PY=\"$(top_srcdir)/misc/num_diff.py\"' >> $@")
        scripts.append(f"\t@echo '{env_prefix}$${{LOG_COMPILER}} ./{routine}.test_%%b%%$${{EXEEXT}} || true' >> $@")
        if tolerance:
            scripts.append(f"\t@echo 'if python3 \"$${{NUM_DIFF_PY}}\" \"$${{srcdir}}/{reference}\" "
                           f"{routine}.txt --tol {tolerance}; then' >> $@")
        else:
            scripts.append(f"\t@echo 'if diff -u $${{srcdir}}/{reference} {routine}.txt > /dev/null 2>&1; then' >> $@")
        scripts += [f"\t@echo '  echo \"PASS: {label}\"; exit 0' >> $@", "\t@echo 'else' >> $@"]
        if tolerance:
            scripts.append(f"\t@echo '  echo \"FAIL: {label}\"; exit 1' >> $@")
        else:
            scripts.append(f"\t@echo '  echo \"FAIL: {label}\"; diff $${{srcdir}}/{reference} {routine}.txt; exit 1' >> $@")
        scripts += ["\t@echo 'fi' >> $@", "\t@chmod +x $@", ""]
    if per_abi:
        scripts.append("EXTRA_DIST = Rlamch_reference32.txt Rlamch_reference64.txt Rlaruv_reference.txt")
    else:
        scripts.append("EXTRA_DIST = Rlaruv_reference.txt Rlamch_reference.txt")
    v["compare_run_scripts"] = scripts
    v["qd_include"] = " -I$(QD_INCLUDEDIR)" if qd else ""
    v["nofma_cxxflags"] = "$(DD_CXXFLAGS)" if nofma else ""
    v["compare_whole_archive_libs"] = "-lquadmath" if "binary128libs" in traits else ""
    v["compare_rpath"] = ["backend_rpath_ldflags = $(QD_RPATH_LDFLAGS)"] if qd else []
    if qd:
        v["compare_mplibs"] = "-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(QD_LIBDIR) -lmpc -lmpfr -lgmp -lqd"
    elif mpfr:
        v["compare_mplibs"] = "-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpc -lmpfr -lgmp"
    else:
        v["compare_mplibs"] = "-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -lmpc -lmpfr -lgmp"

    # benchmark
    if gmp:
        v["bench_cxxflags"] = ("$(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include "
                               "-I$(GMP_INCLUDEDIR) -I$(MPC_INCLUDEDIR) -DMPLAPACK_BUILD_WITH_%%B%%")
        libs = "-L$(GMP_LIBDIR) -lgmp "
    elif mpfr:
        v["bench_cxxflags"] = ("$(GMPFRXX_MKII_CPPFLAGS) $(OPENMP_CXXFLAGS) -I$(top_srcdir)/include "
                               "-I$(GMP_INCLUDEDIR) -I$(MPFR_INCLUDEDIR) -I$(MPC_INCLUDEDIR) "
                               "-DMPLAPACK_BUILD_WITH_%%B%%")
        libs = "-L$(MPC_LIBDIR) -L$(MPFR_LIBDIR) -L$(GMP_LIBDIR) -lmpfr -lmpc -lgmp "
    else:
        v["bench_cxxflags"] = ("$(OPENMP_CXXFLAGS) -I$(top_srcdir)/include"
                               + (" -I$(QD_INCLUDEDIR)" if qd else "")
                               + " -DMPLAPACK_BUILD_WITH_%%B%%"
                               + (" $(DD_CXXFLAGS)" if nofma else ""))
        libs = "-L$(QD_LIBDIR) -lqd " if qd else ""
    quadmath = c.get("bench_quadmath_conditional", "")
    if not quadmath and "binary128libs" in traits:
        quadmath = "BINARY128_USE_QUADMATH"

    def libdepends(extra):
        return [f"%%b%%_libdepends    = $(%%b%%lapack_libdepends) {extra}$(DYLD)",
                f"%%b%%opt_libdepends = -L$(top_builddir)/mpblas/optimized/%%b%% -lmplapack_%%b%%_opt {extra}$(DYLD)"]

    if quadmath:
        v["bench_libdepends"] = (["", f"if {quadmath}"] + libdepends(libs + "-lquadmath ")
                                 + ["else"] + libdepends(libs) + ["endif"])
    else:
        v["bench_libdepends"] = libdepends(libs)
    note = c.get("bench_note", "").strip()
    v["bench_note"] = ([f"# {line}" for line in note.splitlines()] + [""]) if note else []
    return v


def render(template_name, name, v):
    text = (TEMPLATE_DIR / f"{template_name}.am.in").read_text()
    out = []
    for line in text.split("\n"):
        m = re.fullmatch(r"%%block (\w+)%%", line)
        if m:
            out.extend(v[m.group(1)])
        else:
            out.append(line)
    text = "\n".join(out)
    text = text.replace("%%template%%", f"{template_name}.am.in")

    def scalar(m):
        key = m.group(1)
        if key == "b":
            return name
        if key == "B":
            return name.upper()
        value = v[key]
        assert isinstance(value, str), key
        return value

    # Values may themselves contain %%b%% and %%B%%.
    for _ in range(2):
        text = re.sub(r"%%(\w+)%%", scalar, text)
    if "%%" in text:
        sys.exit(f"{template_name}.am.in: unresolved placeholder for {name}")
    return "\n".join(line.rstrip() for line in text.split("\n"))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="report out-of-date files instead of writing them")
    args = parser.parse_args()

    backends = read_backends()
    conf = read_conf({name for name, _ in backends})
    stale = []
    for name, traits in backends:
        v = values(name, traits, conf)
        for template_name, pattern in OUTPUTS:
            path = TOP / pattern.format(b=name)
            text = render(template_name, name, v)
            if path.exists() and path.read_text() == text:
                continue
            if args.check:
                stale.append(path.relative_to(TOP))
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                print(f"wrote {path.relative_to(TOP)}")
    if stale:
        print("out of date (run misc/gen_backend_makefiles.py):", file=sys.stderr)
        for path in stale:
            print(f"  {path}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
