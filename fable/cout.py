import os
import re
from itertools import product
from io import StringIO
import os.path
import math
import tempfile
import typing
from decimal import Decimal, InvalidOperation

def _env_int(name: str, default: int) -> int:
    v = os.environ.get(name)
    if v is None:
        return default
    try:
        return int(str(v).strip())
    except Exception:
        return default

def _env_str(name: str, default: str) -> str:
    v = os.environ.get(name)
    if v is None:
        return default
    return str(v).strip()

def _env_flag(name: str, default: bool = False) -> bool:
    """Parse a boolean environment variable.

    Accepted truthy values:  1, true, yes, on
    Accepted falsy values:   0, false, no, off, (empty)
    Other values fall back to Python truthiness.
    """
    v = os.environ.get(name)
    if v is None:
        return default
    v = str(v).strip().lower()
    if v in ("1", "true", "yes", "on"):
        return True
    if v in ("0", "false", "no", "off", ""):
        return False
    return bool(v)

# Enable the legacy small CHARACTER*n -> char[] optimization only when requested.
# Default is OFF to prefer fem::str<N> (lower maintenance, fewer edge cases).
FABLE_SMALL_CHAR_ENABLED = _env_flag("FABLE_SMALL_CHAR", default=False)

# Compatibility mode for test conversions: when FABLE_SMALL_CHAR is set to "0",
# emit string view types (str_cref / str_ref) for scalar CHARACTER dummy args,
# and keep CHARACTER*1 scalars as fem::str<1> (not plain char).
#
# This avoids mixing plain 'char' with fem::str_ref in generated test code
# (e.g., c1 = aline(i,i)), and still allows MPLAPACK calls to receive raw
# buffers via ' .elems' when needed.
FABLE_SMALL_CHAR_VIEW = (str(os.environ.get("FABLE_SMALL_CHAR", "")).strip() == "0")

# If set, do not emit COMMON/SAVE boilerplate structs into generated C++.
FABLE_SUPPRESS_COMMON = _env_flag("FABLE_SUPPRESS_COMMON", default=False)

def _parse_ident_list_env(name: str) -> typing.Set[str]:
    """Parse an environment variable as a list of identifiers.

    Splits on commas and whitespace, returns lowercase identifiers.
    """
    v = os.environ.get(name)
    if v is None:
        return set()
    s = str(v).strip()
    if not s:
        return set()
    parts = re.split(r"[,\s]+", s)
    return {p.lower() for p in parts if p}

# COMMON scalars that should be treated as externally provided globals
# instead of being accessed as members of `cmn`.
FABLE_EXTERN_COMMON_SCALARS = _parse_ident_list_env("FABLE_EXTERN_COMMON_SCALARS")
if FABLE_SUPPRESS_COMMON:
    # LAPACK test harness commonly externalizes these.
    FABLE_EXTERN_COMMON_SCALARS.update({"infot", "srnamt", "ok", "lerr", "nout"})

def _load_mplapack_signatures():
    """Load mplapack_signatures.py in a robust way.

    The generated signatures module may live in different locations depending
    on the pipeline (e.g. project root vs fable/ directory). Importing it as a
    top-level module is fragile because sys.path may not include the directory
    that contains the generated file.

    Supported discovery order:
      1) Standard import (sys.path)
      2) $MPLAPACK_SIGNATURES_PY (explicit file path)
      3) ./mplapack_signatures.py (cwd)
      4) <this_dir>/mplapack_signatures.py
      5) <parent_dir>/mplapack_signatures.py

    Returns:
      (FUNCTION_SIGNATURES, FUNCTION_RETURNS) as dicts keyed by lowercase name.
    """
    # 1) Standard import.
    try:
        from mplapack_signatures import FUNCTION_SIGNATURES as _FS, FUNCTION_RETURNS as _FR
        fs = dict(_FS) if _FS else {}
        fr = dict(_FR) if _FR else {}
        return {k.lower(): v for k, v in fs.items()}, {k.lower(): v for k, v in fr.items()}
    except Exception:
        pass

    import importlib.util
    import sys

    candidates = []
    env_path = os.environ.get("MPLAPACK_SIGNATURES_PY")
    if env_path:
        candidates.append(env_path)
    candidates.append(os.path.join(os.getcwd(), "mplapack_signatures.py"))
    this_dir = os.path.dirname(__file__)
    candidates.append(os.path.join(this_dir, "mplapack_signatures.py"))
    candidates.append(os.path.join(
        os.path.dirname(this_dir), "mplapack_signatures.py"))

    for p in candidates:
        if not p:
            continue
        p = os.path.abspath(p)
        if not os.path.isfile(p):
            continue
        try:
            spec = importlib.util.spec_from_file_location(
                "mplapack_signatures", p)
            if spec is None or spec.loader is None:
                continue
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            # Make subsequent imports consistent within this Python process.
            sys.modules["mplapack_signatures"] = mod

            fs = getattr(mod, "FUNCTION_SIGNATURES", {}) or {}
            fr = getattr(mod, "FUNCTION_RETURNS", {}) or {}
            fs = {str(k).lower(): v for (k, v) in dict(fs).items()}
            fr = {str(k).lower(): v for (k, v) in dict(fr).items()}

            if os.environ.get("FABLE_DEBUG_SIGNATURES"):
                import sys as _sys
                print(
                    f"[FABLE] Loaded mplapack_signatures from: {p}", file=_sys.stderr)
                print(
                    f"[FABLE] FUNCTION_SIGNATURES={len(fs)} FUNCTION_RETURNS={len(fr)}", file=_sys.stderr)

            return fs, fr
        except Exception:
            continue

    # Not found: disable signature-based adjustments.
    if os.environ.get("FABLE_DEBUG_SIGNATURES"):
        import sys as _sys
        print(
            "[FABLE] mplapack_signatures not found; signatures disabled.", file=_sys.stderr)
    return {}, {}


FUNCTION_SIGNATURES, FUNCTION_RETURNS = _load_mplapack_signatures()

# Cache inferred signatures for EXTERNAL callable dummy arguments.
# Key: id(fproc) -> dict[callable_name_lower] = (ret_type, [arg_type0, ...])
_INFERRED_CALLABLE_SIGNATURES = {}


def _split_actuals(arg_string: str):
    """Split a C++ argument list string on commas, ignoring commas inside parentheses."""
    parts = []
    buf = []
    depth = 0
    for ch in arg_string:
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            buf.append(ch)
        elif ch == "," and depth == 0:
            part = "".join(buf).strip()
            if part:
                parts.append(part)
            buf = []
        else:
            buf.append(ch)
    if buf:
        part = "".join(buf).strip()
        if part:
            parts.append(part)
    return parts


def _is_array_variable(conv_info, name: str) -> bool:
    """Return True if 'name' is declared as an array in the current procedure.

    Uses conv_info.fproc.fdecl_by_identifier to look up the declaration.
    Returns False if conv_info is None or the variable is not found/not an array.
    """
    if conv_info is None:
        return False
    try:
        fdecl = conv_info.fproc.fdecl_by_identifier.get(name.lower())
        if fdecl is None:
            fdecl = conv_info.fproc.fdecl_by_identifier.get(name)
        if fdecl is None:
            return False
        return getattr(fdecl, "dim_tokens", None) is not None
    except (AttributeError, KeyError):
        return False


def _adjust_actuals_using_signature(arg_string: str, signature, conv_info=None, force_elems_call: bool = False) -> str:
    """Adjust actual arguments based on pointer/value signature.

    For PTR_NUMERIC arguments, if the expression looks like an array
    element (e.g. rwork[...]) and is not already passed by address,
    insert '&' in front of it. Optionally normalize '&name[0]' into
    'name' (pointer to the first element) for pointer parameters only.

    For REF_SCALAR arguments (scalar reference like REAL& alpha),
    the expression should be name[0] if it's an array, or name if scalar.
    - If input is '&name[...]', remove '&' -> 'name[...]'
    - If input is 'name[...]', keep as is
    - If input is bare 'name' AND it's an array, convert to 'name[0]'

    For PTR_CHAR_IN / PTR_CHAR_OUT arguments, if the expression is a bare identifier
    (e.g. normin) and not already passed by address, insert '&' in
    front of it. String literals are left unchanged.

    If conv_info is given and the identifier is a CHARACTER dummy
    argument (const char* in the generated interface), it is left
    unchanged because it is already a pointer.
    """
    parts = _split_actuals(arg_string)

    # Be conservative: if the lengths do not match, do nothing.
    if len(parts) != len(signature):
        return arg_string

    new_parts = []

    def _elems_suffix_for_identifier(name: str) -> str:
        """Return correct accessor to obtain a (const) char* buffer.
        - For view-style CHARACTER dummies (fem::str_cref/fem::str_ref), elems is a
          member function: use '.elems()'.
        - For fixed-length fem::str<N>, elems is a data member array: use '.elems'.
        """
        # In view-mode test conversions (FABLE_SMALL_CHAR=0), core MPLAPACK routines
        # expect raw (const) char* buffers. Force '.elems' on fem::str/str_view actuals
        # for those calls to avoid '&fem::str' type mismatches.
        if (conv_info is not None
                and _is_dummy_character_arg(conv_info, name)
                and not _is_plain_character_pointer_dummy(conv_info, name)):
            return ".elems()"
        return ".elems"
    for part, kind in zip(parts, signature):
        s = part.lstrip()

        if kind == "REF_SCALAR":
            # REF_SCALAR: scalar reference (e.g., REAL& alpha).

            if s.startswith("&"):
                # &name[...] -> name[...] (remove the '&')
                m = re.fullmatch(
                    r"&\s*([A-Za-z_][A-Za-z0-9_]*\s*\[[^\]]+\])", s)
                if m:
                    leading = part[:len(part) - len(s)]
                    inner = re.sub(r"\s+", "", m.group(1))
                    part = leading + inner
                    new_parts.append(part)
                    continue

            if "[" in s:
                new_parts.append(part)
                continue

            m = re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", s)
            if m and _is_array_variable(conv_info, s):
                leading = part[:len(part) - len(s)]
                part = leading + s + "[0]"

            new_parts.append(part)
            continue

        if kind == "PTR_NUMERIC":
            # Already passing by address.
            if s.startswith("&"):
                # Normalize '&name[0]' -> 'name' only for pointer params.
                m = re.fullmatch(
                    r"&\s*([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*0\s*\]", s)
                if m:
                    leading = part[:len(part) - len(s)]
                    part = leading + m.group(1)
                new_parts.append(part)
                continue

            # Only handle array element expressions; plain pointer
            # variables (e.g. x) are left unchanged.
            if "[" in s:
                leading = part[:len(part) - len(s)]
                part = leading + "&" + s

        elif kind in ("PTR_CHAR", "PTR_CHAR_IN", "PTR_CHAR_OUT"):
            # fem::str<N> is not implicitly convertible to (const) char*.
            # For character-pointer parameters, pass the underlying buffer.
            m = re.fullmatch(r"&\s*([A-Za-z_][A-Za-z0-9_]*)$", s)
            if m and (
                _is_fem_str_scalar(conv_info, m.group(1))
                or (conv_info is not None and _is_dummy_character_arg(conv_info, m.group(1))
                    and not _is_plain_character_pointer_dummy(conv_info, m.group(1)))
                or (conv_info is not None and _is_scalar_character_fem_str(conv_info, m.group(1)))
            ):
                leading = part[:len(part) - len(s)]
                name = m.group(1)
                part = leading + name + _elems_suffix_for_identifier(name)
                new_parts.append(part)
                continue

            if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*$", s) and (
                _is_fem_str_scalar(conv_info, s)
                or (conv_info is not None and _is_dummy_character_arg(conv_info, s)
                    and not _is_plain_character_pointer_dummy(conv_info, s))
                or (conv_info is not None and _is_scalar_character_fem_str(conv_info, s))
            ):
                leading = part[:len(part) - len(s)]
                part = leading + s + _elems_suffix_for_identifier(s)
                new_parts.append(part)
                continue
            # Already passing by address.
            if s.startswith("&"):
                new_parts.append(part)
                continue

            # Fortran substring on plain CHARACTER dummy (emitted as const char*):
            #   name(i,j)  ->  name + (i-1)
            # This matches LAPACK-style usage (e.g. LSAMEN/LSAME) where the length
            # is explicit or only leading characters are inspected.
            m = re.fullmatch(
                r"([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*([^,()]+?)\s*,\s*([^,()]+?)\s*\)",
                s,
            )
            if m and conv_info is not None and (
                (_is_fem_str_scalar(conv_info, m.group(1)) or _is_scalar_character_fem_str(conv_info, m.group(1)))
                or (_is_dummy_character_arg(conv_info, m.group(1)) and not _is_plain_character_pointer_dummy(conv_info, m.group(1)))
            ):
                base = m.group(1)
                a = m.group(2).strip()
                b = m.group(3).strip()
                leading = part[:len(part) - len(s)]
                # Substring yields a view type (str_ref/str_cref): elems is a function.
                part = leading + f"{base}({a}, {b}).elems()"
                new_parts.append(part)
                continue

            if m and conv_info is not None and _is_plain_character_pointer_dummy(conv_info, m.group(1)):
                base = m.group(1)
                a = m.group(2).strip()  # start index (Fortran 1-based)
                leading = part[:len(part) - len(s)]
                try:
                    off = int(a) - 1
                    if off == 0:
                        part = leading + base
                    else:
                        part = leading + f"{base} + {off}"
                except Exception:
                    part = leading + f"{base} + (({a}) - 1)"
                new_parts.append(part)
                continue

            # String literal: do not add '&'.
            if s.startswith('"') or s.startswith("'"):
                new_parts.append(part)
                continue

            # Array element expression: name[...] -> &name[...]
            # IMPORTANT: only do this when the base identifier is a CHARACTER dummy
            # argument emitted as a plain (const) char* (scalar or CHARACTER*1 array),
            # or a small CHARACTER*n scalar mapped to char[].
            # NOTE: use fullmatch so we do NOT prefix '&' to compound expressions
            # like 'job[0] + compz[0]' (this would create invalid pointer arithmetic).
            m = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*)\s*\[[^\]]+\]\s*$", s)
            if m:
                base = m.group(1)
                if (base in small_char_identifiers
                        or (conv_info is not None and _is_plain_character_pointer_dummy(conv_info, base))):
                    leading = part[:len(part) - len(s)]
                    part = leading + "&" + s
                    new_parts.append(part)
                    continue

            # Bare identifier (e.g. normin) -> &normin,
            # unless it is a CHARACTER dummy argument (already a pointer),
            # or a small CHARACTER*n scalar we mapped to char[].
            if re.match(r"[A-Za-z_][A-Za-z0-9_]*$", s):
                # CHARACTER dummy arguments may be plain (const) char* or view types.
                if conv_info is not None and _is_dummy_character_arg(conv_info, s):
                    if _is_plain_character_pointer_dummy(conv_info, s):
                        new_parts.append(part)
                        continue
                    leading = part[:len(part) - len(s)]
                    part = leading + s + _elems_suffix_for_identifier(s)
                    new_parts.append(part)
                    continue
                # Scalar CHARACTER*n locals emitted as fem::str<N> do NOT
                # implicitly convert to (const) char*. Pass the fixed buffer.
                if conv_info is not None and _is_scalar_character_fem_str(conv_info, s):
                    leading = part[:len(part) - len(s)]
                    part = leading + s + _elems_suffix_for_identifier(s)
                    new_parts.append(part)
                    continue
                # Small CHARACTER*n scalars mapped to char[] decay to char*,
                # so do NOT add '&' here (jbcmpz -> const char*).
                if s in small_char_identifiers:
                    new_parts.append(part)
                    continue
                leading = part[:len(part) - len(s)]
                part = leading + "&" + s

        new_parts.append(part)

    return ", ".join(p.strip() for p in new_parts)


def _resolve_name_map_path(filename: str) -> str:
    """Resolve the path to a name-map file.

    Search order:
      1) current working directory
      2) directory containing this module

    This supports pipelines that `cd` into a tools directory before running
    `python -m fable.command_line.cout`.
    """
    candidates = [
        os.path.join(os.getcwd(), filename),
        os.path.join(os.path.dirname(__file__), filename),
    ]
    for p in candidates:
        try:
            if os.path.isfile(p):
                return p
        except Exception:
            pass
    return candidates[-1]


def _load_mplapack_name_map_file(path: str) -> dict:
    """Load Fortran -> MPLAPACK C++ routine name mapping from a text file.

    Format (whitespace separated, # for comments):

        dgemm           Rgemm
        ztrsv           Ctrsv
        lsame           Mlsame
        xerbla          Mxerbla
        dcabs1          RCabs1
    """
    mapping = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                src, dst = parts[0], parts[1]
                mapping[str(src).lower()] = str(dst)
    except OSError:
        # Map files are optional; fall back to default naming rules.
        pass
    return mapping


def _load_mplapack_name_maps():
    """Load core/testing name maps.

    Policy:
      - View mode (FABLE_SMALL_CHAR=0): read BOTH
          * mplapack_name_map.txt
          * mplapack_testing_name_map.txt
      - Non-view mode (!=0): read ONLY
          * mplapack_name_map.txt

    The maps are kept separate so we can apply special call-site rules
    (e.g., forcing '.elems' for core MPLAPACK routines in view mode).
    """
    core_path = _resolve_name_map_path("mplapack_name_map.txt")
    core = _load_mplapack_name_map_file(core_path)

    testing = {}
    if FABLE_SMALL_CHAR_VIEW:
        testing_path = _resolve_name_map_path("mplapack_testing_name_map.txt")
        testing = _load_mplapack_name_map_file(testing_path)

    combined = dict(core)
    # Keep core precedence if a key accidentally overlaps.
    for k, v in testing.items():
        if k not in combined:
            combined[k] = v
    return core, testing, combined


_MPLAPACK_NAME_MAP_CORE, _MPLAPACK_NAME_MAP_TESTING, _MPLAPACK_NAME_MAP = _load_mplapack_name_maps()

# Reverse lookup: C++ routine name (lowercase) -> Fortran name (lowercase)
_MPLAPACK_CPP_TO_FORTRAN = {
    str(cpp).lower(): str(f_name).lower() for (f_name, cpp) in _MPLAPACK_NAME_MAP.items()
}

# C++-side routine names that belong to the CORE map (lowercase).
# Used to force '.elems' for fem::str/str_view actuals when FABLE_SMALL_CHAR=0.
_MPLAPACK_CORE_CPP_NAMES = {str(cpp).lower() for cpp in _MPLAPACK_NAME_MAP_CORE.values()}


# Track COMPLEX-typed C++ identifiers in the current procedure.
# Names are C++ identifiers after vmapping (e.g. "alpha", "a", "cmn.z").
complex_identifiers = set()
complex_pointer_identifiers = set()
# small fixed-length CHARACTER scalars mapped to char[]
small_char_identifiers = set()
small_char_identifier_lengths = {}  # name -> int length for small char[]

def _fable_small_char_max_len() -> int:
    """Return max CHARACTER*n length to map to a plain C char[].

    Controlled by env var FABLE_SMALL_CHAR.

    Accepted values:
      - unset                 : default 10 (historical behavior)
      - integer string        : use that value (<= 1 disables small-char arrays)
      - 'false'/'off'/'no'/0  : disable
      - 'true'/'on'/'yes'     : default 10
    """
    raw = os.environ.get("FABLE_SMALL_CHAR")
    if raw is None:
        return 10
    s = str(raw).strip().lower()
    if s in ("", "true", "on", "yes"):
        return 10
    if s in ("false", "off", "no", "0"):
        return 0
    try:
        return int(s)
    except ValueError:
        return 10

_FABLE_SMALL_CHAR_MAX_LEN = _fable_small_char_max_len()

# -----------------------------------------------------------------------------
# Machine-constant-style intrinsics in F90 PARAMETER expressions
#
# Some LAPACK F90 sources define scaling constants via RADIX/MINEXPONENT/
# MAXEXPONENT. Our preprocessing sometimes has to drop multi-line PARAMETER
# expressions to keep the legacy parser alive, which previously resulted in
# silent bogus initializers (e.g. 0.0).
#
# Policy: if we detect these intrinsics in a PARAMETER initializer (either
# directly in the converted expression, or via a dropped continuation), emit
# "UNHANDLED" in the generated C++ so the build fails loudly and the human can
# replace it with a correct MPLAPACK-level implementation.
# -----------------------------------------------------------------------------

_MACHINE_CONST_INTRINSICS = ("radix", "minexponent", "maxexponent")

# Lowercased parameter identifiers that must be forced to UNHANDLED.
_FORCE_UNHANDLED_PARAMETER_NAMES = set()


def _contains_machine_const_intrinsics(text: str) -> bool:
    """Return True if 'text' mentions RADIX/MINEXPONENT/MAXEXPONENT."""
    if not text:
        return False
    low = text.lower()
    return any(k in low for k in _MACHINE_CONST_INTRINSICS)


def _mplapack_default_name(name: str) -> str:
    """Fallback rule: s/d -> R, c/z -> C, others unchanged."""
    if not name:
        return name
    lower = name.lower()
    head, tail = lower[0], lower[1:]
    if head in ("s", "d"):
        return "R" + tail
    if head in ("c", "z"):
        return "C" + tail
    return name


def convert_function_name_to_mplapack(name):
    """Convert top-level routine name using mplapack_name_map.txt and default rule."""
    lower = name.lower() if name else name
    mapped = _MPLAPACK_NAME_MAP.get(lower) if name else None
    if mapped is not None:
        return mapped
    return _mplapack_default_name(name)


def _lookup_routine_signature(name: str):
    """Lookup FUNCTION_SIGNATURES with name normalization.

    The signatures file is generated from MPLAPACK headers, so keys are
    typically the C++-side routine names in lowercase (e.g. 'rgemm', 'ctgsy2').
    We also accept Fortran-ish names and apply convert_function_name_to_mplapack().
    """
    if not FUNCTION_SIGNATURES:
        return None
    if name is None:
        return None

    base = str(name).split("::")[-1].strip()
    key1 = base.lower()
    sig = FUNCTION_SIGNATURES.get(key1)
    if sig is not None:
        return sig

    # Try MPLAPACK-mapped name (e.g. dgemm -> Rgemm, xerbla -> Mxerbla).
    mapped = convert_function_name_to_mplapack(base)
    if mapped:
        return FUNCTION_SIGNATURES.get(str(mapped).lower())
    return None


fmt_comma_placeholder = chr(255)


class _AutoType:
    def __repr__(self):
        return "Auto"


Auto = _AutoType()


class group_args:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        args = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"group_args({args})"


class mutable:
    """
    Minimal replacement for libtbx.mutable.

    Usage:
        flag = mutable(value=False)
        flag.value = True
    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        args = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"mutable({args})"


fmt_comma_placeholder = chr(255)


def break_line_if_necessary(callback, line, max_len=80, min_len=70):
    """Break a long C++ line into shorter ones, but never split slice macros.

    Lines containing __SLICE__ or __SLICE2D__ are left untouched so that
    later postprocessing (_postprocess_slice_assignment) can still see
    the full statement on a single physical line.
    """
    def cb_finalize(s):
        callback(s.replace(fmt_comma_placeholder, ","))

    # Do not break lines that contain slice markers.
    # These must stay on a single physical line for pattern matching.
    if "__SLICE2D__" in line or "__SLICE__" in line:
        cb_finalize(line)
        return

    nc = len(line)
    if nc <= max_len:
        cb_finalize(line)
        return

    # Find first non-space character to determine indent.
    for i_start in range(nc):
        if line[i_start] != " ":
            break
    else:
        raise AssertionError

    lsw = line.startswith
    # Do not break full-line comments.
    if lsw("//", i_start):
        cb_finalize(line)
        return

    potential_break_points = []
    ic = i_start
    while ic < nc:
        c = line[ic]
        if "\"'".find(c) >= 0:
            # Skip over string/char literal.
            q = c
            ic += 1
            while ic < nc:
                prev_c = c
                c = line[ic]
                ic += 1
                if c == q and prev_c != "\\":
                    break
            else:
                raise AssertionError
        elif c == "(":
            ic += 1
            potential_break_points.append((0, ic))
        elif lsw(", ", ic):
            ic += 2
            potential_break_points.append((1, ic))
        elif lsw(" = ", ic) or lsw("), ", ic):
            ic += 3
            potential_break_points.append((1, ic))
        elif (lsw(" + ", ic) or lsw(" - ", ic) or
              lsw(" * ", ic) or lsw(" / ", ic)):
            ic += 3
            potential_break_points.append((0, ic))
        elif lsw(" && ", ic) or lsw(" || ", ic):
            ic += 4
            potential_break_points.append((0, ic))
        elif lsw("//", ic):
            # Do not split before a comment.
            break
        else:
            ic += 1
    potential_break_points.append((0, nc))

    n = nc - i_start
    denom = max_len - i_start - 2
    if n <= 0 or denom <= 0:
        l = max_len
    else:
        import math
        blocks = math.ceil(float(n) / float(denom))
        l_est = float(n) / float(blocks)
        l = max(min_len, int(round(l_est)))

    b = 0
    f = 0  # indent for continuation lines

    def break_more_if_necessary(s):
        # This is mostly for very long string literals; keep as-is.
        while f + len(s) > max_len and s.startswith('"'):
            i = max_len - 2 - f
            j = s.rfind(fmt_comma_placeholder, 0, i)
            if j > 4:
                i = j + 1
            else:
                for j in range(i - 1, -1, -1):
                    if s[j] != "\\":
                        if (i - j) % 2 == 0:
                            i -= 1
                        break
                else:
                    raise AssertionError
            cb_finalize(" " * f + s[:i] + '"')
            s = '"' + s[i:]
        cb_finalize(" " * f + s)

    if lsw("if (", i_start):
        indent_width = 4
    else:
        indent_width = 2

    pprio = 0
    pp = 0
    for ip in range(len(potential_break_points)):
        prio, p = potential_break_points[ip]

        def following_point_is_better():
            for jp in range(ip, len(potential_break_points)):
                prio2, p2 = potential_break_points[jp]
                if prio2 == 1 and p2 - b + f <= max_len:
                    return True
            return False

        if (p - b + f > l and b != pp and
                (pprio == 1 or not following_point_is_better())):
            s = line[b:pp].rstrip()
            if f == 0:
                cb_finalize(s)
                f = i_start + indent_width
            else:
                break_more_if_necessary(s=s)
            b = pp
        pprio = prio
        pp = p
    if b < nc:
        break_more_if_necessary(s=line[b:])


def break_lines(cpp_text, prev_line=None):
    prev_line = [prev_line]
    result = []

    def callback(line):
        if (prev_line[0] is None
            or line != prev_line[0]
                or not line.lstrip().startswith("//C")):
            result.append(line)
            prev_line[0] = line
    for line in "\n".join(cpp_text).splitlines():
        break_line_if_necessary(callback=callback, line=line)
    return result


class dynamic_parameter_props(object):

    __slots__ = ["name", "ctype", "default"]

    def __init__(O, name, ctype, default):
        O.name = name
        O.ctype = ctype
        O.default = default


def create_buffer_blocks(
        target_number_of_blocks,
        buffers,
        min_lines_per_block=100):
    if (target_number_of_blocks >= len(buffers)):
        return [[buffer] for buffer in buffers]
    numbers_of_lines = [len(lines) for lines in buffers]
    sum_lines = sum(numbers_of_lines)
    lines_per_block = max(
        min_lines_per_block,
        sum_lines / target_number_of_blocks)
    result = []
    block_sum_lines = 0
    j = 0
    for i, n in enumerate(numbers_of_lines):
        if (block_sum_lines + n > lines_per_block):
            if (i == j or block_sum_lines <= min_lines_per_block):
                result.append(buffers[j:i+1])
                block_sum_lines = 0
                j = i+1
            else:
                result.append(buffers[j:i])
                block_sum_lines = n
                j = i
        else:
            block_sum_lines += n
    if (j < len(buffers)):
        result.append(buffers[j:])
    assert sum([len(block) for block in result]) == len(buffers)
    return result


def show_traceback():
    import traceback
    print(traceback.format_exc(limit=None))


def strip_leading_zeros(string):
    for i in range(len(string)):
        if (string[i] != "0"):
            return string[i:]
    if (len(string) == 0):
        return ""
    return "0"


def escape_string_literal(s):
    return (s
            .replace("\\", "\\\\")
            .replace('"', '\\"')
            .replace("\t", "\\t")
            .replace("??", "\\?\\?"))


def convert_complex_literal(vmap, tok):
    assert len(tok.value) == 4
    cc = []
    for part in tok.value[1:3]:
        c = []
        sign_tok, val_tok = part
        if (sign_tok is not None):
            c.append(sign_tok.value)
        c.append(convert_token(vmap=vmap, leading=None, tok=val_tok))
        cc.append("".join(c))
    # Map Fortran complex literal (a, b) to C++ COMPLEX(a, b)
    return "COMPLEX(%s)" % ", ".join(cc)


def _normalize_leading_dot_float_literal(s: str) -> str:
    """Normalize leading-dot decimals: .5 -> 0.5, -.5 -> -0.5.

    This is a cosmetic normalization to keep generated C++ consistent.
    It intentionally avoids touching member access like `obj.5`.
    """
    return re.sub(r'(?<![0-9A-Za-z_])([+-]?)\.(\d)', r'\g<1>0.\g<2>', s)


def _format_decimal_float_literal(tv: str, *, is_double_precision: bool) -> str:
    """Pretty-print Fortran floating literals for C++.

    - Unifies REAL and DOUBLE PRECISION formatting.
    - Converts Fortran D-exponents to e-exponents.
    - Optionally expands scientific notation (1e-3 -> 0.001) when reasonable.

    Controls (optional env vars):
      FABLE_FLOAT_LITERAL_STYLE = auto|fixed|scientific|original
      FABLE_FLOAT_FIXED_EXP_MAX = max |exp10| to expand in auto mode (default 20)
      FABLE_FLOAT_FIXED_MAXLEN  = max output length for fixed in auto mode (default 120)
    """
    s = tv.strip()

    # Normalize Fortran exponent letters to 'e'.
    # For DOUBLE PRECISION, D (or d) is the canonical exponent marker.
    # For REAL, E (or e) is canonical, but we defensively normalize D as well.
    if is_double_precision:
        s = s.replace('D', 'd').replace('d', 'e')
    else:
        # Only rewrite D/d when it looks like an exponent marker.
        s = re.sub(r'([dD])([+-]?\d+)\s*$', lambda m: 'e' + m.group(2), s)

    style = os.environ.get('FABLE_FLOAT_LITERAL_STYLE', 'auto').strip().lower()
    try:
        fixed_exp_max = int(os.environ.get('FABLE_FLOAT_FIXED_EXP_MAX', '20'))
    except Exception:
        fixed_exp_max = 20
    try:
        fixed_max_len = int(os.environ.get('FABLE_FLOAT_FIXED_MAXLEN', '120'))
    except Exception:
        fixed_max_len = 120

    # Extract base-10 exponent from the literal text (if present).
    m = re.search(r'[eE]([+-]?\d+)\s*$', s)
    exp10 = 0
    if m:
        try:
            exp10 = int(m.group(1))
        except Exception:
            exp10 = 0

    if style == 'original':
        out = s
    else:
        try:
            d = Decimal(s)
        except (InvalidOperation, ValueError):
            out = s
        else:
            if d.is_zero():
                return '0.0'
            if d == 1:
                return '1.0'

            fixed = format(d, 'f')
            if '.' in fixed:
                fixed = fixed.rstrip('0').rstrip('.')

            sci = format(d.normalize(), 'E').replace('E', 'e').replace('e+', 'e')

            if style == 'fixed':
                out = fixed
            elif style == 'scientific':
                out = sci
            else:
                # auto: expand only if exponent is moderate and string does not explode.
                if abs(exp10) <= fixed_exp_max and len(fixed) <= fixed_max_len:
                    out = fixed
                else:
                    out = sci

    out = _normalize_leading_dot_float_literal(out)

    # Ensure it looks like a floating literal in C++.
    if ('e' not in out and 'E' not in out
            and '.' not in out
            and 'nan' not in out.lower()
            and 'inf' not in out.lower()):
        out += '.0'
    return out

def convert_token(vmap, leading, tok, had_str_concat=None):
    tv = tok.value
    if tok.is_identifier():
        # Apply vmap first (Fortran name -> C++ name or other mapped names)
        raw = vmap.get(tv, tv)
        # Case-insensitive lookup in MPLAPACK name map (routine and helper names).
        lname = raw.lower()
        # Special-case LAPACK F90 helper (from LA_XISNAN module).
        # We translate LA_ISNAN(x) into a C++ helper Mla_isnan(x).
        # The caller must provide Mla_isnan for the active REAL type.
        if lname == "la_isnan":
            return "Mla_isnan"
        mapped = _MPLAPACK_NAME_MAP.get(lname)
        if mapped is not None:
            return mapped
        # No special mapping: keep the vmap result as-is.
        return raw

    if (tok.is_op()):
        if (tv == ".not."):
            return "!"
        if (tv == ".and."):
            return " && "
        if (tv == ".or."):
            return " || "
        if (tv == ".eqv."):
            return " == "
        if (tv == ".neqv."):
            return " != "
        if (tv in ["+", "-"]):
            if (leading):
                return tv
            return " "+tv+" "
        if (tv == "*"):
            if (leading):
                return "star "
            return " "+tv+" "
        if (tv == "/"):
            return " "+tv+" "
        if (tv == "//"):
            if (had_str_concat is None):
                tok.raise_internal_error()
            had_str_concat.value = True
            return " + "
        if (tv == ":"):
            if (leading):
                return "1, "
            return ", "
        if (tv == ".eq." or tv == "=="):
            return " == "
        if (tv == ".ne." or tv == "/="):
            return " != "
        if (tv == ".lt." or tv == "<"):
            return " < "
        if (tv == ".le." or tv == "<="):
            return " <= "
        if (tv == ".gt." or tv == ">"):
            return " > "
        if (tv == ".ge." or tv == ">="):
            return " >= "
        tok.raise_not_supported()
    if (tok.is_string()):
        s = '"' + escape_string_literal(tok.value) + '"'
        if (had_str_concat is None or not had_str_concat.value):
            return s
        return "fem::str_cref(%s)" % s
    if (tok.is_logical()):
        if (tv == ".false."):
            return "false"
        return "true"
    if (tok.is_integer()):
        return tv
    if (tok.is_hexadecimal()):
        return "0x"+tv
    if (tok.is_real()):
        return _format_decimal_float_literal(tv, is_double_precision=False)
    if (tok.is_double_precision()):
        return _format_decimal_float_literal(tv, is_double_precision=True)

    if (tok.is_complex()):
        return convert_complex_literal(vmap=vmap, tok=tok)
    tok.raise_not_supported()


class major_types_cache(object):

    __slots__ = ["identifiers"]

    def __init__(self):
        self.identifiers = set()

    def __contains__(self, value):
        # fem::major_types names are ignored in this build
        return False


major_types = major_types_cache()

cpp_keywords = set("""\
and and_eq asm auto bitand bitor bool break case catch char class compl const
const_cast continue default delete do double dynamic_cast else enum explicit
export extern false float for friend goto if inline int long mutable namespace
new not not_eq operator or or_eq private protected public register
reinterpret_cast return short signed sizeof static static_cast struct switch
template this throw true try typedef typeid typename union unsigned using
virtual void volatile wchar_t while xor xor_eq argv argc
""".split())


def prepend_identifier_if_necessary(identifier):
    if (identifier in major_types or identifier in cpp_keywords):
        return "identifier_" + identifier
    return identifier


def convert_to_mplapack_type(ctype):
    """Convert C++ types to MPLAPACK types (INTEGER, REAL, etc.)"""
    # Remove const and & from the type
    ctype_clean = ctype.replace("const", "").replace("&", "").strip()

    # Convert basic types
    if ctype_clean == "int":
        return "INTEGER"
    elif ctype_clean == "double":
        return "REAL"
    elif ctype_clean == "float":
        return "REAL"
#    elif ctype_clean == "bool":
#        return "LOGICAL"
    elif ctype_clean == "std::complex<float>":
        return "COMPLEX"
    elif ctype_clean == "std::complex<double>":
        return "COMPLEX"

    # For complex types like arr_ref<double>, extract the inner type
    if ctype_clean.startswith("arr_ref<"):
        # Extract type from arr_ref<...>
        inner_type = ctype_clean[8:-1].strip()
        if ", " in inner_type:
            inner_type = inner_type.split(",")[0].strip()
        return convert_to_mplapack_type(inner_type)

    return ctype_clean


def produce_comment_given_sl(callback, sl):
    if (sl.stmt_offs is None):
        t = sl.text[1:]
    elif (sl.index_of_exclamation_mark is not None):
        t = sl.stmt[sl.index_of_exclamation_mark+1:]
    else:
        t = None
    if (t is not None):
        callback("//C%s" % t.expandtabs().rstrip())


def produce_comments(callback, ssl_list):
    """Emit all comments without filtering BLAS/LAPACK boilerplate.

    Boilerplate removal is now handled by a separate post-processing script,
    so this function simply converts Fortran comments to C++ line comments.
    """
    for ssl in ssl_list:
        if ssl is None:
            continue
        for sl in ssl.source_line_cluster:
            if sl.stmt_offs is None:
                # Full-line comment: take text after the leading "*"
                t = sl.text[1:]
            elif sl.index_of_exclamation_mark is not None:
                # In-line comment: take text after "!"
                t = sl.stmt[sl.index_of_exclamation_mark + 1:]
            else:
                t = None

            if t is not None:
                callback("//C%s" % t.expandtabs().rstrip())


def flush_comments_if_non_trivial(callback, buffer):
    for line in buffer:
        if (line != "//C"):
            for line in buffer:
                callback(line)
            return


def produce_leading_comments(callback, fproc):
    buffer = []
    produce_comments(callback=buffer.append, ssl_list=fproc.leading_comments)
    if fproc.top_ssl is not None:
        produce_comments(callback=buffer.append, ssl_list=[fproc.top_ssl])
    flush_comments_if_non_trivial(callback=callback, buffer=buffer)


def produce_trailing_comments(callback, fproc):
    buffer = []
    if fproc.end_ssl is not None:
        produce_comments(callback=buffer.append, ssl_list=[fproc.end_ssl])
    produce_comments(
        callback=buffer.append, ssl_list=fproc.trailing_comments)
    flush_comments_if_non_trivial(callback=callback, buffer=buffer)


class comment_manager(object):

    # Boilerplate filtering is done in an external post-processing step.
    # Keep this manager only as an ordered emitter of Fortran comment lines.
    __slots__ = ["sl_list", "index", "first_comment_output"]

    def __init__(O, fproc):
        O.sl_list = []
        for ssl in fproc.body_lines:
            if (ssl is not None):
                for sl in ssl.source_line_cluster:
                    O.sl_list.append(sl)
        O.sl_list.sort(key=lambda source_line: source_line.global_line_index)
        O.index = 0
        O.first_comment_output = False

    def produce(O, callback):
        produce_comment_given_sl(callback=callback, sl=O.sl_list[O.index])
        O.index += 1

    def insert_before(O, executable_info, callback):
        i = executable_info.ssl.source_line_cluster[-1].global_line_index
        # Add leading empty "//" line before the first comment block
        if (not O.first_comment_output and O.index != len(O.sl_list)):
            # Check if there are any comments to output
            has_comment_to_output = False
            for idx in range(O.index, len(O.sl_list)):
                j = O.sl_list[idx].global_line_index
                if (j > i):
                    break
                has_comment_to_output = True
                break
            if has_comment_to_output:
                callback("//")
                O.first_comment_output = True
        while (O.index != len(O.sl_list)):
            j = O.sl_list[O.index].global_line_index
            if (j > i):
                break
            O.produce(callback=callback)

    def flush_remaining(O, callback):
        while (O.index != len(O.sl_list)):
            O.produce(callback=callback)


class conv_hook_info(object):

    __slots__ = [
        "ignore_common_and_save",
        "needs_sve_dynamic_parameters",
        "variant_common_names",
        "needs_is_called_first_time",
        "needs_variant_bind",
        "data_init_after_variant_bind"]

    def __init__(O):
        O.ignore_common_and_save = False
        O.needs_sve_dynamic_parameters = False
        O.variant_common_names = None
        O.needs_is_called_first_time = None
        O.needs_variant_bind = None
        O.data_init_after_variant_bind = None


class global_conversion_info(object):

    __slots__ = [
        "topological_fprocs",
        "dynamic_parameters",
        "fortran_file_comments",
        "fem_do_safe",
        "arr_nd_size_max",
        "inline_all",
        "fprocs_by_name",
        "converted_commons_info",
        "separate_namespaces",
        "data_values_block_size",
        "data_specializations"]

    def __init__(O,
                 topological_fprocs,
                 dynamic_parameters,
                 fortran_file_comments,
                 fem_do_safe,
                 arr_nd_size_max,
                 inline_all,
                 converted_commons_info,
                 separate_namespaces,
                 data_values_block_size,
                 data_specializations):
        O.topological_fprocs = topological_fprocs
        O.dynamic_parameters = dynamic_parameters
        O.fortran_file_comments = fortran_file_comments
        O.fem_do_safe = fem_do_safe
        O.arr_nd_size_max = arr_nd_size_max
        O.inline_all = inline_all
        O.fprocs_by_name = topological_fprocs.all_fprocs.fprocs_by_name()
        O.converted_commons_info = converted_commons_info
        O.separate_namespaces = separate_namespaces
        O.data_values_block_size = data_values_block_size
        O.data_specializations = data_specializations

    def specialized(O, fproc):
        return conversion_info(global_conv_info=O, fproc=fproc)


class conversion_info(global_conversion_info):

    __slots__ = global_conversion_info.__slots__ + [
        "fproc",
        "comment_manager",
        "vmap",
        "data_initializers",
        "array_data_initializers",
        "hoisted_data_array_names",
        "ld_constant_decls",
    ]

    def __init__(O,
                 global_conv_info=None,
                 fproc=None,
                 vmap=None):
        val = None
        for slot in global_conversion_info.__slots__:
            if (global_conv_info is not None):
                val = getattr(global_conv_info, slot)
            setattr(O, slot, val)
        O.fproc = fproc
        if (O.fproc is None):
            O.comment_manager = None
        else:
            O.comment_manager = comment_manager(fproc=O.fproc)
        if (vmap is None):
            O.vmap = {}
        else:
            O.vmap = vmap

        # IMPORTANT: initialize DATA initializer map
        O.data_initializers = None

        # Map: array name (lower) -> (elem_ctype, dims_ints:list[int], init_list:list[str], rank:int)
        # Used to fold DATA statements into a single static array initializer.
        O.array_data_initializers = None

        # Set of array names (lower) that were hoisted and emitted as static initializers.
        O.hoisted_data_array_names = set()

        # Map: auto-generated leading-dimension variable -> initializer expression.
        # This avoids hardcoding constants and repeating expressions in flattened 2D indexing.
        O.ld_constant_decls = {}

    def set_vmap_force_local(O, fdecl):
        identifier = fdecl.id_tok.value
        O.vmap[identifier] = prepend_identifier_if_necessary(identifier)

    def set_vmap_for_callable(O, identifier):
        if (O.separate_namespaces is not None):
            ns = O.separate_namespaces.get(identifier)
            if (ns is not None):
                O.vmap[identifier] \
                    = ns + "::" + prepend_identifier_if_necessary(identifier)
                return True
        if (identifier in ["getargc", "iargc"]):
            O.vmap[identifier] = "cmn.%s" % identifier
            return True
        from fable import intrinsics
        if (identifier in intrinsics.extra_set_lower):
            O.vmap[identifier] = "fem::" + identifier
            return True
        return False

    def set_vmap_from_fdecl(O, fdecl):
        identifier = fdecl.id_tok.value
        if (fdecl.is_common()):
            if (getattr(fdecl, "dim_tokens", None) is None
                    and identifier.lower() in FABLE_EXTERN_COMMON_SCALARS):
                # COMMON scalar is provided externally (declared as an extern global).
                O.vmap[identifier] = prepend_identifier_if_necessary(identifier)
            else:
                O.vmap[identifier] = "cmn." + \
                    prepend_identifier_if_necessary(identifier)
        elif (fdecl.is_save()
              and not (O.fproc is not None
                       and getattr(O.fproc, "conv_hook", None) is not None
                       and O.fproc.conv_hook.ignore_common_and_save)):
            O.vmap[identifier] = "sve." + \
                prepend_identifier_if_necessary(identifier)
        elif (fdecl.is_intrinsic()):
            low = identifier.lower()
            if (low in ["float", "int", "char"]):
                # Normalize case to avoid fem::int vs fem::fint divergence.
                O.vmap[identifier] = "fem::f" + low
            elif (identifier == "iargc"):
                O.vmap[identifier] = "cmn.iargc"
            elif (identifier == "ceiling"):
                O.vmap[identifier] = "iceil"
            else:
                O.vmap[identifier] = "fem::" + identifier
        elif (not O.set_vmap_for_callable(identifier=fdecl.id_tok.value)):
            O.set_vmap_force_local(fdecl=fdecl)
            return False
        return True

    def vmapped(O, fdecl):
        identifier = fdecl.id_tok.value
        result = O.vmap.get(identifier)
        if (result is None):
            O.set_vmap_from_fdecl(fdecl=fdecl)
            result = O.vmap[identifier]
        return result

    def vmapped_callable(O, identifier):
        result = O.vmap.get(identifier)
        if (result is None):
            if (not O.set_vmap_for_callable(identifier=identifier)):
                O.vmap[identifier] = prepend_identifier_if_necessary(
                    identifier)
            result = O.vmap[identifier]
        return result


def called_fproc_needs_cmn(conv_info, called_name):
    # MPLAPACK customization:
    # Do NOT thread `common& cmn` through generated function signatures.
    # The project rewrites COMMON/SAVE handling, so we also suppress the
    # auto-injection of `cmn` as the first actual argument at call sites.
    return False


def cmn_needs_to_be_inserted(conv_info, prev_tok):
    # prev_tok may legitimately be None
    if prev_tok is None:
        return False

    # Only an identifier token can be a procedure name eligible for cmn injection.
    # Be defensive: some token kinds carry list-valued .value, which must not
    # be passed to .lower().
    try:
        if not prev_tok.is_identifier():
            return False
    except Exception:
        return False

    if not isinstance(getattr(prev_tok, "value", None), str):
        return False

    # Intrinsic functions (e.g. HUGE) are NOT variables and must not be
    # queried via get_fdecl().
    from fable import intrinsics
    lname = prev_tok.value.lower()
    if (lname in intrinsics.set_lower
        or lname in intrinsics.extra_set_lower
            or lname in intrinsics.io_set_lower):
        return False

    if (conv_info.fprocs_by_name is not None
            and conv_info.fproc is not None):
        # Some F90 sources call helper procedures brought in via USE
        # statements (e.g. LA_ISNAN) that we may strip during preprocessing.
        # Do not crash if the identifier has no declaration record.
        try:
            fdecl = conv_info.fproc.get_fdecl(id_tok=prev_tok)
        except KeyError:
            return False

        if (fdecl.is_user_defined_callable()):
            if (called_fproc_needs_cmn(
                conv_info=conv_info,
                    called_name=prev_tok.value)):
                return True
    return False


def convert_power(conv_info, tokens):
    # Special-case INTEGER power-of-two:
    #   2 ** k  -> (INTEGER(1) << k)
    # This matches Fortran INTEGER exponentiation and avoids floating pow().
    if len(tokens) == 2:
        base_tok = tokens[0]
        exp_tok = tokens[1]
        if (base_tok is not None and base_tok.is_integer()
                and base_tok.value == "2"):
            exp_str = convert_tokens(
                conv_info=conv_info, tokens=[exp_tok], commas=False).strip()

            def _decl_is_integer(name: str) -> bool:
                if conv_info is None or getattr(conv_info, "fproc", None) is None:
                    return False
                fdecl_map = getattr(
                    conv_info.fproc, "fdecl_by_identifier", None)
                if not fdecl_map:
                    return False
                fd = fdecl_map.get(name.lower()) or fdecl_map.get(name)
                if fd is None:
                    return False
                dt = getattr(fd, "data_type", None)
                code = dt if isinstance(
                    dt, str) else getattr(dt, "value", None)
                return (code or "").lower() == "integer"

            def _decl_exists(name: str) -> bool:
                if conv_info is None or getattr(conv_info, "fproc", None) is None:
                    return False
                fdecl_map = getattr(
                    conv_info.fproc, "fdecl_by_identifier", None) or {}
                return (name.lower() in fdecl_map) or (name in fdecl_map)

            is_integer_exp = False
            if exp_tok.is_integer():
                is_integer_exp = True
            elif exp_tok.is_identifier():
                is_integer_exp = _decl_is_integer(exp_tok.value)
            else:
                # Heuristic: if the exponent expression references only INTEGER
                # declared variables (ignore unknown identifiers/functions),
                # treat it as INTEGER.
                ids = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", exp_str))
                any_decl = False
                ok = True
                for name in ids:
                    if _decl_is_integer(name):
                        any_decl = True
                    elif _decl_exists(name):
                        # Declared but not INTEGER -> reject
                        ok = False
                        break
                # Avoid UB for obvious negative shift like "-k"
                if ok and (any_decl or exp_str.isdigit()) and exp_str and not exp_str.lstrip().startswith("-"):
                    is_integer_exp = True

            if is_integer_exp:
                return f"(INTEGER(1) << ({exp_str}))"

    fun = "fem::pow"
    pow_tok = tokens[1]
    if (pow_tok.is_integer()):
        if (pow_tok.value == "1"):
            fun = ""
            tokens = tokens[:1]
        elif (pow_tok.value in ["2", "3", "4"]):
            fun = "fem::pow%s" % pow_tok.value
            tokens = tokens[:1]
    return fun + "(" + convert_tokens(
        conv_info=conv_info, tokens=tokens, commas=True) + ")"


def _get_return_kind_from_signatures(base_cpp_name: str) -> typing.Optional[str]:
    """Return coarse return-type category using FUNCTION_RETURNS.

    base_cpp_name: bare C++ identifier (no namespace), any case.
    Result is one of: "real", "complex", "integer", "logical", "void", or None.
    """
    if not FUNCTION_RETURNS:
        return None

    name = base_cpp_name.lower()
    # Prefer direct lookup by the C++-side key (what gen_mplapack_signatures.py emits).
    ret = FUNCTION_RETURNS.get(name)
    if ret is None:
        # Fallback: try mapped MPLAPACK name.
        mapped = convert_function_name_to_mplapack(name)
        if mapped:
            ret = FUNCTION_RETURNS.get(str(mapped).lower())
    if ret is None:
        # Last resort: try reverse-mapped Fortran name.
        f_name = _MPLAPACK_CPP_TO_FORTRAN.get(name, name)
        ret = FUNCTION_RETURNS.get(
            f_name) or FUNCTION_RETURNS.get(str(f_name).lower())
    if ret is None:
        return None

    ret = ret.lower()
    if "complex" in ret:
        return "complex"
    if "doubleprecision" in ret or "double precision" in ret:
        return "real"
    if "double" in ret or "real" in ret:
        return "real"
    if "integer" in ret:
        return "integer"
    if "logical" in ret:
        return "logical"
    if "void" in ret:
        return "void"
    return None


def _is_real_valued_function_call(expr: str) -> bool:
    """Return True if expr syntactically looks like a call to a REAL-valued function.

    This is a conservative heuristic used only in coercion logic
    (REAL/DBLE/INT and REAL-LHS assignment fixes). It must never
    claim that a COMPLEX-valued function is REAL.
    """
    s = expr.strip()
    if "(" not in s:
        return False

    # Match leading function name (optionally namespaced), e.g.:
    #   fem::dble(...), Clangb(...), max(...), cabs1(...)
    m = re.match(r"([A-Za-z_][A-Za-z0-9_:]*)\s*\(", s)
    if not m:
        return False

    fname_cpp = m.group(1)
    # Strip namespaces like "fem::aimag"
    base = fname_cpp.split("::")[-1].lower()

    # Obvious real-valued intrinsics / helpers
    if base in ("abs", "fabs", "dabs", "sabs"):
        return True
    if base in ("aimag", "imag", "dimag"):
        # AIMAG/IMAG return REAL even if their argument is COMPLEX
        return True

    # Known real-valued helpers that take COMPLEX arguments
    # but return REAL (matrix norms, 1-norm, etc.).
    if base in ("cabs1", "rcabs1", "rlapy2", "rlapy3"):
        return True

    # Use FUNCTION_RETURNS first if available.
    ret_kind = _get_return_kind_from_signatures(
        base) if FUNCTION_RETURNS else None
    if ret_kind == "real":
        return True
    if ret_kind == "complex":
        return False

    # Fallback: try to map back to the original Fortran-style name
    # and use the usual BLAS/LAPACK prefix heuristic only if we have it.
    f_name = _MPLAPACK_CPP_TO_FORTRAN.get(base, base)
    if not f_name:
        return False

    # BLAS/LAPACK convention (very rough):
    #   s* / d*  → real-valued
    #   c* / z*  → complex-valued
    first = f_name[0].lower()
    if first in ("s", "d"):
        return True
    if first in ("c", "z"):
        return False

    # Unknown prefix: do not assume real-valued
    return False


def _is_simple_lvalue(expr: str) -> bool:
    """Return True if expr is a simple variable or an array element.

    Examples considered simple:
        A
        A[i]
        A[i][j]
        A[(i-1) + (j-1)*lda]
    """
    expr = expr.strip()
    # NAME or NAME[...] or NAME[...]...
    # This is intentionally permissive: any NAME followed by bracketed expressions.
    return re.match(r'^[A-Za-z_][A-Za-z0-9_]*(\[[^\]]*\])*$', expr) is not None


def _rewrite_small_char_substrings(text: str) -> str:
    """Rewrite substring-like calls on small CHARACTER scalars:

       jbcmpz(1, 1) = 'S'  ->  jbcmpz[0] = 'S';
       jbcmpz(2, 2) = 'V'  ->  jbcmpz[1] = 'V';
    """
    if not small_char_identifiers:
        return text

    out = text

    # 1) JBCMPZ(1,1) -> jbcmpz[0]
    for name in small_char_identifiers:
        for idx in range(1, 5):
            pattern = r"\b%s\s*\(\s*%d\s*,\s*%d\s*\)" % (re.escape(name), idx, idx)
            repl = f"{name}[{idx - 1}]"
            out = re.sub(pattern, repl, out)

    # 2) jbcmpz[0] = "S" -> jbcmpz[0] = 'S'
    if small_char_identifiers:
        names_alt = "|".join(re.escape(n) for n in small_char_identifiers)
        pattern = r"(\b(?:%s)\s*\[\d+\]\s*=\s*)\"(.)\"" % names_alt

        def _to_char_literal(m):
            left = m.group(1)
            ch = m.group(2)
            return f"{left}'{ch}'"

        out = re.sub(pattern, _to_char_literal, out)

    return out


def _emit_small_char_string_assignment(curr_scope, clhs: str, crhs: str, conv_info=None) -> bool:
    """Emit element-wise assignments for small CHARACTER*n scalars mapped to char[].

    This is needed because plain C arrays cannot be assigned in C++:
        char c2[2];
        c2 = path(2,3);   // invalid

    We recognize the converted fixed-form substring pattern:
        dst = src(start, end)
    where both dst/src are in small_char_identifier_lengths / small_char_identifiers.

    For constant start/end (integer literals), we emit the exact indices:
        c2[0] = path[(2 - 1)];
        c2[1] = path[(3 - 1)];

    Returns True if we emitted code and the caller should NOT emit "dst = rhs;".
    """
    dst = clhs.strip()
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", dst):
        return False

    dst_len = small_char_identifier_lengths.get(dst)
    if not isinstance(dst_len, int) or dst_len <= 0:
        return False

    rhs = crhs.strip()

    # Parse "src(...)" at the top level.
    m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*(.*)\s*\)\s*$", rhs)
    if m is None:
        return False

    src = m.group(1)
    if src not in small_char_identifiers:
        # Allow substring copies from dummy CHARACTER arguments modeled as plain char*.
        if conv_info is None or src not in _plain_char_ptr_dummy_cpp_names(conv_info):
            return False

    inner = m.group(2)
    args = _split_top_level_commas(inner)
    if len(args) != 2:
        return False

    start_expr = args[0].strip()
    end_expr = args[1].strip()

    def _try_parse_int_literal(expr: str):
        mm = re.match(r"^\(?\s*(\d+)\s*\)?$", expr)
        if not mm:
            return None
        try:
            return int(mm.group(1))
        except ValueError:
            return None

    start_i = _try_parse_int_literal(start_expr)
    end_i = _try_parse_int_literal(end_expr)

    # Constant range: preserve the original indices (good for readability and diffs).
    if start_i is not None and end_i is not None and end_i >= start_i:
        avail = end_i - start_i + 1
        ncopy = min(dst_len, avail)
        for i in range(ncopy):
            pos = start_i + i  # Fortran 1-based
            curr_scope.append(f"{dst}[{i}] = {src}[({pos} - 1)];")
        # Fortran semantics: pad with spaces if rhs is shorter than lhs.
        for i in range(ncopy, dst_len):
            curr_scope.append(f"{dst}[{i}] = ' ';")
        return True

    # Generic range: copy up to dst_len characters, padding with spaces
    # when the substring is shorter (end < start + i).
    for i in range(dst_len):
        if i == 0:
            pos = f"({start_expr})"
        else:
            pos = f"(({start_expr}) + {i})"
        cond = f"({pos}) <= ({end_expr})"
        idx0 = f"(({pos}) - 1)"
        curr_scope.append(f"{dst}[{i}] = ({cond}) ? {src}[{idx0}] : ' ';")
    return True


def _has_top_level_arith_op(expr: str) -> bool:
    """Return True if expr has a top-level +, -, * or / outside parentheses.

    This is used to conservatively decide when an expression is a "real"
    arithmetic combination that may need REAL = COMPLEX coercion.
    """
    s = expr
    depth = 0
    in_squote = False
    in_dquote = False
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "'" and not in_dquote:
            in_squote = not in_squote
        elif ch == '"' and not in_squote:
            in_dquote = not in_dquote
        elif not in_squote and not in_dquote:
            if ch == '(':
                depth += 1
            elif ch == ')':
                if depth > 0:
                    depth -= 1
            elif depth == 0 and ch in "+-*/":
                # We do not try to distinguish unary vs binary here.
                # Any top-level arithmetic suggests a numeric expression.
                return True
        i += 1
    return False


def _extract_base_identifier_from_unary(expr: str):
    """Extract base identifier from unary expression like '-alpha' or '-x[i]'.

    Returns (sign, base_identifier, full_base_expr) or (None, None, None) if not matching.

    Examples:
        '-alpha'    -> ('-', 'alpha', 'alpha')
        '-x[i]'     -> ('-', 'x', 'x[i]')
        '+beta'     -> ('+', 'beta', 'beta')
        'gamma'     -> (None, 'gamma', 'gamma')  # No unary operator
        'a + b'     -> (None, None, None)  # Complex expression
    """
    expr = expr.strip()

    # Check for unary operator at the start
    sign = None
    if expr.startswith('-') or expr.startswith('+'):
        sign = expr[0]
        expr = expr[1:].strip()

    # Now check if remaining is a simple lvalue
    m = re.match(r'^([A-Za-z_][A-Za-z0-9_]*)(\[[^\]]*\])*$', expr)
    if m:
        return (sign, m.group(1), expr)

    return (None, None, None)


def _collect_complex_variable_names(text: str):
    """Collect names of COMPLEX variables and pointers from the C++ translation unit.

    Returns:
        complex_vars: set of variable names that are COMPLEX (scalars or pointers)
        complex_ptrs: subset of complex_vars that are declared as pointers.
    """
    complex_vars = set()
    complex_ptrs = set()

    # Pointer-style declarations and arguments: COMPLEX *a, COMPLEX* a, etc.
    for name in re.findall(r'\bCOMPLEX\s*\*+\s*([A-Za-z_][A-Za-z0-9_]*)', text):
        complex_vars.add(name)
        complex_ptrs.add(name)

    # Scalar or reference-style: COMPLEX alpha; COMPLEX const &alpha; etc.
    for name in re.findall(r'\bCOMPLEX\s+(?:const\s+&\s*)?([A-Za-z_][A-Za-z0-9_]*)', text):
        complex_vars.add(name)

    return complex_vars, complex_ptrs


def _mask_real_valued_subcalls(expr: str) -> str:
    """Replace real-valued function calls in expr with a neutral placeholder.

    This is only used when deciding whether an expression should be treated
    as COMPLEX-valued for REAL/DBLE/INT coercions.
    """
    s = expr
    out = []
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch.isalpha() or ch == "_":
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] in ("_", ":")):
                j += 1
            name = s[i:j]
            k = j
            while k < n and s[k].isspace():
                k += 1
            if k < n and s[k] == "(":
                depth = 0
                end = k
                while end < n:
                    c = s[end]
                    if c == "(":
                        depth += 1
                    elif c == ")":
                        depth -= 1
                        if depth == 0:
                            end += 1
                            break
                    end += 1
                call_str = s[i:end]
                if _is_real_valued_function_call(call_str):
                    # Mask the whole call as a scalar REAL.
                    out.append("0")
                    i = end
                    continue
            out.append(ch)
            i += 1
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def _contains_complex_valued_function_call(expr: str) -> bool:
    """Return True if expr contains a function call known to return COMPLEX.

    This is used to avoid treating expressions as 'obviously real' when they
    still depend on complex-valued subcalls such as CDOTC/ZDOTC.
    """
    s = expr
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch.isalpha() or ch == "_":
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] in ("_", ":")):
                j += 1
            name = s[i:j]
            k = j
            while k < n and s[k].isspace():
                k += 1
            if k < n and s[k] == "(":
                # We found a call; find its end to skip over it.
                depth = 0
                end = k
                while end < n:
                    c = s[end]
                    if c == "(":
                        depth += 1
                    elif c == ")":
                        depth -= 1
                        if depth == 0:
                            end += 1
                            break
                    end += 1
                base = name.split("::")[-1].lower()
                kind = _get_return_kind_from_signatures(
                    base) if FUNCTION_RETURNS else None
                if kind == "complex":
                    return True
                i = end
                continue
        i += 1
    return False


def _looks_complex_expression(arg: str) -> bool:
    """Heuristic: return True if 'arg' syntactically looks COMPLEX-valued.

    Used only for explicit REAL/DBLE/INT handling, when no known COMPLEX
    variable name appears in the expression.
    """
    s = arg.strip()

    # Explicit COMPLEX constructor: COMPLEX(...)
    if re.search(r'\bCOMPLEX\s*\(', s):
        return True

    # std::complex<T>(...) temporaries
    if "std::complex<" in s:
        return True

    # Conjugation intrinsics: conj(...), conjg(...)
    if re.search(r'\bconj\s*\(', s):
        return True
    if re.search(r'\bconjg\s*\(', s):
        return True

    return False


def _real_cast_or_component(arg, complex_identifiers, complex_pointer_identifiers):
    """Implement REAL/DBLE on an arbitrary expression.

    - If the expression is COMPLEX-valued (contains COMPLEX vars or
      COMPLEX-valued function calls), take its real part.
    - Otherwise, treat it as a REAL expression and wrap in castREAL(...).
    """
    arg = arg.strip()
    if not arg:
        return "castREAL(0)"

    # ------------------------------------------------------------------
    # 1) Unary minus on a simple COMPLEX lvalue:
    #    REAL(-alpha) -> -alpha.real()
    #    REAL(-a(i))  -> -a(i).real()
    # ------------------------------------------------------------------
    m_unary = re.match(r'^-\s*([A-Za-z_][A-Za-z0-9_]*(?:\[[^\]]*\])*)$', arg)
    if m_unary:
        base_expr = m_unary.group(1)
        m_name = re.match(r'([A-Za-z_][A-Za-z0-9_]*)', base_expr)
        base = m_name.group(1) if m_name else None
        if base is not None and base in complex_identifiers:
            # COMPLEX* dummy used bare: -a -> -a[0].real()
            if base in complex_pointer_identifiers and "[" not in base_expr:
                return f"-{base}[0].real()"
            # General COMPLEX variable / element: -a(i) -> -a(i).real()
            return f"-{base_expr}.real()"

    # ------------------------------------------------------------------
    # 2) Single function call: use FUNCTION_RETURNS if available.
    #    This covers patterns like DBLE(ZDOTC(...)).
    # ------------------------------------------------------------------
    m_call = re.match(r'^([A-Za-z_][A-Za-z0-9_:]*)\s*\(.*\)\s*$', arg)
    if m_call:
        fname_cpp = m_call.group(1)
        base = fname_cpp.split("::")[-1].lower()
        ret_kind = _get_return_kind_from_signatures(base)
        if ret_kind == "complex":
            # COMPLEX-valued function -> take real part
            #   DBLE(ZDOTC(...)) -> Cdotc(...).real()
            return f"{arg}.real()"
        if ret_kind == "real":
            # REAL-valued function -> keep as real
            return f"castREAL({arg})"
        # If we cannot classify the function, fall through to generic logic.

    # ------------------------------------------------------------------
    # 3) Expressions that are obviously REAL already
    # ------------------------------------------------------------------
    if _is_real_valued_function_call(arg):
        # e.g. DBLE(Rlapy2(...)), cabs1(z), Clangb(...)
        return f"castREAL({arg})"

    # If the expression already uses .real()/.imag() or castREAL,
    # just treat it as real-valued.
    if ".real()" in arg or ".imag()" in arg or "castREAL(" in arg:
        return f"castREAL({arg})"

    # ------------------------------------------------------------------
    # 4) Does the expression contain any known COMPLEX variables?
    # ------------------------------------------------------------------
    complex_vars = set(complex_identifiers) | set(complex_pointer_identifiers)
    owner = None
    for name in complex_vars:
        if arg == name or arg.startswith(name + "[") or arg.startswith(name + "."):
            owner = name
            break
        if re.search(r"\b" + re.escape(name) + r"\b", arg):
            owner = name
            break

    # ------------------------------------------------------------------
    # 5) No known COMPLEX vars: maybe still a COMPLEX expression?
    #    (e.g. COMPLEX(...), conj(z), etc.)
    # ------------------------------------------------------------------
    if owner is None:
        if _looks_complex_expression(arg):
            # If this is a single function call or a simple lvalue,
            # avoid extra parentheses:
            #   func(...).real()       instead of (func(...)).real()
            #   work[i-1].real()       instead of (work[i-1]).real()
            if (re.match(r'^[A-Za-z_][A-Za-z0-9_:]*\s*\(.*\)\s*$', arg)
                    or _is_simple_lvalue(arg)):
                return f"{arg}.real()"
            return f"({arg}).real()"
        # Otherwise: treat as real and cast if needed.
        return f"castREAL({arg})"

    # ------------------------------------------------------------------
    # 6) Expression contains at least one COMPLEX variable
    #    -> take its real part.
    # ------------------------------------------------------------------
    # For pure function calls or simple lvalues, avoid redundant parentheses.
    if (re.match(r'^[A-Za-z_][A-Za-z0-9_:]*\s*\(.*\)\s*$', arg)
            or _is_simple_lvalue(arg)):
        return f"{arg}.real()"
    return f"({arg}).real()"


def _imag_repl(arg: str) -> str:
    """Replacement for AIMAG/IMAG: choose arg.imag() or (arg).imag()."""
    arg = arg.strip()
    if _is_simple_lvalue(arg):
        return f"{arg}.imag()"
    return f"({arg}).imag()"


def _conj_repl(arg: str) -> str:
    """Replacement for CONJG/CONJ: std::conj(arg)."""
    arg = arg.strip()
    return f"conj({arg})"


def _is_intrinsic_name(name: str) -> bool:
    """Return True if 'name' should be treated as a Fortran intrinsic."""
    from fable import intrinsics
    lname = name.lower()
    return (
        lname in intrinsics.set_lower
        or lname in intrinsics.extra_set_lower
        or lname in intrinsics.io_set_lower
    )


def _map_intrinsic_vmap(conv_info, name: str) -> None:
    """Map an intrinsic identifier into conv_info.vmap as fem::<name>."""
    lname = name.lower()
    conv_info.vmap[name] = "fem::" + lname
    conv_info.vmap[lname] = "fem::" + lname


def _huge_repl(_arg: str) -> str:
    """Replacement for HUGE(x): use LAPACK-style overflow constant via Rlamch."""
    # We intentionally ignore the argument and use the machine overflow constant.
    return 'Rlamch("O")'


_dummy_character_args = set()


def _mark_dummy_character_arg(conv_info, name: str) -> None:
    """Remember that (current fproc, name) is a CHARACTER dummy argument."""
    key = (id(conv_info.fproc), name.lower())
    _dummy_character_args.add(key)


def _is_dummy_character_arg(conv_info, name: str) -> bool:
    """Return True if (current fproc, name) is a CHARACTER dummy argument."""
    key = (id(conv_info.fproc), name.lower())
    return key in _dummy_character_args

def _is_fem_str_scalar(conv_info, name: str) -> bool:
    """Return True if 'name' is a scalar CHARACTER mapped to fem::str<N>."""
    if conv_info is None or getattr(conv_info, "fproc", None) is None:
        return False
    if _is_dummy_character_arg(conv_info, name):
        return False
    if name in small_char_identifiers:
        return False
    try:
        fdecl = conv_info.fproc.fdecl_by_identifier.get(name.lower()) or conv_info.fproc.fdecl_by_identifier.get(name)
    except Exception:
        return False
    if fdecl is None:
        return False
    dt = getattr(fdecl, "data_type", None)
    dt_code = dt if isinstance(dt, str) else getattr(dt, "value", None)
    if str(dt_code).lower() != "character":
        return False
    if getattr(fdecl, "dim_tokens", None) is not None:
        return False
    st = getattr(fdecl, "size_tokens", None)
    if st is None:
        return False
    try:
        if len(st) == 1 and getattr(st[0], "is_integer", None) and st[0].is_integer() and st[0].value == "1":
            return False
    except Exception:
        pass
    return True

def _is_plain_character_pointer_dummy(conv_info, name: str) -> bool:
    """Return True if 'name' is a CHARACTER dummy argument emitted as (const) char*.

    We treat these as plain pointers:
      - scalar CHARACTER dummies (no dim_tokens)
      - CHARACTER arrays with element length 1 (CHARACTER*1 X(*))

    We exclude CHARACTER arrays that are emitted as str_arr_*ref (length > 1).
    """
    if conv_info is None or getattr(conv_info, "fproc", None) is None:
        return False
    if not _is_dummy_character_arg(conv_info, name):
        return False
    try:
        fdecl = conv_info.fproc.fdecl_by_identifier.get(name.lower())
        if fdecl is None:
            fdecl = conv_info.fproc.fdecl_by_identifier.get(name)
    except Exception:
        return False
    if fdecl is None:
        return False
    dt = getattr(fdecl, "data_type", None)
    dt_code = getattr(dt, "value", None) if dt is not None else None
    if dt_code != "character":
        return False
    # Scalar CHARACTER dummy: emitted as (const) char* unless view mode is enabled.
    if getattr(fdecl, "dim_tokens", None) is None:
        return (not FABLE_SMALL_CHAR_VIEW)
    # Array CHARACTER dummy: only CHARACTER*1 arrays are treated as plain pointer.
    st = getattr(fdecl, "size_tokens", None)
    if st is None:
        return True
    try:
        if len(st) == 1 and st[0].is_integer() and st[0].value == "1":
            return True
    except Exception:
        pass
    return False

def _is_scalar_character_fem_str(conv_info, name: str) -> bool:
    """True if 'name' is scalar CHARACTER emitted as fem::str<...> (not char/char[]/dummy)."""
    if conv_info is None or getattr(conv_info, "fproc", None) is None:
        return False
    if _is_dummy_character_arg(conv_info, name):
        return False
    if name in small_char_identifiers:
        return False
    try:
        fdecl = conv_info.fproc.fdecl_by_identifier.get(name.lower()) or \
                conv_info.fproc.fdecl_by_identifier.get(name)
    except Exception:
        return False
    if fdecl is None:
        return False
    dt = getattr(fdecl, "data_type", None)
    dt_code = dt if isinstance(dt, str) else getattr(dt, "value", None)
    if dt_code != "character":
        return False
    if getattr(fdecl, "dim_tokens", None) is not None:
        return False
    st = getattr(fdecl, "size_tokens", None)
    # In view mode, CHARACTER*1 scalars are emitted as fem::str<1> (not plain char).
    if st is None:
        return bool(FABLE_SMALL_CHAR_VIEW)
    try:
        if len(st) == 1 and st[0].is_integer() and st[0].value == "1":
            return bool(FABLE_SMALL_CHAR_VIEW)
    except Exception:
        pass
    try:
        size_expr = convert_tokens(conv_info=conv_info, tokens=st, commas=False).strip()
        if size_expr == "1":
            return bool(FABLE_SMALL_CHAR_VIEW)
    except Exception:
        pass
    return True

def _plain_char_ptr_dummy_cpp_names(conv_info):
    """Return a cached set of C++ names for CHARACTER dummy arguments modeled as plain (const) char*."""
    if conv_info is None or getattr(conv_info, "fproc", None) is None:
        return set()
    names = set()
    for id_tok in getattr(conv_info.fproc, "args", []) or []:
        try:
            if getattr(id_tok, "value", None) == "*":
                continue
            if _is_plain_character_pointer_dummy(conv_info, id_tok.value):
                cpp_name = conv_info.vmap.get(
                    id_tok.value, prepend_identifier_if_necessary(id_tok.value))
                names.add(cpp_name)
        except Exception:
            continue
    return names

def _rewrite_plain_char_ptr_unit_substrings(text: str, conv_info) -> str:
    """Rewrite Fortran-style unit-length substrings on plain char* dummies.

    Example:
        path(1, 1)  ->  path[(1) - 1]

    This is only applied to dummy CHARACTER arguments emitted as (const) char*.
    We intentionally only rewrite unit-length substrings (start == end) because
    longer substrings generally require an explicit temporary buffer.
    """
    names = _plain_char_ptr_dummy_cpp_names(conv_info)
    if not names:
        return text

    int_re = re.compile(r'^[+-]?[0-9]+$')
    name_re = re.compile(r'^[A-Za-z_]\w*$')

    out = text
    for name in names:
        # Match: name(arg1, arg2) where arg1/arg2 have no nested parentheses.
        pat = re.compile(rf"\b{re.escape(name)}\s*\(\s*([^,()]+?)\s*,\s*([^,()]+?)\s*\)")

        def _repl(m):
            a = m.group(1).strip()
            b = m.group(2).strip()
            if a.replace(" ", "") != b.replace(" ", ""):
                return m.group(0)
            # Prefer clean C++ for simple cases:
            #   path(1,1) -> path[0]
            #   path(i,i) -> path[i - 1]
            a_key = a.replace(" ", "")
            if int_re.fullmatch(a_key):
                try:
                    return f"{name}[{int(a_key) - 1}]"
                except Exception:
                    # Fallback (should be rare)
                    return f"{name}[{a_key} - 1]"
            if name_re.fullmatch(a_key):
                return f"{name}[{a_key} - 1]"
            # Complex expression: keep parentheses for safety
            return f"{name}[({a}) - 1]"

        out = pat.sub(_repl, out)
    return out

def _rewrite_unary_intrinsic(text: str, func_name: str, repl_func):
    """Rewrite occurrences of func_name(arg) using repl_func(arg).

    This handles nested parentheses in arg by scanning until matching ')'.
    """
    out = []
    i = 0
    n = len(text)
    fname_len = len(func_name)

    while i < n:
        j = text.find(func_name, i)
        if j < 0:
            out.append(text[i:])
            break

        out.append(text[i:j])

        # Find '(' after function name
        k = text.find("(", j + fname_len)
        if k < 0:
            # Malformed: stop rewriting
            out.append(text[j:])
            break

        # Scan to matching ')'
        depth = 0
        p = k
        while p < n:
            c = text[p]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            p += 1
        if depth != 0:
            # Unbalanced; give up
            out.append(text[j:])
            break

        # Extract argument
        arg = text[k + 1:p]
        out.append(repl_func(arg))

        i = p + 1

    return "".join(out)


def _rewrite_max_min_calls(text: str) -> str:
    """Rewrite fem::max/min calls: drop fem:: and cast integer literal arguments.

    Examples:
      fem::max(1, m)   -> max((INTEGER)1, m)
      fem::min(-2, n)  -> min((INTEGER)-2, n)
      max(m, 1)        -> max(m, (INTEGER)1)
    """
    # First drop fem:: namespace for max/min.
    # This is cheap even if not present.
    text = text.replace("fem::max", "max")
    text = text.replace("fem::min", "min")

    # Case 1: integer literal as FIRST argument
    #   max(1, m)   -> max((INTEGER)1, m)
    #   min(-2, n)  -> min((INTEGER)-2, n)
    pattern_first = r'\b(max|min)\(\s*([+-]?[0-9]+)\s*,'

    def repl_first(m):
        func = m.group(1)
        lit = m.group(2)
        return f"{func}((INTEGER){lit}, "

    text = re.sub(pattern_first, repl_first, text)

    # Case 2: integer literal as SECOND argument
    #   max(m, 1)   -> max(m, (INTEGER)1)
    #   min(k, -2)  -> min(k, (INTEGER)-2)
    pattern_second = r'\b(max|min)\(\s*([^,]+?),\s*([+-]?[0-9]+)\s*\)'

    def repl_second(m):
        func = m.group(1)
        first = m.group(2)
        lit = m.group(3)
        return f"{func}({first}, (INTEGER){lit})"

    text = re.sub(pattern_second, repl_second, text)

    return text


_single_char_string_assign_re = re.compile(
    r'(\b[A-Za-z_][A-Za-z0-9_]*\b(?:\s*\[\s*\d+\s*\])?)\s*=\s*"([^"\\])"\s*;'
)
_single_char_string_eq_re = re.compile(
    r'(\b[A-Za-z_][A-Za-z0-9_]*\b)\s*==\s*"([^"\\])"'
)
_single_char_string_ne_re = re.compile(
    r'(\b[A-Za-z_][A-Za-z0-9_]*\b)\s*!=\s*"([^"\\])"'
)


def rewrite_single_char_string_literals(line: str) -> str:
    """Rewrite x = "A"; into x = 'A'; for single-character variables.

    Comparisons (x == "A", x != "A") are left as-is, because in the
    translated code such variables are often const char* coming from
    Fortran CHARACTER*1 and should be handled via Mlsame or explicit
    dereference rather than a naive char literal compare.
    """
    line = _single_char_string_assign_re.sub(r"\1 = '\2';", line)
    return line


def get_lower_bound(conv_info, fdecl, dim_index):
    """Return the lower bound expression (as C++ string) for a given array dimension.

    For a declaration like A(L:U), this returns convert_tokens(L).
    If no explicit lower bound is present (e.g. A(N) or A(U)), this returns "1"
    which corresponds to the usual Fortran 1-based convention.
    """
    dt = getattr(fdecl, "dim_tokens", None)
    if dt is None:
        return "1"
    if dim_index < 0 or dim_index >= len(dt):
        return "1"

    # Each item in dim_tokens has a .value field which is
    # a token list for that dimension.
    tokens = dt[dim_index].value

    # Look for explicit "lower:upper" form.
    for i, tok in enumerate(tokens):
        if tok.is_op_with(value=":"):
            # ":upper" -> default lower bound = 1
            if i == 0:
                return "1"
            # Convert lower-bound tokens using the normal token converter.
            lb_code = convert_tokens(conv_info=conv_info, tokens=tokens[:i])
            return lb_code.strip()

    # No ":" found -> implicit 1-based lower bound.
    return "1"


def _try_extract_first_dim_extent_identifier(conv_info, dim0_tokens):
    """Try to return a simple identifier for the extent of the first dimension.

    This avoids emitting verbose extents like '((ldb - 1) - 0 + 1)' when it is
    obviously equal to 'ldb' for common Fortran bounds such as '0:ldb-1' or '1:ldb'.
    This helper is intentionally conservative.
    """
    # Find "lower:upper" split.
    colon_idx = None
    for i, tok in enumerate(dim0_tokens):
        if hasattr(tok, "is_op_with") and tok.is_op_with(value=":"):
            colon_idx = i
            break
    if colon_idx is None:
        return None

    lb_tokens = dim0_tokens[:colon_idx]
    ub_tokens = dim0_tokens[colon_idx + 1:]
    if not ub_tokens:
        return None

    def _is_int_literal(toks, v):
        return (len(toks) == 1 and hasattr(toks[0], "is_integer")
                and toks[0].is_integer() and toks[0].value == str(v))

    def _is_identifier(toks):
        return (len(toks) == 1 and hasattr(toks[0], "is_identifier")
                and toks[0].is_identifier())

    # Fortran default lower bound is 1 if omitted.
    if not lb_tokens:
        lb = 1
    elif _is_int_literal(lb_tokens, 0):
        lb = 0
    elif _is_int_literal(lb_tokens, 1):
        lb = 1
    else:
        return None

    # Pattern: 1:ID  -> extent = ID
    if lb == 1 and _is_identifier(ub_tokens):
        return convert_tokens(conv_info=conv_info, tokens=ub_tokens).strip()

    # Pattern: 0:(ID-1) -> extent = ID
    if (lb == 0 and len(ub_tokens) == 3
            and ub_tokens[0].is_identifier()
            and ub_tokens[1].is_op_with(value="-")
            and ub_tokens[2].is_integer()
            and ub_tokens[2].value == "1"):
        return convert_tokens(conv_info=conv_info, tokens=[ub_tokens[0]]).strip()

    return None


def get_leading_dimension_expr(conv_info, fdecl, default=None):
    """Return leading dimension (extent of the first dimension) for a 2D array.

    Prefer the declaration's first dimension, e.g.:
        DOUBLE PRECISION A(LDA, * )      -> "lda"
        COMPLEX*16 WORK13(LDWORK, NBMAX) -> "ldwork"

    If the first dimension is explicit-bounds (L:U), return (U - L + 1).
    For common bounds (0:ld-1, 1:ld), return just "ld" to keep code clean.

    If anything is missing/unknown, fall back to 'default'.
    """
    dt = getattr(fdecl, "dim_tokens", None)
    if dt is None or len(dt) == 0:
        return default

    try:
        dim0 = dt[0].value

        # Try to reduce common "lower:upper" bounds to a simple identifier.
        simple = _try_extract_first_dim_extent_identifier(conv_info, dim0)
        if simple:
            return simple

        # Generic: extent of the first dimension (not the bound itself).
        ld = convert_dim_to_static_size(conv_info, tokens=dim0).strip()
        if not ld:
            return default
        return parenthesize_if_necessary(ld)

    except Exception:
        return default


def _strip_outer_parens_balanced(expr: str) -> str:
    """Strip outer parentheses if they wrap the entire expression."""
    s = expr.strip()
    while s.startswith("(") and s.endswith(")"):
        depth = 0
        balanced_outer = True
        for pos, ch in enumerate(s):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and pos != len(s) - 1:
                    balanced_outer = False
                    break
        if not balanced_outer or depth != 0:
            break
        s = s[1:-1].strip()
    return s


_int_lit_re = re.compile(r"^[0-9]+$")
_ident_re = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _maybe_use_ld_variable(conv_info, ldexpr: str, default_ldname: str) -> str:
    """Prefer a symbolic leading-dimension variable over an inlined expression.

    Rules:
      - If ldexpr is a simple identifier (e.g. "lda", "ldwork"), keep it.
      - Otherwise, use default_ldname (e.g. "ldwork") and remember a declaration
        "INTEGER ldwork = <ldexpr>;" to be emitted near the top of the generated function.

    This avoids hardcoding constants (e.g. "2") and also avoids repeating
    expressions (e.g. "n + nb + 1") in flattened 2D indexing.
    """
    if conv_info is None or default_ldname is None:
        return ldexpr

    core = _strip_outer_parens_balanced(ldexpr)

    # Already symbolic (single token): keep as-is.
    if _ident_re.fullmatch(core):
        return core

    # If the default name already exists as a Fortran identifier, prefer it
    # (we assume the original code assigns it appropriately).
    fproc = getattr(conv_info, "fproc", None)
    fdecl_map = getattr(fproc, "fdecl_by_identifier",
                        None) if fproc is not None else None
    if fdecl_map is not None:
        if default_ldname.lower() in fdecl_map or default_ldname in fdecl_map:
            return default_ldname

    # Otherwise, create a local ld variable initialized with the expression.
    if getattr(conv_info, "ld_constant_decls", None) is None:
        conv_info.ld_constant_decls = {}

    # Keep the first initializer we saw (should be stable for a given array).
    conv_info.ld_constant_decls.setdefault(default_ldname, core)
    return default_ldname


def _emit_constant_ld_decls(top_scope, conv_info) -> None:
    """Emit leading-dimension declarations recorded during conversion."""
    if conv_info is None or top_scope is None:
        return
    ld_map = getattr(conv_info, "ld_constant_decls", None)
    if not ld_map:
        return
    # Ensure we have a valid insertion point.
    if getattr(top_scope, "insert_point", None) is None:
        top_scope.remember_insert_point()

    for ldname, val in ld_map.items():
        # Skip if already declared.
        already = False
        for obj in getattr(top_scope, "data", []):
            if not isinstance(obj, str):
                continue
            if re.match(rf"\s*(?:const\s+)?INTEGER\s+(?:const\s+)?{re.escape(ldname)}\b", obj):
                already = True
                break
        if already:
            continue
        top_scope.top_append(f"INTEGER {ldname} = {val};")


def rewrite_intrinsics(text: str) -> str:
    """Rewrite fem:: intrinsics to C++-friendly forms."""
    def _rewrite_complex_ctor_literals(s: str) -> str:
        """Normalize COMPLEX(real, imag) constructor arguments.

        Policy:
          - If an argument is a pure integer literal (e.g. 0, -1, +2),
            promote it to a real literal (0.0, -1.0, +2.0).
          - Non-literal expressions are left unchanged.

        This fixes patterns such as:
          COMPLEX(lwkopt, 0) -> COMPLEX(lwkopt, 0.0)
        """
        pat = re.compile(r"\bCOMPLEX\s*\(")
        out = []
        i = 0
        n = len(s)
        while True:
            m = pat.search(s, i)
            if not m:
                out.append(s[i:])
                break

            out.append(s[i:m.start()])
            paren_open = m.end() - 1  # points to '('

            # Find matching ')'
            depth = 0
            p = paren_open
            while p < n:
                ch = s[p]
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        break
                p += 1
            if depth != 0:
                # Unbalanced: give up and append rest verbatim
                out.append(s[m.start():])
                break

            inner = s[paren_open + 1:p]
            args = _split_actuals(inner)
            if len(args) == 2:
                new_args = []
                for a in args:
                    a0 = a.strip()
                    if re.fullmatch(r"[+-]?[0-9]+", a0):
                        if a0 in ("0", "+0", "-0"):
                            a0 = "0.0"
                        else:
                            a0 = a0 + ".0"
                    new_args.append(a0)
                out.append("COMPLEX(" + ", ".join(new_args) + ")")
            else:
                # Not a 2-arg ctor: keep as-is
                out.append(s[m.start():p + 1])

            i = p + 1
        return "".join(out)

    def _real_dispatch(arg: str) -> str:
        # Use COMPLEX information collected for the current procedure
        return _real_cast_or_component(
            arg,
            complex_identifiers,
            complex_pointer_identifiers,
        )

    def _int_dispatch(arg: str) -> str:
        """Replacement for INT/FINT.

        Policy:
          - If the expression involves COMPLEX variables, use the same logic
            as DBLE/REAL to get a REAL-valued expression, then castINTEGER.
            (e.g. INT(work(1)) -> castINTEGER(work[0].real()))
          - If the expression does not involve COMPLEX variables, assume it is
            already REAL/INTEGER and use castINTEGER(expr) directly, without
            going through castREAL.
        """
        arg_stripped = arg.strip()

        # Detect whether the expression mentions any COMPLEX identifiers.
        owner = None
        for name in complex_identifiers:
            # Match as an identifier, not as a substring of another name.
            if re.search(r"\b" + re.escape(name) + r"\b", arg_stripped):
                owner = name
                break

        if owner is None:
            # No COMPLEX variables involved: treat as a real/integer expression.
            return f"castINTEGER({arg_stripped})"

        # COMPLEX variables involved: map to a REAL-valued expression first,
        # then cast to INTEGER.
        real_expr = _real_cast_or_component(
            arg,
            complex_identifiers,
            complex_pointer_identifiers,
        )
        return f"castINTEGER({real_expr})"

    # --------------------------------------------------------------
    # 1) Unary intrinsics that need smart handling of parentheses
    # --------------------------------------------------------------
    unary_intrinsics = [
        ("fem::dble",   _real_dispatch),
        ("fem::real",   _real_dispatch),
        ("fem::fint",   _int_dispatch),
        # Some paths emit fem::int instead of fem::fint.
        # Keep behavior consistent: always castINTEGER(...).
        ("fem::int",    _int_dispatch),
        ("fem::INT",    _int_dispatch),
        ("fem::aimag",  _imag_repl),
        ("fem::imag",   _imag_repl),
        ("fem::dimag",  _imag_repl),
        ("fem::conjg",  _conj_repl),
        ("fem::conj",   _conj_repl),
        ("fem::dconjg", _conj_repl),
        ("fem::huge",   _huge_repl),
    ]

    # Apply unary intrinsics first (they keep parentheses balance themselves)
    for func_name, repl in unary_intrinsics:
        if func_name in text:
            text = _rewrite_unary_intrinsic(text, func_name, repl)

    # --------------------------------------------------------------
    # 2) Simple name substitutions
    #    fem::foo(...) -> foo(...) or a common implementation
    # --------------------------------------------------------------
    simple_map = {
        "fem::pow2":  "pow2",
        "fem::mod":   "mod",

        # elementary math
        "fem::cos":   "cos",
        "fem::sin":   "sin",
        "fem::sqrt":  "sqrt",
        "fem::atan2": "atan2",
        "fem::dsqrt": "sqrt",

        # absolute value
        "fem::dabs":  "abs",
        "fem::abs":   "abs",
        "fem::cdabs": "abs",

        # complex constructor
        "fem::dcmplx": "COMPLEX",

        # sign / dsign -> common sign(a, b) helper
        # (C++ side must provide sign(a, b) implementation)
        "fem::sign":  "sign",
        "fem::dsign": "sign",

        # array intrinsics
        "fem::maxloc": "Mmaxloc",
        "fem::maxval": "Mmaxval",
        "fem::minval": "Mminval",
    }

    for src, dst in simple_map.items():
        if src in text:
            text = text.replace(src, dst)

    # --------------------------------------------------------------
    # 3) MAX / MIN need special handling for integer literals, etc.
    # --------------------------------------------------------------
    text = _rewrite_max_min_calls(text)

    # 4) COMPLEX(...) constructor: promote integer literals to real literals.
    text = _rewrite_complex_ctor_literals(text)

    return text


UNARY_BRACKET_PARENS_RE = re.compile(
    r"\[\s*\(\s*(?:0\s*\+\s*)?([A-Za-z_][A-Za-z0-9_]*|[0-9]+)(?:\s*\+\s*0)?\s*\)\s*\]"
)


def rewrite_unary_bracket_parens(text: str) -> str:
    """Remove redundant parentheses for unary subscripts only.

    Examples:
      a[(k)]   -> a[k]
      a[(lda)] -> a[lda]
      a[(0)]   -> a[0]

    Non-unary expressions are preserved:
      a[(i+j)]   (unchanged)
      a[(k-1)]   (unchanged)
    """
    lines = text.splitlines(keepends=True)
    out = []
    for line in lines:
        idx = line.find("//")
        if idx >= 0:
            code, comment = line[:idx], line[idx:]
        else:
            code, comment = line, ""
        code = UNARY_BRACKET_PARENS_RE.sub(r"[\1]", code)
        out.append(code + comment)
    return "".join(out)




def _convert_parentheses_postfix(conv_info, prev_tok, tok, had_str_concat):
    """Return the postfix string for a parentheses token.

    This centralizes the translation of Fortran-style parentheses that follow
    an identifier:
      - array element / slice: a(i), a(i,j), a(i:j), z(i:j, col)
      - function call: foo(...)
      - plain grouping: (expr)

    The same logic must be used by both convert_tokens() and convert_io_loop()
    so that array flattening (lower-bound shifts and leading-dimension
    multiplication) is defined in exactly one place.
    """

    # Detect array reference: a(i) or a(i,j)
    is_array_ref = False
    fdecl = None

    if (prev_tok is not None
            and getattr(prev_tok, 'is_identifier', lambda: False)()
            and conv_info is not None
            and getattr(conv_info, 'fproc', None) is not None):
        try:
            fdecl = conv_info.fproc.get_fdecl(id_tok=prev_tok)
        except Exception:
            fdecl = None

    if (fdecl is not None
            and getattr(fdecl, 'dim_tokens', None) is not None
            and not fdecl.is_user_defined_callable()):
        # Non-callable with dimensions -> array reference
        is_array_ref = True

    if is_array_ref:
        # Check if this is an array slice (contains ':' operator)
        def _contains_colon_op(tokens):
            """Return True if token list contains a colon operator (array slice)."""
            for t in tokens:
                if hasattr(t, 'is_op_with') and t.is_op_with(value=':'):
                    return True
                if hasattr(t, 'value') and isinstance(t.value, (list, tuple)):
                    if _contains_colon_op(t.value):
                        return True
            return False

        is_array_slice = _contains_colon_op(tok.value)

        # Convert indices inside parentheses to C++ array indexing
        idx_str = convert_tokens(
            conv_info=conv_info,
            tokens=tok.value,
            commas=True,
            had_str_concat=had_str_concat)

        # Split by top-level commas, ignoring commas inside nested parentheses
        parts = _split_actuals(idx_str)

        if is_array_slice and len(parts) == 2:
            # 1D slice: a(start:end) -> [__SLICE__(start, end)]
            start_expr = parts[0].strip()
            end_expr = parts[1].strip()
            return f'[__SLICE__({start_expr}, {end_expr})]'

        if is_array_slice and len(parts) == 3:
            # Two possible meanings:
            #   (A) 1D slice with stride: a(start:end:step)
            #   (B) 2D slice on first dimension: z(start:end, col)
            # Disambiguate using declared rank.
            rank = 0
            try:
                rank = len(getattr(fdecl, 'dim_tokens', None) or [])
            except Exception:
                rank = 0

            if rank == 1:
                start_expr = parts[0].strip()
                end_expr = parts[1].strip()
                step_expr = parts[2].strip()
                return f'[__SLICE__({start_expr}, {end_expr}, {step_expr})]'

            # 2D slice on first dimension: [__SLICE2D__(start, end, col, ldname)]
            start_expr = parts[0].strip()
            end_expr = parts[1].strip()
            col_expr = parts[2].strip()
            default_ldname = 'ld' + prev_tok.value.lower()
            ldexpr = get_leading_dimension_expr(conv_info, fdecl, default=default_ldname)
            ldexpr = _maybe_use_ld_variable(conv_info, ldexpr, default_ldname)
            return f'[__SLICE2D__({start_expr}, {end_expr}, {col_expr}, {ldexpr})]'

        if len(parts) == 1:
            # One-dimensional array: a(i)
            i_expr = parts[0].strip()
            lb_expr = get_lower_bound(conv_info, fdecl, 0).strip()

            int_re = re.compile(r'^[0-9]+$')
            name_re = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')
            array_elem_re = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*\s*\[[^\[\]]+\]\s*$')

            def canonical_simple_index(s):
                """Return canonical simple index (identifier/int/array/func) or None."""
                s = s.strip()
                # Strip outer layers of parentheses if they enclose the whole expr.
                while s.startswith('(') and s.endswith(')'):
                    depth = 0
                    balanced_outer = True
                    for pos, ch in enumerate(s):
                        if ch == '(':
                            depth += 1
                        elif ch == ')':
                            depth -= 1
                        if depth == 0 and pos != len(s) - 1:
                            balanced_outer = False
                            break
                    if not balanced_outer:
                        break
                    s = s[1:-1].strip()

                if name_re.fullmatch(s) or int_re.fullmatch(s):
                    return s
                if array_elem_re.fullmatch(s):
                    return s

                # Recognize NAME( ... ) with balanced parentheses.
                if not s:
                    return None
                if not (s[0].isalpha() or s[0] == '_'):
                    return None
                pos = 0
                while pos < len(s) and (s[pos].isalnum() or s[pos] == '_'):
                    pos += 1
                name = s[:pos]
                while pos < len(s) and s[pos].isspace():
                    pos += 1
                if pos >= len(s) or s[pos] != '(':
                    return None

                start_paren = pos
                depth = 0
                end_paren = None
                for i in range(start_paren, len(s)):
                    ch = s[i]
                    if ch == '(':
                        depth += 1
                    elif ch == ')':
                        depth -= 1
                        if depth == 0:
                            end_paren = i
                            break
                if end_paren is None:
                    return None

                pos = end_paren + 1
                while pos < len(s) and s[pos].isspace():
                    pos += 1
                if pos != len(s):
                    return None

                inner = s[start_paren + 1:end_paren]
                return f'{name}({inner})'

            # Keep "1 - 1" (do not fold into 0) for auditability.
            canon_i = canonical_simple_index(i_expr)

            if lb_expr == '0':
                index_expr = canon_i if canon_i is not None else f'({i_expr})'
            else:
                lb_simple = bool(name_re.fullmatch(lb_expr) or int_re.fullmatch(lb_expr))
                idx_simple = canon_i is not None
                if idx_simple and lb_simple:
                    index_expr = f'{canon_i} - {lb_expr}'
                elif idx_simple and not lb_simple:
                    index_expr = f'{canon_i} - ({lb_expr})'
                elif not idx_simple and lb_simple:
                    index_expr = f'({i_expr}) - {lb_expr}'
                else:
                    index_expr = f'({i_expr}) - ({lb_expr})'

            return '[' + index_expr + ']'

        if len(parts) == 2:
            # Two-dimensional array: a(i, j)
            i_expr, j_expr = parts
            i_expr = i_expr.strip()
            j_expr = j_expr.strip()

            default_ldname = 'ld' + prev_tok.value.lower()
            ldexpr = get_leading_dimension_expr(conv_info, fdecl, default=default_ldname)
            ldexpr = _maybe_use_ld_variable(conv_info, ldexpr, default_ldname)

            lb1 = get_lower_bound(conv_info, fdecl, 0).strip()
            lb2 = get_lower_bound(conv_info, fdecl, 1).strip()

            simple_name = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')
            simple_int = re.compile(r'^[0-9]+$')

            def canonical_simple_2(s):
                """Return canonical simple index (identifier/int/array element) or None."""
                s = s.strip()
                while s.startswith('(') and s.endswith(')'):
                    depth = 0
                    balanced_outer = True
                    for pos, ch in enumerate(s):
                        if ch == '(':
                            depth += 1
                        elif ch == ')':
                            depth -= 1
                        if depth == 0 and pos != len(s) - 1:
                            balanced_outer = False
                            break
                    if not balanced_outer:
                        break
                    s = s[1:-1].strip()

                if simple_name.fullmatch(s) or simple_int.fullmatch(s):
                    return s

                # Recognize NAME[ ... ] with balanced brackets.
                if not s:
                    return None
                if not (s[0].isalpha() or s[0] == '_'):
                    return None
                pos = 0
                while pos < len(s) and (s[pos].isalnum() or s[pos] == '_'):
                    pos += 1
                name = s[:pos]
                while pos < len(s) and s[pos].isspace():
                    pos += 1
                if pos >= len(s) or s[pos] != '[':
                    return None

                start_bracket = pos
                depth = 0
                end_bracket = None
                for i in range(start_bracket, len(s)):
                    ch = s[i]
                    if ch == '[':
                        depth += 1
                    elif ch == ']':
                        depth -= 1
                        if depth == 0:
                            end_bracket = i
                            break
                if end_bracket is None:
                    return None

                pos = end_bracket + 1
                while pos < len(s) and s[pos].isspace():
                    pos += 1
                if pos != len(s):
                    return None

                inner = s[start_bracket + 1:end_bracket]
                return f'{name}[{inner}]'

            if lb1 == '1' and lb2 == '1':
                core_i = canonical_simple_2(i_expr)
                if core_i is not None:
                    i_term = f'{core_i} - 1'
                else:
                    if simple_name.match(i_expr) or simple_int.match(i_expr):
                        i_term = f'{i_expr} - 1'
                    else:
                        i_term = f'({i_expr}) - 1'

                core_j = canonical_simple_2(j_expr)
                if core_j is not None:
                    j_term = f'({core_j} - 1)'
                else:
                    if simple_name.match(j_expr) or simple_int.match(j_expr):
                        j_term = f'({j_expr} - 1)'
                    else:
                        j_term = f'(({j_expr}) - 1)'

                index_expr = f'({i_term}) + {j_term}*{ldexpr}'
                return '[' + index_expr + ']'

            # General lower-bound aware path
            def make_offset(idx_expr: str, lb_expr: str) -> str:
                idx_expr = idx_expr.strip()
                lb_expr = lb_expr.strip()
                core = canonical_simple_2(idx_expr)

                if lb_expr == '0':
                    return core if core is not None else f'({idx_expr})'

                lb_simple = bool(simple_name.fullmatch(lb_expr) or simple_int.fullmatch(lb_expr))

                if core is not None:
                    return f'{core} - {lb_expr}' if lb_simple else f'{core} - ({lb_expr})'

                return f'({idx_expr}) - {lb_expr}' if lb_simple else f'({idx_expr}) - ({lb_expr})'

            row_off = make_offset(i_expr, lb1)
            col_off = make_offset(j_expr, lb2)
            index_expr = f'{row_off} + {col_off}*{ldexpr}'
            return '[' + index_expr + ']'

        # Fallback: treat like a normal call/parenthesized expression
        if (cmn_needs_to_be_inserted(conv_info=conv_info, prev_tok=prev_tok)):
            if (len(tok.value) != 0) and (len(tok.value[0].value) != 0):
                op = '(cmn, '
            else:
                op = '(cmn'
        else:
            op = '('
        inner = convert_tokens(
            conv_info=conv_info,
            tokens=tok.value,
            commas=True,
            had_str_concat=had_str_concat)
        return op + inner + ')'

    # Normal function call or parenthesized expression
    inner = convert_tokens(
        conv_info=conv_info,
        tokens=tok.value,
        commas=True,
        had_str_concat=had_str_concat)

    # If the previous token is an identifier and corresponds to a known signature,
    # adjust actual arguments according to pointer/value information.
    if (prev_tok is not None and getattr(prev_tok, 'is_identifier', lambda: False)()):
        sig = _lookup_routine_signature(prev_tok.value)
        if sig is not None:
            force_elems = False
            if FABLE_SMALL_CHAR_VIEW:
                base_name = str(prev_tok.value).split("::")[-1].strip()
                mapped_name = convert_function_name_to_mplapack(base_name)
                force_elems = (str(mapped_name).lower() in _MPLAPACK_CORE_CPP_NAMES)
            inner = _adjust_actuals_using_signature(inner, sig, conv_info, force_elems_call=force_elems)

    # Decide whether we need to inject "cmn" as the first argument.
    if (cmn_needs_to_be_inserted(conv_info=conv_info, prev_tok=prev_tok)):
        if (len(tok.value) != 0) and (len(tok.value[0].value) != 0):
            op = '(cmn, '
        else:
            op = '(cmn'
    else:
        op = '('

    return op + inner + ')'

def convert_tokens(conv_info, tokens, commas=False, had_str_concat=None):
    result = []
    rapp = result.append
    prev_tok = None
    if (had_str_concat is None):
        had_str_concat = mutable(value=False)
    from fable.tokenization import group_power
    for tok in group_power(tokens=tokens):
        if (tok.is_seq()):
            if (len(tok.value) == 2
                    and tok.value[0].is_op_with(value="*")
                    and tok.value[1].is_integer()):
                rapp("star /* %s UNHANDLED */" % tok.value[1].value)
            else:
                rapp(convert_tokens(
                    conv_info=conv_info,
                    tokens=tok.value,
                    commas=False,
                    had_str_concat=had_str_concat))
        elif (tok.is_parentheses()):
            rapp(_convert_parentheses_postfix(
                conv_info=conv_info,
                prev_tok=prev_tok,
                tok=tok,
                had_str_concat=had_str_concat))
        elif (tok.is_implied_do()):
            raise AssertionError
        elif (tok.is_power()):
            rapp(convert_power(conv_info=conv_info, tokens=tok.value))
        else:
            rapp(convert_token(
                vmap=conv_info.vmap,
                leading=(len(result) == 0),
                tok=tok,
                had_str_concat=had_str_concat))
        prev_tok = tok

    # Build base string
    if (commas):
        s = ", ".join(result)
    else:
        s = "".join(result)

    # Rewrite intrinsic calls such as fem::dble(), fem::conjg(), etc.
    s = rewrite_intrinsics(s)
    # Clean redundant parentheses for unary subscripts only, e.g. a[(k)] -> a[k]
    s = rewrite_unary_bracket_parens(s)
    # Rewrite unit-length substrings on plain char* dummy arguments: path(1,1) -> path[(1)-1]
    s = _rewrite_plain_char_ptr_unit_substrings(s, conv_info)
    return s


def convert_to_int_literal(tokens):
    assert tokens is not None and len(tokens) != 0
    if (len(tokens) != 1 or not tokens[0].is_integer()):
        tokens[0].raise_not_supported()
    return int(strip_leading_zeros(string=tokens[0].value))


def convert_data_type(conv_info, fdecl, crhs):
    if (fdecl.data_type is None):
        assert conv_info.fproc is not None
        fdecl.data_type = conv_info.fproc.implicit.get(fdecl.id_tok.value[0])
        if (fdecl.data_type is None):
            raise fdecl.id_tok.raise_semantic_error(msg="Missing data type")
    if (isinstance(fdecl.data_type, str)):
        data_type_code = fdecl.data_type
    else:
        data_type_code = fdecl.data_type.value
    size_tokens = fdecl.size_tokens
    dim_tokens = fdecl.dim_tokens

    if data_type_code == "character":
        if size_tokens is None:
            csize = "1"
        else:
            csize = convert_tokens(conv_info=conv_info, tokens=size_tokens)

        # CHARACTER*1 scalar → plain char (no implicit initialization)
        if csize.strip() == "1" and dim_tokens is None and not FABLE_SMALL_CHAR_VIEW:
            ctype = "char"
            # Do not set crhs here if there is no explicit initializer.
            # If the Fortran code had "CHARACTER*1 normin /'N'/", that will
            # already be reflected in crhs before this point.
        else:
            # For CHARACTER*n (n > 1) scalars we use fem::str<n>.
            # IMPORTANT: do not inject an implicit initializer here.
            #
            # Reason:
            #   If we set crhs="0" at this stage, later declaration emission
            #   cannot see and apply a DATA initializer like:
            #       DATA intstr / '0123456789' /
            #
            # The actual initializer is resolved in convert_declaration() via
            # build_scalar_data_initializers(), falling back to a zero shortcut
            # only when no DATA initializer exists.
            ctype = f"fem::str<{csize}>"

    else:
        def convert_to_ctype_with_size(ctype):
            if (size_tokens is None):
                return ctype
            return "fem::%s_star_%d" % (data_type_code, convert_to_int_literal(
                tokens=size_tokens))
        if (data_type_code == "logical"):
            ctype = convert_to_ctype_with_size(ctype="bool")
        elif (data_type_code == "integer"):
            ctype = convert_to_ctype_with_size(ctype="int")
        elif (data_type_code == "real"):
            ctype = convert_to_ctype_with_size(ctype="float")
        elif (data_type_code == "doubleprecision"):
            if (size_tokens is not None):
                size_tokens[0].raise_syntax_error()
            ctype = "double"
        elif (data_type_code == "byte"):
            if (size_tokens is not None):
                size_tokens[0].raise_syntax_error()
            ctype = "char"
        elif (data_type_code == "complex"):
            if (size_tokens is None):
                ctype = "std::complex<float>"
                if (crhs is None):
                    crhs = "0.0"
            else:
                sz = convert_to_int_literal(tokens=size_tokens)
                if (sz == 8):
                    ctype = "std::complex<float>"
                    if (crhs is None):
                        crhs = "0.0"
                elif (sz == 16):
                    ctype = "std::complex<double>"
                    if (crhs is None):
                        crhs = "0.0"
                elif (sz == 32):
                    ctype = "std::complex<long double>"
                    if (crhs is None):
                        crhs = "0.0"
                else:
                    size_tokens[0].raise_not_supported()
        elif (data_type_code == "doublecomplex"):
            if (size_tokens is not None):
                size_tokens[0].raise_not_supported()
            ctype = "std::complex<double>"
            if (crhs is None):
                crhs = "0.0"
        else:
            raise RuntimeError(
                "Not implemented: data_type_code = %s" % data_type_code)
    return ctype, crhs


def convert_dims(conv_info, dim_tokens):
    need_origin = False
    for tokens in dim_tokens:
        for tok in tokens.value:
            if (tok.is_op_with(value=":")):
                need_origin = True
                break
    dims = []
    for i_dim, tokens in enumerate(dim_tokens):
        cdim = convert_tokens(conv_info=conv_info, tokens=tokens.value)
        if (cdim == "star "):
            cdim = "star"
        else:
            cdim = cdim.replace(",  * ", ", star")  # XXX
        if (need_origin):
            dims.append("dim%d(%s)" % (i_dim+1, cdim))
        else:
            dims.append(cdim)
    if (need_origin):
        result = ".".join(dims)
    else:
        result = "dimension(" + ", ".join(dims) + ")"
    return result


def parenthesize_if_necessary(expr):
    from fable import unsigned_integer_scan, identifier_scan
    if (unsigned_integer_scan(code=expr) == len(expr)
            or identifier_scan(code=expr) == len(expr)):
        return expr
    return "(" + expr + ")"


def convert_dim_to_static_size(conv_info, tokens):
    def conv(toks):
        return convert_tokens(conv_info=conv_info, tokens=toks)
    for i, tok in enumerate(tokens):
        if (tok.is_op_with(value=":")):
            f = parenthesize_if_necessary(expr=conv(toks=tokens[:i]))
            l = parenthesize_if_necessary(expr=conv(toks=tokens[i+1:]))
            return "%s-%s+1" % (l, f)
    return "%s" % conv(toks=tokens)


def convert_dims_to_static_size(conv_info, dim_tokens):
    terms = []
    for tokens in dim_tokens:
        terms.append(convert_dim_to_static_size(
            conv_info, tokens=tokens.value))
    if (len(terms) == 1):
        return terms[0]
    return " * ".join([parenthesize_if_necessary(expr=expr) for expr in terms])


def convert_data_type_and_dims(conv_info, fdecl, crhs, force_arr=False):
    ctype, crhs = convert_data_type(
        conv_info=conv_info, fdecl=fdecl, crhs=crhs)
    dt = fdecl.dim_tokens
    cdims = None
    cfill0 = "fem::fill0"
    if (dt is not None):
        atype = None
        if (not force_arr
                and conv_info.arr_nd_size_max is not None
                and len(dt) <= 3):
            vals = conv_info.fproc.eval_dimensions_simple(
                dim_tokens=dt, allow_power=False)
            if (vals.count(None) == 0):
                sz = math.prod(vals)
                if (sz <= abs(conv_info.arr_nd_size_max)):
                    from fable.read import dimensions_are_simple
                    if (dimensions_are_simple(dim_tokens=dt)):
                        cdims = Auto
                        t = convert_tokens(
                            conv_info=conv_info, tokens=dt, commas=True)
                    else:
                        t = ", ".join(["%d" % v for v in vals])
                    t = "%s, %s" % (t, ctype)
                    if (t.endswith(">")):
                        templs = " "
                    else:
                        templs = ""
                    atype = "arr_%dd<%s%s>" % (len(dt), t, templs)
                    if (conv_info.arr_nd_size_max < 0):
                        cfill0 = "fem::no_fill0"
        if (atype is None):
            if (len(dt) != 1):
                t = "%s, %d" % (ctype, len(dt))
            else:
                t = ctype
            if (t.endswith(">")):
                templs = " "
            else:
                templs = ""
            atype = "arr<%s%s>" % (t, templs)
        ctype = atype
        if (cdims is None):
            cdims = convert_dims(conv_info=conv_info, dim_tokens=dt)
    return ctype, cdims, crhs, cfill0


def ad_hoc_change_arr_to_arr_ref(ctype, cconst=""):
    return ctype.replace("arr<", "arr_%sref<" % cconst, 1)


def zero_shortcut_if_possible(ctype):
    # Convert to simple zero values for basic types
    if ctype == "INTEGER" or ctype == "int":
        return "0"
    elif ctype == "REAL" or ctype == "double":
        return "0.0"
    elif ctype == "float":
        return "0.0"
    elif ctype == "bool" or ctype == "LOGICAL":
        return "false"
    elif ctype == "char":
        return "0"
    elif ctype.startswith("std::complex<"):
        # Complex types in fem:: world; keep old behavior.
        return "0.0"
    elif ctype == "COMPLEX":
        # MPLAPACK scalar COMPLEX: explicit complex zero
        return "COMPLEX(0.0, 0.0)"

    # Original behavior for fem:: types
    if (ctype.startswith("fem::")):
        if (ctype.endswith(">")):
            s = " "
        else:
            s = ""
        return "fem::zero<%s%s>()" % (ctype, s)
    return "fem::%s0" % ctype


def convert_declaration(rapp, conv_info, fdecl, crhs, const):
    ctype, cdims, crhs, cfill0 = convert_data_type_and_dims(
        conv_info=conv_info, fdecl=fdecl, crhs=crhs)
    vname = conv_info.vmapped(fdecl=fdecl)

    # Track COMPLEX locals / COMMON / SAVE by their C++ name.
    dt_code = None
    if fdecl.data_type is not None:
        if isinstance(fdecl.data_type, str):
            dt_code = fdecl.data_type
        else:
            dt_code = getattr(fdecl.data_type, "value", None)
        if dt_code is not None:
            dt_code = dt_code.lower()
    if dt_code in ("complex", "doublecomplex"):
        complex_identifiers.add(vname)

    if (cdims is None):
        # Scalar declaration.

        # Small CHARACTER*n scalars (n > 1) -> plain char arrays.
        # Example:
        #   CHARACTER*2 NAME
        # becomes
        #   char name[2];
        #
        # We only apply this to short fixed-length scalars (length <= _FABLE_SMALL_CHAR_MAX_LEN).
        if dt_code == "character" and fdecl.size_tokens is not None:
            size_expr = convert_tokens(
                conv_info=conv_info,
                tokens=fdecl.size_tokens,
                commas=False,
            ).strip()
            length_is_small = False
            try:
                length_val = int(size_expr)
                if 1 < length_val <= _FABLE_SMALL_CHAR_MAX_LEN:
                    length_is_small = True
            except ValueError:
                # Non-numeric length expression -> fall back to fem::str<...>
                length_is_small = False

            if length_is_small and FABLE_SMALL_CHAR_ENABLED:
                def const_qualifier():
                    if const:
                        return "const "
                    return ""
                # For CHARACTER we map to plain 'char' on the C++ side.
                mplapack_ctype = convert_to_mplapack_type("char")
                rapp("%s%s %s[%s];" % (
                    const_qualifier(), mplapack_ctype, vname, size_expr))
                # Remember that this scalar CHARACTER*n was mapped to a small char[].
                small_char_identifiers.add(vname)
                small_char_identifier_lengths[vname] = length_val
                # If this declaration had an initializer (typically from a merged assignment),
                # we cannot initialize a plain C array in the declaration. Ask the caller
                # to emit an assignment statement that will be expanded into element copies.
                if crhs is not None:
                    return True
                return False

        # For plain CHARACTER*1 (mapped to C++ 'char') without an explicit
        # Fortran initializer, do not inject an artificial '= 0;'.
        # This keeps:
        #   CHARACTER          NORMIN
        # as:
        #   char normin;
        # and relies on the subsequent assignments in the Fortran logic.
        if crhs is None and ctype == "char":
            def const_qualifier():
                if const:
                    return "const "
                return ""
            # Use MPLAPACK-style scalar types (INTEGER, REAL, etc.)
            mplapack_ctype = convert_to_mplapack_type(ctype)
            rapp("%s%s %s;" % (const_qualifier(), mplapack_ctype, vname))
            return False

        if crhs is None:
            # Try DATA-based initializer for simple scalar DATA
            if hasattr(conv_info, "data_initializers"):
                if conv_info.data_initializers is None:
                    conv_info.data_initializers = build_scalar_data_initializers(
                        conv_info)
                init = conv_info.data_initializers.get(
                    fdecl.id_tok.value.lower())
                if init is not None:
                    crhs = init

            # For fixed-length CHARACTER scalars emitted as fem::str<N>, avoid injecting
            # an artificial initializer (no fem::zero<> fallback). If there is no DATA
            # initializer, keep the declaration uninitialized to match Fortran semantics.
            #
            # This preserves DATA initialization when it exists, e.g.:
            #   DATA threq / 2.0d0 / , intstr / '0123456789' /
            if crhs is None and ctype.startswith("fem::str<") and not const:
                def const_qualifier():
                    if const:
                        return "const "
                    return ""
                mplapack_ctype = convert_to_mplapack_type(ctype)
                rapp("%s%s %s;" % (const_qualifier(), mplapack_ctype, vname))
                return False

            # Fallback: zero-initialize as before
            if crhs is None:
                crhs = zero_shortcut_if_possible(ctype=ctype)

        def const_qualifier():
            if const:
                return "const "
            return ""
        # Use MPLAPACK-style scalar types (INTEGER, REAL, etc.)
        mplapack_ctype = convert_to_mplapack_type(ctype)
        rapp("%s%s %s = %s;" % (const_qualifier(), mplapack_ctype, vname, crhs))
        return False

    # Array declaration (local fixed-size arrays).
    # MPLAPACK reference code prefers plain C arrays for small work arrays,
    # e.g. "INTEGER isave[3];".
    # If this array was already emitted as a hoisted static DATA initializer,
    # suppress the plain declaration here.
    try:
        if (hasattr(conv_info, "hoisted_data_array_names")
                and fdecl.id_tok.value.lower() in conv_info.hoisted_data_array_names):
            return False
    except Exception:
        pass
    def const_qualifier():
        if (const):
            return "const "
        return ""

    # Determine element type (scalar) and map it to MPLAPACK typedefs.
    elem_ctype = convert_data_type(
        conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
    mplapack_elem_ctype = convert_to_mplapack_type(elem_ctype)
    # Prefer plain C arrays when all extents are compile-time constants.
    # If any extent depends on runtime values (e.g. MAX(M,N)), fall back to a
    # heap allocation managed by std::unique_ptr<T[]> and expose a raw pointer
    # with the original variable name for downstream indexing.
    vals = conv_info.fproc.eval_dimensions_simple(
        dim_tokens=fdecl.dim_tokens, allow_power=False)
    is_dynamic = (vals is None or vals.count(None) != 0)

    # Prefer symbolic size expressions for local work arrays.
    # For 1D arrays like RWORK(MAXDIM) or WORK(4*MAXDIM) we want:
    #   REAL rwork[maxdim];
    #   COMPLEX work[4 * maxdim];
    # For higher-rank arrays like WORK13(LDWORK,NBMAX) we flatten but keep
    # the product of extents symbolically:
    #   COMPLEX work13[ldwork * nbmax];
    dim_expr = convert_tokens(
        conv_info=conv_info,
        tokens=fdecl.dim_tokens,
        commas=True,
    ).strip()

    # Build a single "number of elements" expression.
    # Policy: never fold extents into numeric literals, even if they are compile-time constants.
    dt_rank = len(fdecl.dim_tokens) if fdecl.dim_tokens is not None else 0

    if dt_rank == 1:
        size_expr = dim_expr
    else:
        parts = _split_actuals(dim_expr)
        if parts and len(parts) == dt_rank:
            size_expr = " * ".join(parts)
        else:
            # Conservative fallback: compute a symbolic product of extents.
            size_expr = convert_dims_to_static_size(
                conv_info=conv_info, dim_tokens=fdecl.dim_tokens)

    if is_dynamic:
        # Runtime-sized local array: allocate on the heap and expose a raw pointer.
        # Example:
        #   DOUBLE PRECISION WNRM(MAX(M,N))
        # ->
        #   std::unique_ptr<REAL[]> fable_wnrm_storage(new REAL[max(m, n)]);
        #   REAL *wnrm = __wnrm_storage.get();
        storage_name = f"__{vname}_storage"
        rapp("%sstd::unique_ptr<%s[]> %s(new %s[%s]);" % (
            const_qualifier(), mplapack_elem_ctype, storage_name,
            mplapack_elem_ctype, size_expr))
        rapp("%s *%s = %s.get();" % (
            mplapack_elem_ctype, vname, storage_name))
        return False

    # Constant-size local array: keep as a plain C array.
    rapp("%s%s %s[%s];" % (
        const_qualifier(), mplapack_elem_ctype, vname, size_expr))
    return False


class scope(object):

    __slots__ = [
        "parent",
        "opening_text",
        "auto_close_parent",
        "closing_text",
        "data",
        "insert_point",
        "trailing_statement_label_index",
        "tail"]

    def __init__(O, parent, opening_text=None, auto_close_parent=False):
        O.parent = parent
        O.opening_text = opening_text
        O.auto_close_parent = auto_close_parent
        O.closing_text = None
        O.data = []
        O.insert_point = None
        O.trailing_statement_label_index = None
        O.tail = None

    def current_point(O):
        return len(O.data)

    def point_is_current(O, point):
        return (point == len(O.data))

    def remember_insert_point(O):
        O.insert_point = len(O.data)

    def insert_point_is_current(O):
        return O.point_is_current(point=O.insert_point)

    def top_append(O, obj):
        assert O.insert_point is not None
        # Inserting into O.data shifts indices; keep trailing label index consistent.
        if (O.trailing_statement_label_index is not None):
            if (O.insert_point <= O.trailing_statement_label_index):
                O.trailing_statement_label_index += 1
            else:
                # Insert after the label: it is no longer trailing.
                O.trailing_statement_label_index = None
        O.data.insert(O.insert_point, obj)
        O.insert_point += 1

    def append(O, obj):
        O.data.append(obj)
        O.trailing_statement_label_index = None

    def append_statement_label(O, label):
        O.trailing_statement_label_index = len(O.data)
        O.data.append("statement_%s:" % label)

    def append_comment(O, line):
        O.data.append(line)

    def open_nested_scope(O, opening_text, auto_close_parent=False):
        O.trailing_statement_label_index = None
        return scope(
            parent=O, opening_text=opening_text, auto_close_parent=auto_close_parent)

    def finalize(O):
        if (O.trailing_statement_label_index is not None):
            O.data[O.trailing_statement_label_index] += ";"

    def close_nested_scope(O):
        assert O.opening_text is not None
        assert O.closing_text is None
        assert O.tail is None
        O.finalize()
        O.closing_text = ["}"]
        head = O
        while (head.parent.tail is head):
            head = head.parent
        head.parent.data.append(head)
        if (head.auto_close_parent):
            return head.parent.close_nested_scope()
        return head.parent

    def attach_tail(O, opening_text):
        assert O.opening_text is not None
        assert O.closing_text is None
        assert O.tail is None
        O.finalize()
        O.closing_text = ["}"]
        O.tail = scope(parent=O, opening_text=opening_text)
        return O.tail

    def collect_text(O, callback, indent="  "):
        for obj in O.data:
            if (isinstance(obj, scope)):
                curr = obj
                while (curr is not None):
                    for text in curr.opening_text:
                        callback(indent+text)
                    curr.collect_text(callback=callback, indent=indent+"  ")
                    for text in curr.closing_text:
                        callback(indent+text)
                    curr = curr.tail
            else:
                callback(indent+obj)


def convert_io_statement_with_err(
        conv_info,
        curr_scope,
        io_function,
        io_function_specialization,
        io_call_args,
        iolist):
    if (iolist.err is None):
        io_scope = curr_scope
    else:
        io_scope = curr_scope.open_nested_scope(opening_text=["try {"])
    io_scope.append("cmn.io.%s%s(%s)" % (
        io_function, io_function_specialization, io_call_args))
    for slot in iolist.chain:
        tokens = getattr(iolist, slot)
        if (tokens is not None):
            carg = convert_tokens(conv_info=conv_info, tokens=tokens)
            io_scope.append("  .%s(%s)" % (slot, carg))
    io_scope.data[-1] += ";"
    if (io_scope is not curr_scope):
        io_scope.close_nested_scope()
        catch_scope = curr_scope.open_nested_scope(
            opening_text=["catch (fem::io_err const&) {"])
        clabel = convert_tokens(conv_info=conv_info, tokens=iolist.err)
        catch_scope.append("goto statement_%s;" % clabel)
        catch_scope.close_nested_scope()


def is_simple_do_last(tokens):
    i = 0
    if (len(tokens) == 2 and tokens[0].is_unary_plus_or_minus()):
        i = 1
    if (i+1 == len(tokens)):
        tok = tokens[i]
        return tok.is_identifier() or tok.is_integer()
    return False


# this part is very dangerous. It can change the meanings of the loops.
def convert_to_fem_do(conv_info, parent_scope, i_tok, fls_tokens):
    assert 2 <= len(fls_tokens) <= 3
    i = convert_token(vmap=conv_info.vmap, leading=True, tok=i_tok)
    f = convert_tokens(conv_info=conv_info, tokens=fls_tokens[0].value)
    l = convert_tokens(conv_info=conv_info, tokens=fls_tokens[1].value)
    if (len(fls_tokens) == 3):
        s = convert_tokens(conv_info=conv_info, tokens=fls_tokens[2].value)
        if (s.lstrip('+-').isdigit()):
            if int(s) >= 0:
                return parent_scope.open_nested_scope(
                    opening_text=["for(%s=%s; %s<=%s; %s=%s+%s) {" % (i, f, i, l, i, i, s)])
            else:
                return parent_scope.open_nested_scope(
                    opening_text=["for(%s=%s; %s>=%s; %s=%s%s) {" % (i, f, i, l, i, i, s)])
        else:
            if s.lstrip().startswith('-'):
                # Step starts with '-', definitely negative.
                return parent_scope.open_nested_scope(
                    opening_text=["for(%s=%s; %s>=%s; %s=%s%s) {" % (i, f, i, l, i, i, s)])
            else:
                # Step is a complex expression; sign unknown at compile time.
                # Use ternary operator to select the correct loop condition.
                return parent_scope.open_nested_scope(
                    opening_text=["for(%s=%s; %s<=%s; %s=%s+%s) {" % (i, f, i, l, i, i, s)])

    if (conv_info.fem_do_safe):
        return parent_scope.open_nested_scope(
            opening_text=["for(%s=%s; %s<=%s; %s=%s+1) {" % (i, f, i, l, i, i)])
    if (is_simple_do_last(tokens=fls_tokens[1].value)):
        return parent_scope.open_nested_scope(
            opening_text=["FEM_DO(%s, %s, %s) {" % (i, f, l)])
    scope_for_last = parent_scope.open_nested_scope(opening_text=["{"])
    scope_for_last.append("int fem_do_last = %s;" % l)
    return scope_for_last.open_nested_scope(
        opening_text=["FEM_DO(%s, %s, fem_do_last) {" % (i, f)],
        auto_close_parent=True)


def find_implied_dos(result, tokens):
    assert isinstance(tokens, list)
    for i, tok in enumerate(tokens):
        if (tok.is_seq_or_parentheses()):
            find_implied_dos(result=result, tokens=tok.value)
        elif (tok.is_implied_do()):
            result.append(tok)


def convert_io_loop(
        io_scope, io_op, conv_info, tokens, cbuf=None, had_str_concat=None):
    class cbuffer(object):
        __slots__ = ["strings", "leading"]

        def __init__(O):
            O.strings = []
            O.leading = True

        def append(O, string):
            O.strings.append(string)
            O.leading = False

        def append_comma(O):
            if (len(O.strings) != 0 and O.strings[-1] != ", "):
                O.strings.append(", ")
            O.leading = True

        def remove_trailing_comma(O):
            if (len(O.strings) != 0):
                assert O.strings[-1] == ", "
                O.strings.pop()
                assert O.leading
                O.leading = False

        def append_opening_parenthesis(O):
            O.strings.append("(")
            O.leading = True

        def append_closing_parenthesis(O):
            assert len(O.strings) != 0
            if (O.strings[-1] == ", "):
                O.remove_trailing_comma()
            O.strings.append(")")

        def flush(O):
            O.remove_trailing_comma()
            if (len(O.strings) != 0):
                io_scope.append("%s, %s;" % (io_op, "".join(O.strings)))
                O.strings = []
    if (cbuf is None):
        cbuf = cbuffer()
        owning_cbuf = True
    else:
        owning_cbuf = False
    prev_tok = None
    if (had_str_concat is None):
        had_str_concat = mutable(value=False)
    from fable.tokenization import group_power
    for tok in group_power(tokens=tokens):
        if (tok.is_seq()):
            convert_io_loop(
                io_scope,
                io_op,
                conv_info,
                tokens=tok.value,
                cbuf=cbuf,
                had_str_concat=had_str_concat)
            cbuf.append_comma()
        elif (tok.is_parentheses()):
            cbuf.append(_convert_parentheses_postfix(
                conv_info=conv_info,
                prev_tok=prev_tok,
                tok=tok,
                had_str_concat=had_str_concat))
        elif (tok.is_implied_do()):
            cbuf.flush()
            from fable.tokenization import implied_do_info
            idi = implied_do_info(tokens=tok.value)
            do_scope = convert_to_fem_do(
                conv_info=conv_info,
                parent_scope=io_scope,
                i_tok=idi.id_tok,
                fls_tokens=idi.fls_tokens)
            convert_io_loop(
                io_scope=do_scope,
                io_op=io_op,
                conv_info=conv_info,
                tokens=tok.value[:idi.dlist_size],
                had_str_concat=had_str_concat)
            do_scope.close_nested_scope()
            return
        elif (tok.is_power()):
            cbuf.append(convert_power(conv_info=conv_info, tokens=tok.value))
        else:
            cbuf.append(convert_token(
                vmap=conv_info.vmap,
                leading=cbuf.leading,
                tok=tok,
                had_str_concat=had_str_concat))
        prev_tok = tok
    if (owning_cbuf):
        cbuf.flush()


def equivalence_align_with_arg(conv_info, top_scope, identifier, tok_seq):
    assert tok_seq.is_seq()
    tokens = tok_seq.value
    assert len(tokens) > 0
    if (len(tokens) == 1):
        return ""
    cindices = []
    for i in range(1, len(tokens)):
        tok = tokens[i]
        if (i == 3 or not tok.is_parentheses()):
            tok.raise_semantic_error()
        declare_identifiers_parameter_recursion(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=top_scope,
            tokens=tokens[i].value)
        cindices.append(convert_tokens(
            conv_info=conv_info, tokens=tokens[i].value, commas=True))
    fdecl = conv_info.fproc.fdecl_by_identifier[identifier]
    if (len(cindices) == 1):
        if (fdecl.dim_tokens is not None):
            return "arr_index(%s)" % cindices[0]
        if (fdecl.data_type.value != "character"):
            tok_seq.raise_semantic_error()
        return "str_index(%s)" % cindices[0]
    if (len(cindices) != 2
        or fdecl.data_type.value != "character"
            or fdecl.dim_tokens is None):
        tok_seq.raise_semantic_error()
    return "arr_index(%s)(%s)" % tuple(cindices)


def cconst(fdecl, short):
    if (fdecl.is_modified):
        return ""
    if (short):
        return "c"
    return " const"


def convert_to_mbr_bind(
        conv_info,
        top_scope,
        variant_bind_chain,
        mbr_buffer,
        bind_buffer,
        identifier):
    fdecl = conv_info.fproc.fdecl_by_identifier[identifier]
    ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
    if (fdecl.dim_tokens is None):
        cdims_parens = ""
    else:
        declare_identifiers_parameter_recursion(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=top_scope,
            tokens=fdecl.dim_tokens)
        cdims = convert_dims(conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
        cdims_parens = "(" + cdims + ")"
    identifier = fdecl.id_tok.value
    ctype_targ = ctype
    if (ctype_targ.endswith(">")):
        ctype_targ += " "
    mbr_buffer.append("mbr<%s> %s%s;" % (ctype_targ, identifier, cdims_parens))
    conv_info.set_vmap_force_local(fdecl=fdecl)
    vname = conv_info.vmapped(fdecl=fdecl)
    if (fdecl.use_count == 0):
        pr = "/* "
        eq = "*/"
        clm = " */ "
        prm = " /* "
    else:
        pr = ""
        eq = "="
        clm = ""
        prm = ""
    if (fdecl.dim_tokens is None):
        if (fdecl.data_type.value == "character"):
            binding = "%sstr_%sref %s %s %s.bind_str();" % (
                pr, cconst(fdecl=fdecl, short=True), vname, eq, variant_bind_chain)
        else:
            binding = "%s%s%s& %s %s %s.bind<%s>();" % (
                pr, ctype, cconst(fdecl=fdecl, short=False), vname,
                eq, variant_bind_chain, ctype_targ)
    else:
        if (fdecl.data_type.value == "character"):
            ref_dim = "%d" % len(fdecl.dim_tokens)
            if (ref_dim == "1"):
                ref_dim = ""
            binding = "%sstr_arr_%sref<%s> %s(%s%s.bind_str()%s, %s)%s;" % (
                pr, cconst(fdecl=fdecl, short=True), ref_dim, vname,
                clm, variant_bind_chain, prm, cdims, clm)
        else:
            ref_dim = ", %d" % len(fdecl.dim_tokens)
            if (ref_dim == ", 1"):
                ref_dim = ""
            binding = "%sarr_%sref<%s%s> %s(%s%s.bind<%s>()%s, %s)%s;" % (
                pr, cconst(fdecl=fdecl, short=True), ctype, ref_dim, vname,
                clm, variant_bind_chain, ctype, prm, cdims, clm)
    bind_buffer.append(binding)


def assemble_allocate_line_lists(
        conv_info,
        top_scope,
        variant_bind_chain,
        mbr_buffer,
        bind_buffer,
        allocate_line_lists,
        equiv_tok_cluster,
        identifier):
    if (allocate_line_lists[-1] == [" "]):
        allocate_line_lists.pop()
    i_mbr_by_identifer = {identifier: 0}
    eq_identifiers = [identifier]
    i_block = len(allocate_line_lists)
    allocate_line_lists.append(None)
    for equiv_tok in equiv_tok_cluster:
        align_with = ".align"
        for tok_seq in equiv_tok.value:
            eq_identifier = tok_seq.value[0].value
            i_mbr = i_mbr_by_identifer.get(eq_identifier)
            if (i_mbr is None):
                i_mbr_by_identifer[eq_identifier] = i_mbr = len(eq_identifiers)
                eq_identifiers.append(eq_identifier)
                if (bind_buffer is not None):
                    convert_to_mbr_bind(
                        conv_info=conv_info,
                        top_scope=top_scope,
                        variant_bind_chain=variant_bind_chain,
                        mbr_buffer=mbr_buffer,
                        bind_buffer=bind_buffer,
                        identifier=eq_identifier)
            allocate_line_lists.append([
                "    %s<%d>(%s)" % (
                    align_with,
                    i_mbr+1,
                    equivalence_align_with_arg(
                        conv_info=conv_info,
                        top_scope=top_scope,
                        identifier=eq_identifier,
                        tok_seq=tok_seq))])
            align_with = " .with"
    allocate_line_lists[i_block] = \
        ["  equivalence(%s)" % (", ".join(eq_identifiers))]
    allocate_line_lists[-1][-1] += ","
    allocate_line_lists.append([" "])


def add_allocate_lines_to_mbr_scope(allocate_line_lists, mbr_buffer):
    if (allocate_line_lists[-1] == [" "]):
        allocate_line_lists.pop()
    if (allocate_line_lists[-1][-1][-1] == ","):
        allocate_line_lists[-1][-1] = allocate_line_lists[-1][-1][:-1]
    if (len(allocate_line_lists) == 1):
        allocate_line_lists[-1][-1] += ";"
    else:
        allocate_line_lists.append([";"])
    for line_list in allocate_line_lists:
        mbr_buffer.append(" ".join(line_list))


def convert_variant_allocate_and_bindings(conv_info, top_scope):
    result_buffers = group_args(
        first_time=[],
        loc_equivalences=[],
        bindings=[])
    equiv_info = conv_info.fproc.equivalence_info()
    equiv_tok_clusters = equiv_info.equiv_tok_clusters
    for common_name, common_fdecl_list in conv_info.fproc.common.items():
        vcn = conv_info.fproc.conv_hook.variant_common_names
        if (vcn is None or common_name not in vcn):
            continue
        top_scope.append(
            "common_variant %s(cmn.common_%s, sve.%s_bindings);" % (
                (common_name,)*3))
        mbr_buffer = []
        result_buffers.first_time.append(mbr_buffer)
        allocate_line_lists = [["%s.allocate()," % common_name]]
        for fdecl in common_fdecl_list:
            identifier = fdecl.id_tok.value
            convert_to_mbr_bind(
                conv_info=conv_info,
                top_scope=top_scope,
                variant_bind_chain=common_name,
                mbr_buffer=mbr_buffer,
                bind_buffer=result_buffers.bindings,
                identifier=identifier)
            equiv_tok_cluster = equiv_info.equiv_tok_cluster_by_identifier.get(
                identifier)
            if (equiv_tok_cluster is None):
                allocate_line_lists[-1].append(identifier+",")
            else:
                assemble_allocate_line_lists(
                    conv_info=conv_info,
                    top_scope=top_scope,
                    variant_bind_chain=common_name,
                    mbr_buffer=mbr_buffer,
                    bind_buffer=result_buffers.bindings,
                    allocate_line_lists=allocate_line_lists,
                    equiv_tok_cluster=equiv_tok_cluster,
                    identifier=identifier)
        add_allocate_lines_to_mbr_scope(
            allocate_line_lists=allocate_line_lists, mbr_buffer=mbr_buffer)
    #
    cei = conv_info.fproc.classified_equivalence_info()
    if (len(cei.save.equiv_tok_clusters) != 0):
        top_scope.append(
            "save_equivalences sve_equivalences(sve.save_equivalences);")
        mbr_buffer = []
        result_buffers.first_time.append(mbr_buffer)
        mbr_defined_already = set()
        for equiv_tok_cluster in cei.save.equiv_tok_clusters:
            for equiv_tok in equiv_tok_cluster:
                for tok_seq in equiv_tok.value:
                    identifier = tok_seq.value[0].value
                    if (identifier in mbr_defined_already):
                        continue
                    mbr_defined_already.add(identifier)
                    convert_to_mbr_bind(
                        conv_info=conv_info,
                        top_scope=top_scope,
                        variant_bind_chain="sve_equivalences",
                        mbr_buffer=mbr_buffer,
                        bind_buffer=result_buffers.bindings,
                        identifier=identifier)
        allocate_line_lists = [["sve_equivalences.allocate(),"]]
        for equiv_tok_cluster in cei.save.equiv_tok_clusters:
            assemble_allocate_line_lists(
                conv_info=conv_info,
                top_scope=top_scope,
                variant_bind_chain=None,
                mbr_buffer=mbr_buffer,
                bind_buffer=None,
                allocate_line_lists=allocate_line_lists,
                equiv_tok_cluster=equiv_tok_cluster,
                identifier=equiv_tok_cluster[0].value[0].value[0].value)
        add_allocate_lines_to_mbr_scope(
            allocate_line_lists=allocate_line_lists, mbr_buffer=mbr_buffer)
    #
    if (len(cei.local.equiv_tok_clusters) != 0):
        mbr_defined_already = set()
        for equiv_tok_cluster in cei.local.equiv_tok_clusters:
            for equiv_tok in equiv_tok_cluster:
                for tok_seq in equiv_tok.value:
                    identifier = tok_seq.value[0].value
                    if (identifier in mbr_defined_already):
                        continue
                    mbr_defined_already.add(identifier)
                    convert_to_mbr_bind(
                        conv_info=conv_info,
                        top_scope=top_scope,
                        variant_bind_chain="loc_equivalences",
                        mbr_buffer=result_buffers.loc_equivalences,
                        bind_buffer=result_buffers.bindings,
                        identifier=identifier)
        allocate_line_lists = [["loc_equivalences.allocate(),"]]
        for equiv_tok_cluster in cei.local.equiv_tok_clusters:
            assemble_allocate_line_lists(
                conv_info=conv_info,
                top_scope=top_scope,
                variant_bind_chain=None,
                mbr_buffer=result_buffers.loc_equivalences,
                bind_buffer=None,
                allocate_line_lists=allocate_line_lists,
                equiv_tok_cluster=equiv_tok_cluster,
                identifier=equiv_tok_cluster[0].value[0].value[0].value)
        add_allocate_lines_to_mbr_scope(
            allocate_line_lists=allocate_line_lists,
            mbr_buffer=result_buffers.loc_equivalences)
    #
    return result_buffers


def convert_data(conv_info, data_init_scope):
    for nlist, clist in conv_info.fproc.data:
        ccs = []
        have_repetitions = False
        tok_types = set()
        for repetition_tok, ctoks in clist:
            i = 0
            if (ctoks[0].is_unary_plus_or_minus() and len(ctoks) > 1):
                i = 1
            tok_types.add(ctoks[i].type())
            cc = convert_tokens(conv_info=conv_info, tokens=ctoks)
            if (repetition_tok is not None):
                have_repetitions = True
                cr = convert_tokens(conv_info=conv_info,
                                    tokens=[repetition_tok])
                cc = "%s*datum(%s)" % (cr, cc)
            ccs.append(cc)
        homogeneous_ctype = None
        if (conv_info.data_specializations
                and not have_repetitions
                and len(tok_types) == 1):
            homogeneous_ctype = {
                "integer": "int",
                "hexadecimal": "int",
                "real": "float",
                "double_precision": "double",
                "logical": "bool",
                "string": "char*",
                "complex": None  # TODO
            }.get(list(tok_types)[0])

        def data_values_blocked():
            data_scope.append("fem::data_values data;")
            for i_block in range(0, len(ccs), conv_info.data_values_block_size):
                data_scope.append(
                    "data.values, %s;"
                    % ", ".join(ccs[i_block:i_block+conv_info.data_values_block_size]))

        def values_for_data_of_type():
            data_scope.append(
                "static const %s values[] = {" % homogeneous_ctype)
            data_scope.append("  %s" % ", ".join(ccs))
            data_scope.append("};")
            if (homogeneous_ctype == "char*"):
                s = "_str"
            else:
                s = "<%s>" % homogeneous_ctype
            return s

        def have_no_array_targets():
            for tok_seq in nlist:
                if (len(tok_seq.value) == 1):
                    fdecl = conv_info.fproc.get_fdecl(id_tok=tok_seq.value[0])
                    if (fdecl.dim_tokens is not None):
                        return False
            return True
        implied_dos = []
        find_implied_dos(result=implied_dos, tokens=nlist)
        if (len(implied_dos) == 0):
            if (conv_info.data_specializations
                    and len(nlist) == len(ccs)
                    and have_no_array_targets()):
                for tok_seq, cc in zip(nlist, ccs):
                    cn = convert_tokens(conv_info=conv_info,
                                        tokens=tok_seq.value)
                    data_init_scope.append("%s = %s;" % (cn, cc))
            else:
                cn = convert_tokens(conv_info=conv_info,
                                    tokens=nlist, commas=True)
                if (homogeneous_ctype is None):
                    if (len(ccs) <= conv_info.data_values_block_size):
                        data_init_scope.append(
                            "fem::data((values, %s)), %s;" % (", ".join(ccs), cn))
                    else:
                        data_scope = data_init_scope.open_nested_scope(
                            opening_text=["{"])
                        data_values_blocked()
                        data_scope.append("data, %s;" % cn)
                        data_scope.close_nested_scope()
                elif (len(ccs) != 1):
                    if (len(conv_info.fproc.data) == 1
                            and len(conv_info.fproc.conv_hook.variant_common_names) == 0):
                        data_scope = data_init_scope
                    else:
                        data_scope = data_init_scope.open_nested_scope(
                            opening_text=["{"])
                    s = values_for_data_of_type()
                    data_scope.append(
                        "fem::data_of_type%s(FEM_VALUES_AND_SIZE)," % s)
                    data_scope.append("  %s;" % cn)
                    if (data_scope is not data_init_scope):
                        data_scope.close_nested_scope()
                else:
                    data_init_scope.append("%s = %s;" % (cn, ccs[0]))
        else:
            if (len(conv_info.fproc.data) == 1
                    and len(conv_info.fproc.conv_hook.variant_common_names) == 0
                    and len(ccs) <= conv_info.data_values_block_size):
                data_scope = data_init_scope
            else:
                data_scope = data_init_scope.open_nested_scope(
                    opening_text=["{"])
            if (homogeneous_ctype is None):
                if (len(ccs) <= conv_info.data_values_block_size):
                    data_scope.append(
                        "fem::data_values data((values, %s));" % ", ".join(ccs))
                else:
                    data_values_blocked()
            else:
                s = values_for_data_of_type()
                data_scope.append(
                    "fem::data_of_type%s data(FEM_VALUES_AND_SIZE);" % s)
            convert_io_loop(
                io_scope=data_scope,
                io_op="data",
                conv_info=conv_info,
                tokens=nlist)
            if (data_scope is not data_init_scope):
                data_scope.close_nested_scope()


def build_scalar_data_initializers(conv_info):
    """Collect simple scalar DATA initializers: name -> C++ literal.

    We only handle the simple case:
      DATA A,B,C / v1, v2, v3 /
    with no implied DO loops, no repetitions, and scalar targets.
    """
    init = {}
    for nlist, clist in conv_info.fproc.data:
        # Build ccs exactly as in convert_data
        ccs = []
        have_repetitions = False
        tok_types = set()
        for repetition_tok, ctoks in clist:
            i = 0
            if ctoks and ctoks[0].is_unary_plus_or_minus() and len(ctoks) > 1:
                i = 1
            tok_types.add(ctoks[i].type())
            cc = convert_tokens(conv_info=conv_info, tokens=ctoks)
            if repetition_tok is not None:
                # Repeated DATA: skip in this helper for now.
                have_repetitions = True
            ccs.append(cc)

        if have_repetitions:
            continue

        # Skip implied DO loops
        implied_dos = []
        find_implied_dos(result=implied_dos, tokens=nlist)
        if implied_dos:
            continue

        if len(nlist) != len(ccs):
            continue

        # Collect scalar identifiers only
        scalar_ids = []
        for tok_seq in nlist:
            toks = tok_seq.value
            if len(toks) != 1 or not toks[0].is_identifier():
                scalar_ids = []
                break
            id_tok = toks[0]
            fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
            if fdecl is None or fdecl.dim_tokens is not None:
                # Array target or unknown: give up
                scalar_ids = []
                break
            scalar_ids.append(id_tok.value.lower())

        if not scalar_ids:
            continue

        for name, cc in zip(scalar_ids, ccs):
            init[name] = cc

    return init


def build_array_data_initializers(conv_info):
    """Collect constant DATA initializers for local arrays.

    Supported patterns (safe subset):
      - Single target array name
      - Any rank >= 1 with compile-time constant extents
      - No implied DO loops
      - Optional repetition factors (e.g. 5*4) where the repetition is an integer literal

    Returns:
      dict name_lower -> (elem_ctype:str, dims_ints:list[int], init_list:list[str], rank:int)

    Notes:
      - Multi-dimensional DATA is flattened in Fortran storage order as it appears in the DATA list.
      - We require that the expanded DATA list length equals the product of extents.
    """
    init = {}
    fproc = getattr(conv_info, "fproc", None)
    if fproc is None:
        return init

    for nlist, clist in getattr(fproc, "data", []) or []:
        # Skip implied DO loops
        implied_dos = []
        find_implied_dos(result=implied_dos, tokens=nlist)
        if implied_dos:
            continue

        # Only handle single simple identifier target
        if len(nlist) != 1:
            continue
        toks = getattr(nlist[0], "value", None)
        if not toks or len(toks) != 1 or not toks[0].is_identifier():
            continue
        id_tok = toks[0]

        try:
            fdecl = fproc.get_fdecl(id_tok=id_tok)
        except Exception:
            fdecl = None
        if fdecl is None or getattr(fdecl, "dim_tokens", None) is None:
            continue

        rank = len(getattr(fdecl, "dim_tokens", []) or [])
        if rank <= 0:
            continue

        vals = fproc.eval_dimensions_simple(dim_tokens=fdecl.dim_tokens, allow_power=False)
        if vals is None or vals.count(None) != 0:
            continue
        try:
            dims_ints = [int(v) for v in vals]
        except Exception:
            continue
        if any(v <= 0 for v in dims_ints):
            continue

        try:
            n_elems = int(math.prod(dims_ints))
        except Exception:
            continue
        if n_elems <= 0:
            continue

        expanded = []
        ok = True
        for repetition_tok, ctoks in clist:
            cc = convert_tokens(conv_info=conv_info, tokens=ctoks)
            rep = 1
            if repetition_tok is not None:
                rep_s = convert_tokens(conv_info=conv_info, tokens=[repetition_tok]).strip()
                try:
                    rep = int(rep_s)
                except Exception:
                    ok = False
                    break
                if rep <= 0:
                    ok = False
                    break
            expanded.extend([cc] * rep)
            if len(expanded) > n_elems:
                ok = False
                break

        if not ok:
            continue
        if len(expanded) != n_elems:
            continue

        try:
            elem_ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
        except Exception:
            elem_ctype = None

        init[str(id_tok.value).lower()] = (elem_ctype, dims_ints, expanded, rank)

    return init


def declare_identifiers_parameter_recursion(
        conv_info, top_scope, curr_scope, tokens):
    from fable.tokenization import extract_identifiers
    for id_tok in extract_identifiers(tokens=tokens):
        if (id_tok.value in conv_info.vmap):
            continue
        # Intrinsics used in PARAMETER expressions must NOT trigger declaration lookup.
        # However, Fortran allows local variables/parameters to shadow intrinsic names
        # (e.g. INTEGER NINT, INTEGER LEN, REAL DMIN1).
        # Only treat the token as an intrinsic if there is NO declaration record
        # for it in the current procedure.
        if _is_intrinsic_name(id_tok.value):
            fdecl_map = getattr(conv_info.fproc, "fdecl_by_identifier", None)
            if fdecl_map is None:
                _map_intrinsic_vmap(conv_info, id_tok.value)
                continue
            if (id_tok.value.lower() not in fdecl_map
                    and id_tok.value not in fdecl_map):
                _map_intrinsic_vmap(conv_info, id_tok.value)
                continue
        # Tentatively mark as seen to avoid recursion loops.
        # IMPORTANT: Do not store None here; convert_token() may call .lower() on vmap results.
        conv_info.vmap[id_tok.value] = prepend_identifier_if_necessary(
            id_tok.value)
        try:
            fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
        except KeyError:
            # Unknown identifier in PARAMETER expressions: skip recursion.
            continue

        # Some identifiers inside PARAMETER expressions are not parameters
        # (e.g., external/specification functions). Do not crash.
        try:
            req_tokens = fdecl.required_parameter_assignment_tokens()
        except Exception:
            req_tokens = []
        if req_tokens:
            declare_identifiers_parameter_recursion(
                conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
                tokens=req_tokens)

        declare_identifier(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=curr_scope,
            id_tok=id_tok)


def declare_size_dim_identifiers(conv_info, top_scope, curr_scope, fdecl):
    for sd_tokens in [fdecl.size_tokens, fdecl.dim_tokens]:
        if (sd_tokens is None):
            continue
        from fable.tokenization import extract_identifiers
        sd_id_tokens = extract_identifiers(tokens=sd_tokens)
        for sd_id_tok in sd_id_tokens:
            if (sd_id_tok.value == fdecl.id_tok.value):
                sd_id_tok.raise_semantic_error(msg="Recursion in declaration")
            if (sd_id_tok.value in conv_info.vmap):
                continue
            sd_fdecl = conv_info.fproc.get_fdecl(id_tok=sd_id_tok)
            if (sd_fdecl.parameter_assignment_tokens is None):
                sd_crhs = None
            else:
                declare_identifiers_parameter_recursion(
                    conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
                    tokens=sd_fdecl.parameter_assignment_tokens)
                if (sd_id_tok.value in conv_info.fproc.dynamic_parameters):
                    sd_crhs = "cmn.dynamic_params." + sd_id_tok.value
                else:
                    sd_crhs = convert_tokens(
                        conv_info=conv_info, tokens=sd_fdecl.parameter_assignment_tokens)
            if (not conv_info.set_vmap_from_fdecl(fdecl=sd_fdecl)):
                have_goto = (
                    len(conv_info.fproc.target_statement_labels()) != 0)
                if (have_goto):
                    rapp = top_scope.top_append
                else:
                    rapp = top_scope.append
                convert_declaration(
                    rapp=rapp,
                    conv_info=conv_info,
                    fdecl=sd_fdecl,
                    crhs=sd_crhs,
                    const=True)


def simple_equivalence(
        conv_info,
        top_scope,
        curr_scope,
        target_fdecl,
        equiv_tok_cluster):
    assert len(equiv_tok_cluster) != 0
    for equiv_tok in equiv_tok_cluster:
        target_tok_seq = None
        source_tok_seq = None
        source_fdecl = None
        for tok_seq in equiv_tok.value:
            identifier = tok_seq.value[0].value
            if (identifier == target_fdecl.id_tok.value):
                assert target_tok_seq is None
                target_tok_seq = tok_seq
            else:
                fdecl = conv_info.fproc.get_fdecl(id_tok=tok_seq.value[0])
                if (fdecl is not None and fdecl.is_common()):
                    assert source_tok_seq is None
                    source_tok_seq = tok_seq
                    source_fdecl = fdecl
        if (target_tok_seq is not None
                and source_tok_seq is not None):
            break
    else:
        raise AssertionError
    conv_info.set_vmap_force_local(fdecl=target_fdecl)
    if (conv_info.vmap.get(source_fdecl.id_tok.value) is None):
        declare_identifier(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=curr_scope,
            id_tok=source_fdecl.id_tok)
    crhs = convert_tokens(conv_info=conv_info, tokens=[source_tok_seq])
    se = "// SIMPLE EQUIVALENCE"
    if (target_fdecl.data_type.value == "character"):
        clen = convert_tokens(
            conv_info=conv_info, tokens=target_fdecl.size_tokens)
        if (target_fdecl.dim_tokens is None):
            return "str_%sref %s(%s, %s); %s" % (
                cconst(fdecl=target_fdecl, short=True),
                target_fdecl.id_tok.value, crhs, clen, se)
        cdims = convert_dims(
            conv_info=conv_info, dim_tokens=target_fdecl.dim_tokens)
        return "str_arr_%sref<%d> %s(%s, %s, %s); %s" % (
            cconst(fdecl=target_fdecl, short=True),
            len(target_fdecl.dim_tokens),
            target_fdecl.id_tok.value, crhs, clen, cdims, se)
    ctype, cdims = convert_data_type_and_dims(
        conv_info=conv_info, fdecl=target_fdecl, crhs=None, force_arr=True)[:2]
    if (cdims is None):
        return "%s%s& %s = %s; %s" % (
            ctype, cconst(fdecl=target_fdecl, short=False),
            target_fdecl.id_tok.value, crhs, se)
    return "%s %s(%s, %s); %s" % (
        ad_hoc_change_arr_to_arr_ref(
            ctype=ctype, cconst=cconst(fdecl=target_fdecl, short=True)),
        target_fdecl.id_tok.value, crhs, cdims, se)


def declare_identifier(conv_info, top_scope, curr_scope, id_tok, crhs=None):
    # Intrinsic names can be shadowed by local variables/parameters.
    # Only treat this identifier as an intrinsic if there is NO declaration
    # record for it in the current procedure.
    try:
        fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
    except KeyError:
        fdecl = None

    if fdecl is None and _is_intrinsic_name(id_tok.value):
        _map_intrinsic_vmap(conv_info, id_tok.value)
        return crhs is not None

    if fdecl is None:
        # Unknown identifier: do not crash. Keep it as a plain identifier.
        # (This may still lead to a compile error downstream, but conversion continues.)
        conv_info.vmap[id_tok.value] = prepend_identifier_if_necessary(
            id_tok.value)
        return crhs is not None

    conv_info.set_vmap_from_fdecl(fdecl=fdecl)
    have_goto = (len(conv_info.fproc.target_statement_labels()) != 0)

    # Skip CHARACTER declarations for dummy arguments.
    # Their types are already handled in the generated C++ function signature.
    if _is_dummy_character_arg(conv_info, id_tok.value):
        return crhs is not None

    def get_rapp():
        if (have_goto):
            return top_scope.top_append
        return top_scope.append

    if (not fdecl.is_common()):
        equiv_tok_cluster = (
            conv_info.fproc.equivalence_info()
            .equiv_tok_cluster_by_identifier.get(id_tok.value)
        )
        if (equiv_tok_cluster is not None):
            rapp = get_rapp()
            rapp(simple_equivalence(
                conv_info=conv_info,
                top_scope=top_scope,
                curr_scope=curr_scope,
                target_fdecl=fdecl,
                equiv_tok_cluster=equiv_tok_cluster))
            return crhs is not None
    # check if this procedure wants to ignore COMMON/SAVE
    ignore_cs = (
        conv_info.fproc is not None
        and getattr(conv_info.fproc, "conv_hook", None) is not None
        and conv_info.fproc.conv_hook.ignore_common_and_save
    )

    if (fdecl is not None
        and (fdecl.is_local()
             or fdecl.is_parameter()
             or (fdecl.is_save() and ignore_cs))):
        const = False
        have_crhs = (crhs is not None)
        if (have_goto or curr_scope != top_scope):
            crhs = None
        if (crhs is None):
            if (fdecl.parameter_assignment_tokens is not None):
                declare_identifiers_parameter_recursion(
                    conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
                    tokens=fdecl.parameter_assignment_tokens)
                if (id_tok.value in conv_info.fproc.dynamic_parameters):
                    crhs = "cmn.dynamic_params." + prepend_identifier_if_necessary(
                        id_tok.value)
                else:
                    crhs = convert_tokens(
                        conv_info=conv_info, tokens=fdecl.parameter_assignment_tokens)
                const = True
                # If the PARAMETER initializer involves machine-constant intrinsics
                # (RADIX/MINEXPONENT/MAXEXPONENT) or was dropped during preprocessing,
                # do not silently emit a bogus numeric value. Force an explicit
                # UNHANDLED initializer so the C++ build fails loudly.
                if (fdecl.is_parameter()
                        and (fdecl.id_tok.value.lower() in _FORCE_UNHANDLED_PARAMETER_NAMES
                             or _contains_machine_const_intrinsics(crhs))):
                    crhs = "UNHANDLED"
        elif (fdecl.parameter_assignment_tokens is not None):
            id_tok.raise_semantic_error(
                msg="Assignment to PARAMETER %s" % id_tok.value)
        rapp = get_rapp()
        declare_size_dim_identifiers(
            conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
            fdecl=fdecl)
        result = convert_declaration(
            rapp=rapp,
            conv_info=conv_info,
            fdecl=fdecl,
            crhs=crhs,
            const=const)
        if (have_crhs and (have_goto or curr_scope != top_scope)):
            result = True
        return result
    identifier = id_tok.value

    def get_common_name_if_cast_is_needed():
        if (not fdecl.is_common()):
            return None
        if (conv_info.converted_commons_info is None):
            return None
        common_names = conv_info.converted_commons_info.member_registry.get(
            identifier)
        if (common_names is None):
            return None
        if (len(common_names) < 2):
            return None
        return conv_info.fproc.common_name_by_identifier().get(identifier)

    common_name = get_common_name_if_cast_is_needed()
    if (common_name is not None):
        src_var = "static_cast<common_%s&>(cmn).%s" % (
            common_name, prepend_identifier_if_necessary(identifier))
    else:
        src_var = conv_info.vmap[identifier]

    # COMMON scalar handled externally: keep it as a plain identifier and
    # do not emit local reference aliases like `T& x = cmn.x;`.
    if (fdecl.is_common()
            and getattr(fdecl, "dim_tokens", None) is None
            and identifier.lower() in FABLE_EXTERN_COMMON_SCALARS):
        conv_info.vmap[identifier] = prepend_identifier_if_necessary(identifier)
        return crhs is not None

    # For COMMON or SAVE variables, we sometimes create local references
    # like "REAL& x = cmn.x;" or "REAL& x = sve.x;".
    save_ok = (
        fdecl.is_save()
        and not (conv_info.fproc is not None
                 and getattr(conv_info.fproc, "conv_hook", None) is not None
                 and conv_info.fproc.conv_hook.ignore_common_and_save)
    )

    if (fdecl.dim_tokens is not None):
        conv_info.vmap[identifier] = prepend_identifier_if_necessary(
            identifier)
        if (fdecl.data_type.value == "character"):
            ctype = "str_arr_%sref<%d>" % (
                cconst(fdecl=fdecl, short=True), len(fdecl.dim_tokens))
        else:
            ctype = convert_data_type_and_dims(
                conv_info=conv_info, fdecl=fdecl, crhs=None, force_arr=True)[0]
            ctype = ad_hoc_change_arr_to_arr_ref(
                ctype=ctype, cconst=cconst(fdecl=fdecl, short=True))
        declare_size_dim_identifiers(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=top_scope,
            fdecl=fdecl)
        cdims = convert_dims(conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
        if (common_name is not None):
            src_var = "static_cast<common_%s&>(cmn).%s" % (
                common_name, prepend_identifier_if_necessary(identifier))
        rapp = get_rapp()
        rapp("%s %s(%s, %s);" % (
            ctype, prepend_identifier_if_necessary(identifier), src_var, cdims))
    elif (common_name is not None
          or (fdecl.use_count > 1
              and (fdecl.is_common() or save_ok))):
        conv_info.vmap[identifier] = prepend_identifier_if_necessary(
            identifier)
        ctype = convert_data_type(
            conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
        rapp = get_rapp()
        rapp("%s& %s = %s;" % (
            ctype, prepend_identifier_if_necessary(identifier), src_var))

    if (crhs is not None):
        return True
    return False


def convert_executable(
        callback, conv_info, args_fdecl_with_dim=None, blockdata=None):
    top_scope = scope(parent=None)
    # User policy: do not emit FEM_CMN_SVE(...) boilerplate.
    top_scope.remember_insert_point()
    curr_scope = top_scope
    if (args_fdecl_with_dim is not None):
        for fdecl in args_fdecl_with_dim:
            declare_size_dim_identifiers(
                conv_info=conv_info,
                top_scope=top_scope,
                curr_scope=top_scope,
                fdecl=fdecl)
            cdims = convert_dims(conv_info=conv_info,
                                 dim_tokens=fdecl.dim_tokens)
            top_scope.append("%s(%s);" % (fdecl.id_tok.value, cdims))
    if (blockdata is not None):
        for fproc in blockdata:
            callback("  %s(cmn);" % fproc.name.value)
    if (conv_info.fproc.uses_read):
        top_scope.append("common_read read(cmn);")
    if (conv_info.fproc.uses_write):
        top_scope.append("common_write write(cmn);")
    top_scope_point_before_common = top_scope.current_point()
    for common_name, fdecl_list in conv_info.fproc.common.items():
        if (common_name in conv_info.fproc.conv_hook.variant_common_names):
            continue
        top_scope.remember_insert_point()
        for common_fdecl in fdecl_list:
            fdecl = conv_info.fproc.fdecl_by_identifier.get(
                common_fdecl.id_tok.value)
            if (fdecl.use_count != 0
                    and fdecl.id_tok.value not in conv_info.vmap):
                declare_identifier(
                    conv_info=conv_info,
                    top_scope=top_scope,
                    curr_scope=top_scope,
                    id_tok=fdecl.id_tok)
        if (not top_scope.insert_point_is_current()):
            top_scope.top_append("// COMMON %s" % common_name)
    if (not top_scope.point_is_current(point=top_scope_point_before_common)):
        top_scope.append("//")
    top_scope.remember_insert_point()

    def declare_identifiers(id_tokens):
        for id_tok in id_tokens:
            if (id_tok.value not in conv_info.vmap):
                declare_identifier(
                    conv_info=conv_info,
                    top_scope=top_scope,
                    curr_scope=curr_scope,
                    id_tok=id_tok)
    from fable.tokenization import extract_identifiers

    # Hoist constant DATA initializers for local arrays as static declarations.
    # This avoids runtime DATA initialization blocks and keeps the static tables
    # grouped at the top of the function.
    try:
        if conv_info.array_data_initializers is None:
            conv_info.array_data_initializers = build_array_data_initializers(conv_info)
        if conv_info.array_data_initializers:
            for nlist, _clist in conv_info.fproc.data:
                if not nlist or len(nlist) != 1:
                    continue
                toks = getattr(nlist[0], "value", None)
                if not toks or len(toks) != 1 or not toks[0].is_identifier():
                    continue
                tgt = toks[0].value.lower()
                if tgt in conv_info.hoisted_data_array_names:
                    continue
                rec = conv_info.array_data_initializers.get(tgt)
                if rec is None:
                    continue
                elem_ctype, dims_ints, init_list, rank = rec
                fdecl = conv_info.fproc.get_fdecl(id_tok=toks[0])
                if fdecl is None:
                    continue
                # Only hoist local arrays.
                if not (fdecl.is_local() or fdecl.is_save() or fdecl.is_parameter()):
                    continue
                # DATA-initialized local arrays are safe to treat as function-local static.
                # Force local mapping to avoid invalid names like `sve.X` in declarations.
                conv_info.set_vmap_force_local(fdecl=fdecl)
                vname = conv_info.vmapped(fdecl=fdecl)
                mplapack_elem_ctype = convert_to_mplapack_type(elem_ctype)
                if rank == 1 and dims_ints and len(dims_ints) == 1:
                    dim_part = f"[{int(dims_ints[0])}]"
                else:
                    # For multi-dimensional DATA (e.g. IPIVOT(4,4)), emit an unsized initializer.
                    dim_part = "[]"
                top_scope.append(
                    "static %s %s%s = {%s};" % (
                        mplapack_elem_ctype,
                        vname,
                        dim_part,
                        ", ".join(init_list),
                    )
                )
                conv_info.hoisted_data_array_names.add(tgt)
    except Exception:
        pass

    # Hoist constant DATA initializers for local scalars as static declarations.
    # This covers patterns like:
    #   DATA threq / 2.0d0 / , intstr / '0123456789' /
    # even when runtime DATA init blocks are suppressed.
    try:
        if conv_info.data_initializers is None:
            conv_info.data_initializers = build_scalar_data_initializers(conv_info)
        if conv_info.data_initializers:
            for fdecl in conv_info.fproc.fdecl_by_identifier.values():
                if fdecl is None:
                    continue
                # Only scalars.
                if getattr(fdecl, "dim_tokens", None) is not None:
                    continue
                # Only local/SAVE/PARAMETER (skip COMMON).
                if not (fdecl.is_local() or fdecl.is_save() or fdecl.is_parameter()):
                    continue
                name_lower = fdecl.id_tok.value.lower()
                init_expr = conv_info.data_initializers.get(name_lower)
                if init_expr is None:
                    continue
                # Avoid redeclaration if already mapped/declared.
                if fdecl.id_tok.value in conv_info.vmap:
                    continue
                # Skip dummy CHARACTER arguments (already in signature).
                if _is_dummy_character_arg(conv_info, fdecl.id_tok.value):
                    continue

                # DATA implies SAVE for local variables: emit as function-local static.
                conv_info.set_vmap_force_local(fdecl=fdecl)
                vname = conv_info.vmapped(fdecl=fdecl)

                ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
                mplapack_ctype = convert_to_mplapack_type(ctype)
                # If the target is a plain char (CHARACTER*1), convert "X" -> 'X'.
                if mplapack_ctype == "char":
                    m = re.match(r'^"((?:\\.)|[^\\])"$', init_expr)
                    if m:
                        inner = m.group(1)
                        # Escape single quote inside the char literal.
                        if inner == "'":
                            inner = "\\'"
                        init_expr = "'" + inner + "'"

                top_scope.append(f"static {mplapack_ctype} {vname} = {init_expr};")
    except Exception:
        pass

    variant_buffers = convert_variant_allocate_and_bindings(
        conv_info=conv_info, top_scope=top_scope)
    top_scope.remember_insert_point()
    cei = conv_info.fproc.classified_equivalence_info()
    sve_equivalences = cei.save.equiv_tok_cluster_by_identifier
    for identifier in sorted(conv_info.fproc.fdecl_by_identifier.keys()):
        fdecl = conv_info.fproc.fdecl_by_identifier[identifier]
        if (fdecl.is_save()
                and fdecl.use_count > 1
                and not identifier in sve_equivalences
                and not identifier in conv_info.vmap):
            declare_identifier(
                conv_info=conv_info,
                top_scope=top_scope,
                curr_scope=top_scope,
                id_tok=fdecl.id_tok)
    if (not top_scope.insert_point_is_current()):
        top_scope.top_append("// SAVE")
        top_scope.append("//")
    # User policy: do not emit runtime DATA initialization blocks.
    if (False and conv_info.fproc.conv_hook.needs_is_called_first_time):
        first_time_scope = top_scope.open_nested_scope(
            opening_text=["if (is_called_first_time) {"])
        if (len(variant_buffers.first_time) != 0):
            first_time_scope.append(
                "using fem::mbr; // member of variant common or equivalence")
            for mbr_buffer in variant_buffers.first_time:
                mbr_scope = first_time_scope.open_nested_scope(
                    opening_text=["{"])
                for line in mbr_buffer:
                    mbr_scope.append(line)
                mbr_scope.close_nested_scope()
        for nlist, clist in conv_info.fproc.data:
            declare_identifiers(
                id_tokens=extract_identifiers(tokens=nlist))
            for repetition_tok, ctoks in clist:
                if (repetition_tok is not None):
                    declare_identifiers(
                        id_tokens=extract_identifiers(tokens=[repetition_tok]))
                declare_identifiers(
                    id_tokens=extract_identifiers(tokens=ctoks))
        if (not conv_info.fproc.conv_hook.data_init_after_variant_bind):
            convert_data(conv_info=conv_info, data_init_scope=first_time_scope)
        first_time_scope.close_nested_scope()
        top_scope.remember_insert_point()
    if (len(variant_buffers.loc_equivalences) != 0):
        top_scope.append("local_equivalences loc_equivalences;")
        mbr_scope = top_scope.open_nested_scope(opening_text=["{"])
        mbr_scope.append("using fem::mbr; // member")
        for line in variant_buffers.loc_equivalences:
            mbr_scope.append(line)
        mbr_scope.close_nested_scope()
    for line in variant_buffers.bindings:
        top_scope.append(line)
    top_scope.remember_insert_point()
    if (False and conv_info.fproc.conv_hook.data_init_after_variant_bind):
        data_init_scope = top_scope.open_nested_scope(
            opening_text=["if (is_called_first_time) {"])
        convert_data(conv_info=conv_info, data_init_scope=data_init_scope)
        data_init_scope.close_nested_scope()
    top_scope.remember_insert_point()
    from fable.tokenization import fmt_tokens_as_string

    def get_cfmt_from_format(stmt_label):
        fmt_tokens = conv_info.fproc.format.get(stmt_label)
        if (fmt_tokens is None):
            tok.raise_semantic_error(
                "Unknown FORMAT statement label: %s" % tok.value)
        return '"(' + escape_string_literal(fmt_tokens_as_string(
            tokens=fmt_tokens, comma=fmt_comma_placeholder)) + ')"'
    fmt_counts_by_statement_label = \
        conv_info.fproc.fmt_counts_by_statement_label()
    for stmt_label in sorted(fmt_counts_by_statement_label.keys()):
        if (fmt_counts_by_statement_label[stmt_label] > 1):
            cfmt = get_cfmt_from_format(stmt_label=stmt_label)
            top_scope.append(
                "static const char* format_%s = %s;" % (stmt_label, cfmt))

    def curr_scope_append_return_function():
        curr_scope.append(
            "return %s;" % conv_info.vmap[conv_info.fproc.name.value])
    close_scope_after_next_executable = False
    dos_to_close_by_label = {}
    from fable.read import Error
    from fable import SemanticError
    for ei in conv_info.fproc.executable:
        conv_info.comment_manager.insert_before(
            executable_info=ei, callback=curr_scope.append_comment)
        lbl = ei.ssl.label
        if (lbl is not None
                and lbl in conv_info.fproc.target_statement_labels()
                and not close_scope_after_next_executable):
            curr_scope.append_statement_label(label=lbl)

        def search_for_id_tokens_and_declare_identifiers():
            id_tokens = []

            def callback(tok, next_tok):
                id_tokens.append(tok)
            ei.search_for_id_tokens(callback=callback)
            declare_identifiers(id_tokens=id_tokens)
            return id_tokens
        try:
            if (ei.key == "assignment"):
                lhs_id_tokens = extract_identifiers(tokens=ei.lhs_tokens)
                assert len(lhs_id_tokens) != 0
                rhs_id_tokens = extract_identifiers(tokens=ei.rhs_tokens)
                for id_tokens in lhs_id_tokens[1:], rhs_id_tokens:
                    declare_identifiers(id_tokens=id_tokens)
                crhs = convert_tokens(
                    conv_info=conv_info, tokens=ei.rhs_tokens)
                # If LHS is REAL/DOUBLEPRECISION or INTEGER and RHS is a simple
                # COMPLEX or REAL variable/element, map according to Fortran
                # implicit conversion rules:
                #   REAL  = COMPLEX  -> REAL  = real(COMPLEX)
                #   INTEGER = COMPLEX -> INTEGER = castINTEGER(real(COMPLEX))
                #   INTEGER = REAL    -> INTEGER = castINTEGER(REAL)
                lhs_fdecl = conv_info.fproc.get_fdecl(id_tok=lhs_id_tokens[0])
                lhs_dt_code = None
                if lhs_fdecl is not None and lhs_fdecl.data_type is not None:
                    dt = lhs_fdecl.data_type
                    if isinstance(dt, str):
                        lhs_dt_code = dt
                    else:
                        lhs_dt_code = getattr(dt, "value", None)
                    if lhs_dt_code is not None:
                        lhs_dt_code = lhs_dt_code.lower()
                lhs_is_real = lhs_dt_code in ("real", "doubleprecision")
                lhs_is_integer = lhs_dt_code == "integer"
                lhs_is_character = lhs_dt_code == "character"
                lhs_is_complex = lhs_dt_code in ("complex", "doublecomplex")
                rhs_is_simple = _is_simple_lvalue(crhs)
                rhs_is_complex = False
                rhs_is_real = False

                # Determine the base identifier for type lookup
                # Handle both simple lvalues and unary-negated expressions
                rhs_base_id = None
                if rhs_is_simple and rhs_id_tokens:
                    rhs_base_id = rhs_id_tokens[0]
                elif rhs_id_tokens:
                    # Try to extract base identifier from unary expression
                    # e.g., "-savealpha" -> "savealpha"
                    sign, base_name, base_expr = _extract_base_identifier_from_unary(
                        crhs)
                    if base_name is not None:
                        # Find the corresponding id_token
                        for tok in rhs_id_tokens:
                            if tok.value.lower() == base_name.lower():
                                rhs_base_id = tok
                                break
                # If LHS is a scalar CHARACTER (mapped to 'char') and RHS is a
                # scalar CHARACTER dummy argument (mapped to (const) char*),
                # dereference RHS so that:
                #   char compq2 = compq;  ->  char compq2 = *compq;
                if lhs_is_character and rhs_base_id is not None:
                    if _is_dummy_character_arg(conv_info, rhs_base_id.value):
                        rhs_expr = crhs.strip()
                        # Only rewrite a bare identifier; do not touch substrings/expressions.
                        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", rhs_expr):
                            crhs = "*" + rhs_expr

                # Use fdecl information when available to classify RHS type
                if (lhs_is_real or lhs_is_integer) and rhs_base_id is not None:
                    rhs_fdecl = conv_info.fproc.get_fdecl(id_tok=rhs_base_id)
                    rhs_dt_code = None
                    if rhs_fdecl is not None and rhs_fdecl.data_type is not None:
                        dt = rhs_fdecl.data_type
                        if isinstance(dt, str):
                            rhs_dt_code = dt
                        else:
                            rhs_dt_code = getattr(dt, "value", None)
                        if rhs_dt_code is not None:
                            rhs_dt_code = rhs_dt_code.lower()
                            rhs_is_complex = rhs_dt_code in (
                                "complex", "doublecomplex")
                            rhs_is_real = rhs_dt_code in (
                                "real", "doubleprecision")

                # Additional heuristic: detect simple complex-valued expressions of the form
                #   a[...] op ...
                # where 'a' is a COMPLEX variable and op is one of +,-,*,/.
                # This is needed for patterns such as:
                #   d11 = a[(k - 1) + (k - 1)*lda] / d;
                # where d11 is REAL and a is COMPLEX.
                if (lhs_is_real or lhs_is_integer) and not rhs_is_complex:
                    rhs_expr = crhs.strip()
                    m_leading = re.match(
                        r"([A-Za-z_][A-Za-z0-9_]*)\s*(\[[^\]]*\])?\s*([-+*/])",
                        rhs_expr,
                    )
                    if m_leading:
                        base_name = m_leading.group(1)
                        if (base_name in complex_identifiers
                                or base_name in complex_pointer_identifiers):
                            rhs_is_complex = True
                if lhs_is_real and rhs_is_complex:
                    rhs_expr = crhs.strip()
                    # Avoid double-wrapping if someone already wrote .real() / .imag()
                    if ".real()" not in rhs_expr and ".imag()" not in rhs_expr:
                        if _is_simple_lvalue(rhs_expr):
                            crhs = f"{rhs_expr}.real()"
                        else:
                            crhs = f"({rhs_expr}).real()"

                elif lhs_is_integer and rhs_is_complex:
                    rhs_expr = crhs.strip()
                    # INTEGER = COMPLEX -> INTEGER = castINTEGER(real(COMPLEX))
                    if ".real()" not in rhs_expr and ".imag()" not in rhs_expr:
                        if _is_simple_lvalue(rhs_expr):
                            real_expr = f"{rhs_expr}.real()"
                        else:
                            real_expr = f"({rhs_expr}).real()"
                    else:
                        # Already has .real()/.imag(), just reuse
                        real_expr = rhs_expr
                    crhs = f"castINTEGER({real_expr})"

                elif lhs_is_integer and rhs_is_real:
                    rhs_expr = crhs.strip()
                    # INTEGER = REAL -> INTEGER = castINTEGER(REAL)
                    if "castINTEGER(" in rhs_expr:
                        # Already converted to INTEGER somewhere else
                        pass
                    else:
                        if _is_simple_lvalue(rhs_expr):
                            crhs = f"castINTEGER({rhs_expr})"
                        else:
                            crhs = f"castINTEGER({rhs_expr})"

                # Do not apply any additional coercion for REAL left-hand sides here.
                # All REAL/DBLE handling should be driven by:
                #   - explicit REAL/DBLE intrinsics (handled in rewrite_intrinsics)
                #   - simple REAL = COMPLEX assignments handled above.

                    if rhs_expr:
                        # Expressions that are clearly real-valued already:
                        #   - known real-valued function calls
                        #   - explicit castREAL(...)
                        #   - explicit .real() / .imag() projections
                        rhs_obviously_real = (
                            _is_real_valued_function_call(rhs_expr)
                            or "castREAL(" in rhs_expr
                            or ".real()" in rhs_expr
                            or ".imag()" in rhs_expr
                        )

                        # If the RHS still is not obviously REAL and it contains
                        # at least one COMPLEX-valued function call (e.g. CDOTC,
                        # ZDOTC), emulate Fortran's implicit REAL = COMPLEX
                        # coercion by taking the real part of the entire RHS.
                        if (not rhs_obviously_real
                                and _contains_complex_valued_function_call(rhs_expr)):
                            crhs = _real_cast_or_component(
                                rhs_expr,
                                complex_identifiers,
                                complex_pointer_identifiers,
                            )
                # Promote pure integer literals to real literals when assigning
                # to REAL/DOUBLEPRECISION *or* COMPLEX variables, so that:
                #   t[...] = 0; becomes t[...] = 0.0;
                # This helps overload resolution for MPLAPACK scalar types.
                if lhs_is_real or lhs_is_complex:
                    rhs_expr = crhs.strip()
                    # Match a pure integer literal (e.g. 0, -1, +2)
                    if re.fullmatch(r"[+-]?[0-9]+", rhs_expr):
                        # Normalize zero to 0.0
                        if rhs_expr in ("0", "+0", "-0"):
                            crhs = "0.0"
                        else:
                            # e.g. 1 -> 1.0, -2 -> -2.0
                            crhs = rhs_expr + ".0"

                id_tok = lhs_id_tokens[0]

                assign_here = id_tok.value in conv_info.vmap
                if (not assign_here):
                    assign_here = declare_identifier(
                        conv_info=conv_info,
                        top_scope=top_scope,
                        curr_scope=curr_scope,
                        id_tok=id_tok,
                        crhs=crhs)
                clhs = convert_tokens(
                    conv_info=conv_info, tokens=ei.lhs_tokens)

                # For scalar CHARACTER dummy arguments emitted as plain (const) char*,
                # assignments should dereference the pointer:
                #   EQUED = 'N'  ->  *equed = 'N';
                #
                # In view mode (FABLE_SMALL_CHAR=0), scalar CHARACTER dummies are
                # emitted as fem::str_ref / fem::str_cref, so dereferencing would
                # be incorrect (e.g. dist = 'S' must stay as `dist = "S";`).
                if (
                    lhs_is_character
                    and lhs_fdecl is not None
                    and lhs_fdecl.dim_tokens is None
                    and _is_dummy_character_arg(conv_info, id_tok.value)
                    and _is_plain_character_pointer_dummy(conv_info, id_tok.value)
                ):
                    # clhs is typically just the mapped name (e.g. "equed")
                    if not clhs.startswith("*"):
                        clhs = "*" + clhs

                if (assign_here):
                    def in_place_op_left():
                        if (not crhs.startswith(clhs)):
                            return False
                        i = len(clhs)
                        if (i == len(crhs)):
                            return False
                        if (crhs[i] != " "):
                            return False
                        i += 1
                        if (i == len(crhs)):
                            return False
                        op = crhs[i]
                        if (op != "+"):
                            return False
                        i += 1
                        if (i == len(crhs)):
                            return False
                        if (crhs[i] != " "):
                            return False
                        i += 1
                        if (i == len(crhs)):
                            return False
                        if (crhs[i:] == "1"):
                            curr_scope.append("%s++;" % clhs)
                        else:
                            curr_scope.append("%s %s= %s;" %
                                              (clhs, op, crhs[i:]))
                        return True

                    def in_place_op_right():
                        if (not crhs.endswith(clhs)):
                            return False
                        i = len(crhs) - len(clhs)
                        if (i == 0):
                            return False
                        i -= 1
                        if (crhs[i] != " "):
                            return False
                        if (i == 0):
                            return False
                        i -= 1
                        op = crhs[i]
                        if (op != "+"):
                            return False
                        if (i == 0):
                            return False
                        i -= 1
                        if (crhs[i] != " "):
                            return False
                        if (i == 0):
                            return False
                        if (crhs[:i] == "1"):
                            curr_scope.append("%s++;" % clhs)
                        else:
                            curr_scope.append("%s %s= %s;" %
                                              (clhs, op, crhs[:i]))
                        return True
                    if (not in_place_op_left() and not in_place_op_right()):
                        if _emit_small_char_string_assignment(
                                curr_scope=curr_scope, clhs=clhs, crhs=crhs, conv_info=conv_info):
                            pass
                        else:
                            # Rewrite substring-style calls on small CHARACTER*n scalars
                            # that were mapped to char arrays, e.g.:
                            #   jbcmpz(1,1) = 'S'  ->  jbcmpz[0] = 'S';
                            clhs_fixed = _rewrite_small_char_substrings(clhs)
                            crhs_fixed = _rewrite_small_char_substrings(crhs)
                            curr_scope.append("%s = %s;" %
                                              (clhs_fixed, crhs_fixed))

            elif (ei.key == "inquire"):
                search_for_id_tokens_and_declare_identifiers()
                iuflist = ei.iuflist
                if (iuflist.unit is not None):
                    if (iuflist.file is not None):
                        ei.ssl.raise_semantic_error(
                            "Conflicting UNIT vs. FILE in INQUIRE statement"
                            " (exactly one is needed)", i=ei.start)
                    io_function_specialization = "_unit"
                    uf_tokens = iuflist.unit
                elif (iuflist.file is not None):
                    io_function_specialization = "_file"
                    uf_tokens = iuflist.file
                else:
                    ei.ssl.raise_semantic_error(
                        "Missing UNIT or FILE in INQUIRE statement", i=ei.start)
                io_call_args = convert_tokens(
                    conv_info=conv_info, tokens=uf_tokens)
                convert_io_statement_with_err(
                    conv_info=conv_info,
                    curr_scope=curr_scope,
                    io_function="inquire",
                    io_function_specialization=io_function_specialization,
                    io_call_args=io_call_args,
                    iolist=iuflist)
            elif (ei.key == "file_positioning"):
                search_for_id_tokens_and_declare_identifiers()
                io_call_args = convert_tokens(
                    conv_info=conv_info, tokens=ei.alist.unit)
                convert_io_statement_with_err(
                    conv_info=conv_info,
                    curr_scope=curr_scope,
                    io_function=ei.io_function,
                    io_function_specialization="",
                    io_call_args=io_call_args,
                    iolist=ei.alist)
            elif (ei.key == "open"):
                search_for_id_tokens_and_declare_identifiers()
                olist = ei.olist
                if (olist.unit is None):
                    ei.ssl.raise_semantic_error(
                        "Missing UNIT in OPEN statement", i=ei.start)
                cunit = convert_tokens(conv_info=conv_info, tokens=olist.unit)
                if (olist.file is None):
                    cfile = "fem::file_not_specified"
                else:
                    cfile = convert_tokens(
                        conv_info=conv_info, tokens=olist.file)
                convert_io_statement_with_err(
                    conv_info=conv_info,
                    curr_scope=curr_scope,
                    io_function="open",
                    io_function_specialization="",
                    io_call_args="%s, %s" % (cunit, cfile),
                    iolist=olist)
            elif (ei.key == "close"):
                search_for_id_tokens_and_declare_identifiers()
                cllist = ei.cllist
                if (cllist.unit is None):
                    ei.ssl.raise_semantic_error(
                        "Missing UNIT in CLOSE statement", i=ei.start)
                cunit = convert_tokens(conv_info=conv_info, tokens=cllist.unit)
                convert_io_statement_with_err(
                    conv_info=conv_info,
                    curr_scope=curr_scope,
                    io_function="close",
                    io_function_specialization="",
                    io_call_args=cunit,
                    iolist=cllist)
            elif (ei.key in ["read", "write", "print"]):
                search_for_id_tokens_and_declare_identifiers()
                cilist = ei.cilist
                if (ei.key == "print"):
                    work_key = "write"
                    cunit = "6"
                else:
                    work_key = ei.key
                    assert cilist.unit is not None
                    cunit = convert_tokens(
                        conv_info=conv_info, tokens=cilist.unit)
                    if (cunit == "star "):
                        cunit = "6"

                def conv_fmt():
                    if (ei.fmt_tokens is not None):
                        return '"(' + escape_string_literal(fmt_tokens_as_string(
                            tokens=ei.fmt_tokens, comma=fmt_comma_placeholder)) + ')"'
                    tl = cilist.fmt
                    if (tl is None):
                        return None
                    if (len(tl) == 1):
                        tok = tl[0]
                        if (tok.is_op_with(value="*")):
                            return "star"
                        if (tok.is_integer()):
                            stmt_label = tok.value
                            if (fmt_counts_by_statement_label[stmt_label] > 1):
                                return "format_%s" % stmt_label
                            return get_cfmt_from_format(stmt_label=stmt_label)
                    return convert_tokens(conv_info=conv_info, tokens=tl)
                cfmt = conv_fmt()
                cchain = []
                has_iostat = False
                for slot in ["rec", "iostat"]:
                    tokens = getattr(cilist, slot)
                    if (tokens is not None):
                        cchain.append("%s(%s)" % (
                            slot, convert_tokens(conv_info=conv_info, tokens=tokens)))
                        if slot == "iostat":
                            has_iostat = True
                if (len(cchain) == 0):
                    cchain = ""
                else:
                    cchain = "." + ".".join(cchain)
                iolist_id_tokens = extract_identifiers(tokens=ei.iolist)
                declare_identifiers(id_tokens=iolist_id_tokens)
                if (cfmt is None):
                    cargs = "%s, fem::unformatted" % cunit
                else:
                    cargs = "%s, %s" % (cunit, cfmt)
                implied_dos = []
                find_implied_dos(result=implied_dos, tokens=ei.iolist)
                if (len(implied_dos) == 0):
                    if (cilist.end is None
                            and cilist.err is None
                            and not has_iostat):
                        io_scope = curr_scope
                    else:
                        io_scope = curr_scope.open_nested_scope(
                            opening_text=["try {"])
                    io_op = "%s(%s)%s" % (work_key, cargs, cchain)
                    if (len(ei.iolist) == 0):
                        io_scope.append(io_op+";")
                    else:
                        convert_io_loop(
                            io_scope=io_scope,
                            io_op=io_op,
                            conv_info=conv_info,
                            tokens=ei.iolist)
                else:
                    is_internal_file = False
                    if (cilist.unit is not None):
                        unit_id_tokens = extract_identifiers(
                            tokens=cilist.unit)
                        if (len(unit_id_tokens) >= 1):
                            unit_fdecl = conv_info.fproc.get_fdecl(
                                id_tok=unit_id_tokens[0])
                            if (unit_fdecl.data_type is not None
                                    and unit_fdecl.data_type.value == "character"):
                                is_internal_file = True
                    if (cilist.end is not None
                        or cilist.err is not None
                            or has_iostat):
                        opening_line = "try {"
                    else:
                        opening_line = "{"
                    io_scope = curr_scope.open_nested_scope(
                        opening_text=[opening_line])
                    if (is_internal_file):
                        cmn = ""
                    else:
                        cmn = "cmn, "
                    io_scope.append("%s_loop %sloop(%s%s);" % (
                        work_key, work_key[0], cmn, cargs))
                    if (len(cchain) != 0):
                        io_scope.append("%sloop%s;" % (work_key[0], cchain))
                    convert_io_loop(
                        io_scope=io_scope,
                        io_op=work_key[0]+"loop",
                        conv_info=conv_info,
                        tokens=ei.iolist)
                if (io_scope is not curr_scope):
                    io_scope.close_nested_scope()
                for slot, cexception in [("end", "read_end"), ("err", "io_err")]:
                    tokens = getattr(cilist, slot)
                    if (tokens is not None):
                        catch_scope = curr_scope.open_nested_scope(
                            opening_text=["catch (fem::%s const&) {" % cexception])
                        clabel = convert_tokens(
                            conv_info=conversion_info(), tokens=tokens)
                        catch_scope.append("goto statement_%s;" % clabel)
                        catch_scope.close_nested_scope()
                    else:
                        if has_iostat:
                            catch_scope = curr_scope.open_nested_scope(
                                opening_text=["catch (fem::%s const&) {" % cexception])
                            catch_scope.close_nested_scope()
            elif (ei.key == "do"):
                if (ei.id_tok.value not in conv_info.vmap):
                    declare_identifier(
                        conv_info=conv_info,
                        top_scope=top_scope,
                        curr_scope=curr_scope,
                        id_tok=ei.id_tok)
                for token in ei.tokens:
                    declare_identifiers(
                        id_tokens=extract_identifiers(tokens=token.value))
                curr_scope = convert_to_fem_do(
                    conv_info=conv_info,
                    parent_scope=curr_scope,
                    i_tok=ei.id_tok,
                    fls_tokens=ei.tokens)
                if (ei.label is not None):
                    dos_to_close_by_label.setdefault(
                        ei.label, []).append(curr_scope)
            elif (ei.key == "dowhile"):
                declare_identifiers(
                    id_tokens=extract_identifiers(tokens=ei.cond_tokens))
                c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
                curr_scope = curr_scope.open_nested_scope(
                    opening_text=["while %s {" % c])
                if (ei.label is not None):
                    dos_to_close_by_label.setdefault(
                        ei.label, []).append(curr_scope)
            elif (ei.key == "cycle"):
                curr_scope.append("continue;")
            elif (ei.key == "exit"):
                curr_scope.append("break;")
            elif (ei.key == "enddo"):
                if (dos_to_close_by_label.get(ei.ssl.label) is None):
                    curr_scope = curr_scope.close_nested_scope()
            elif (ei.key == "if"):
                declare_identifiers(
                    id_tokens=extract_identifiers(tokens=ei.cond_tokens))
                c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
                curr_scope = curr_scope.open_nested_scope(
                    opening_text=["if (%s) {" % c])
                close_scope_after_next_executable = True
                continue
            elif (ei.key == "if_then"):
                declare_identifiers(
                    id_tokens=extract_identifiers(tokens=ei.cond_tokens))
                c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
                curr_scope = curr_scope.open_nested_scope(
                    opening_text=["if (%s) {" % c])
            elif (ei.key == "elseif_then"):
                declare_identifiers(
                    id_tokens=extract_identifiers(tokens=ei.cond_tokens))
                c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)

                # Robust handling: the upstream reader can occasionally misclassify
                # an IF(... ) THEN block as a single-statement IF, so ELSEIF may
                # appear without an open IF scope here. Do not crash; recover and
                # continue generating output.
                def _scope_is_if_or_elseif(sc):
                    if sc is None or sc.opening_text is None or len(sc.opening_text) == 0:
                        return False
                    head = sc.opening_text[0].lstrip()
                    return head.startswith("if (") or head.startswith("else if (")

                attach_scope = curr_scope
                while (attach_scope is not None
                       and (not _scope_is_if_or_elseif(attach_scope)
                            or attach_scope.closing_text is not None
                            or attach_scope.tail is not None)):
                    attach_scope = attach_scope.parent

                if attach_scope is None:
                    curr_scope.append(
                        "// UNHANDLED: ELSEIF without matching IF; treating as IF")
                    curr_scope = curr_scope.open_nested_scope(
                        opening_text=["if (%s) {" % c])
                else:
                    while (curr_scope is not attach_scope
                           and curr_scope.parent is not None
                           and curr_scope.opening_text is not None):
                        curr_scope = curr_scope.close_nested_scope()
                    curr_scope = attach_scope.attach_tail(
                        opening_text=["else if (%s) {" % c])

            elif (ei.key == "else"):
                # ELSE must attach to a currently-open IF/ELSEIF scope.
                def _scope_is_if_or_elseif(sc):
                    if sc is None or sc.opening_text is None or len(sc.opening_text) == 0:
                        return False
                    head = sc.opening_text[0].lstrip()
                    return head.startswith("if (") or head.startswith("else if (")

                attach_scope = curr_scope
                while (attach_scope is not None
                       and (not _scope_is_if_or_elseif(attach_scope)
                            or attach_scope.closing_text is not None
                            or attach_scope.tail is not None)):
                    attach_scope = attach_scope.parent

                if attach_scope is None:
                    curr_scope.append(
                        "// UNHANDLED: ELSE without matching IF; emitting as a standalone block")
                    curr_scope = curr_scope.open_nested_scope(
                        opening_text=["{"])
                else:
                    while (curr_scope is not attach_scope
                           and curr_scope.parent is not None
                           and curr_scope.opening_text is not None):
                        curr_scope = curr_scope.close_nested_scope()
                    curr_scope = attach_scope.attach_tail(
                        opening_text=["else {"])

            elif (ei.key == "endif"):
                # Be tolerant of malformed IF/ENDIF structure.
                if (curr_scope.parent is None):
                    curr_scope.append(
                        "// UNHANDLED: ENDIF without matching IF")
                else:
                    curr_scope = curr_scope.close_nested_scope()

            elif (ei.key == "if_arithmetic"):
                declare_identifiers(
                    id_tokens=extract_identifiers(tokens=ei.cond_tokens))
                c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
                curr_scope = curr_scope.open_nested_scope(
                    opening_text=["switch (fem::if_arithmetic(%s)) {" % c])

                def lbl(i): return "statement_" + ei.labels[i].value
                curr_scope.append("case -1: goto %s;" % lbl(0))
                curr_scope.append("case  0: goto %s;" % lbl(1))
                curr_scope.append("default: goto %s;" % lbl(2))
                curr_scope = curr_scope.close_nested_scope()
            elif (ei.key == "call"):
                fdecl = conv_info.fproc.get_fdecl(id_tok=ei.subroutine_name)
                if (fdecl.is_intrinsic()):
                    from fable import intrinsics
                    cmn = ""
                    if (ei.subroutine_name.value == "getarg"):
                        called = "cmn.getarg"
                    elif (ei.subroutine_name.value in intrinsics.io_set_lower):
                        called = "cmn.io.%s" % ei.subroutine_name.value
                    else:
                        called = "fem::%s" % ei.subroutine_name.value
                else:
                    if (called_fproc_needs_cmn(
                        conv_info=conv_info,
                            called_name=ei.subroutine_name.value)):
                        cmn = "cmn"
                    else:
                        cmn = ""
                    called = conv_info.vmapped_callable(
                        identifier=ei.subroutine_name.value)
                    # Rename BLAS/LAPACK helper calls (case-insensitive)
                    simple_name = called.split("::")[-1]
                    lower_name = simple_name.lower()
                    prefix = called[:-len(simple_name)]
                    if lower_name == "lsame":
                        called = prefix + "Mlsame"
                    elif lower_name == "xerbla":
                        called = prefix + "Mxerbla"
                    else:
                        # General MPLAPACK renaming for BLAS/LAPACK-style routines
                        mpl_name = convert_function_name_to_mplapack(
                            simple_name)
                        called = prefix + mpl_name

                if (ei.arg_token is None):
                    curr_scope.append("%s(%s);" % (called, cmn))
                else:
                    id_tokens = extract_identifiers(tokens=ei.arg_token.value)
                    declare_identifiers(id_tokens=id_tokens)
                    a = convert_tokens(
                        conv_info=conv_info, tokens=ei.arg_token.value, commas=True)

                    # Use the callee name (without namespace) as the lookup key.
                    # For example, "Rlasr" -> "rlasr".
                    callee_key = called.split("::")[-1]
                    sig = _lookup_routine_signature(callee_key)
                    if sig is not None:
                        force_elems = bool(FABLE_SMALL_CHAR_VIEW and callee_key.lower() in _MPLAPACK_CORE_CPP_NAMES)
                        a = _adjust_actuals_using_signature(a, sig, conv_info, force_elems_call=force_elems)

                    def cmn_a():
                        if (len(cmn) == 0):
                            return a
                        if (len(a) == 0):
                            return cmn
                        return cmn + ", " + a
                    curr_scope.append("%s(%s);" % (called, cmn_a()))

            elif (ei.key == "return"):
                if (conv_info.fproc.fproc_type == "function"):
                    curr_scope_append_return_function()
                elif (ei is not conv_info.fproc.executable[-1]):
                    curr_scope.append("return;")
            elif (ei.key == "continue"):
                pass
            elif (ei.key == "goto"):
                curr_scope.append("goto statement_%s;" % ei.label.value)
            elif (ei.key == "goto_computed"):
                search_for_id_tokens_and_declare_identifiers()
                ccond = convert_tokens(conv_info=conv_info, tokens=ei.tokens)
                switch_scope = curr_scope.open_nested_scope(
                    opening_text=["switch (%s) {" % ccond])
                for i, label in enumerate(ei.labels):
                    switch_scope.append(
                        "case %d: goto statement_%s;" % (i+1, label.value))
                switch_scope.append("default: break;")
                switch_scope.close_nested_scope()
            elif (ei.key == "stop"):
                if (ei.arg_token is None):
                    cmsg = "0"
                elif (ei.arg_token.is_integer()):
                    cmsg = strip_leading_zeros(string=ei.arg_token.value)
                else:
                    cmsg = convert_token(
                        vmap={}, leading=True, tok=ei.arg_token)
                curr_scope.append("FEM_STOP(%s);" % cmsg)
            elif (ei.key == "entry"):
                curr_scope.append(
                    "// UNHANDLED: ENTRY %s" % ei.ssl.code_with_strings()[5:])
            else:
                curr_scope.append(
                    'FEM_THROW_UNHANDLED("executable %s: %s");' % (
                        ei.key, ei.ssl.code_with_strings()))
            if (close_scope_after_next_executable):
                close_scope_after_next_executable = False
                curr_scope = curr_scope.close_nested_scope()
            if (ei.ssl.label is not None):
                dos_to_close = dos_to_close_by_label.get(ei.ssl.label)
                if (dos_to_close is not None):
                    for do_scope in reversed(dos_to_close):
                        curr_scope = do_scope.close_nested_scope()
        except (Error, SemanticError):
            raise
        except Exception:
            print("*"*80)
            print(ei.ssl.format_error(
                i=None,
                msg="Sorry: fable internal error"))
            print("*"*80)
            print()
            raise
    # Emit any constant leading-dimension variables requested during conversion.
    _emit_constant_ld_decls(top_scope, conv_info)

    assert curr_scope.parent is None
    if (conv_info.fproc.fproc_type == "function"
            and len(conv_info.fproc.executable) != 0
            and conv_info.fproc.executable[-1].key != "return"):
        curr_scope_append_return_function()
    conv_info.comment_manager.flush_remaining(
        callback=curr_scope.append_comment)
    curr_scope.finalize()
    curr_scope.collect_text(callback=callback)


def export_save_struct(callback, conv_info):
    # User policy: when COMMON/SAVE boilerplate is suppressed, do not emit
    # any auto-generated SAVE structs. These will be provided manually.
    if FABLE_SUPPRESS_COMMON:
        return
    cci = conv_info.converted_commons_info
    if (cci is not None):
        buffer = cci.save_struct_buffers.get(conv_info.fproc.name.value)
        if (buffer is not None):
            for line in buffer:
                callback(line)


def produce_fortran_file_comment(conv_info, callback):
    if (conv_info.fortran_file_comments):
        callback("// Fortran file: %s"
                 % conv_info.fproc.body_lines[0].source_line_cluster[0].file_name)


def _mark_call_actuals_used(conv_info):
    """Bump use_count for variables that appear as actual arguments in CALL.

    Some local workspace arrays may only be used as actual arguments to
    external routines, and their use_count may be zero from the expression
    analysis. This pass ensures that any such variable is treated as "used"
    at least once so that its declaration is not dropped.
    """
    from fable.tokenization import extract_identifiers
    import sys

    # conv_info.fproc is the current Fortran procedure
    fproc = getattr(conv_info, "fproc", None)
    if fproc is None:
        return

    # Executable statements are stored in fproc.executable
    execs = getattr(fproc, "executable", None)
    if execs is None:
        return

    for ei in execs:
        # We only care about CALL statements
        if getattr(ei, "key", None) != "call":
            continue

        # The tokens for the actual argument list live in ei.arg_token.value
        arg_tok = getattr(ei, "arg_token", None)
        if arg_tok is None:
            continue

        id_tokens = extract_identifiers(tokens=arg_tok.value)
        for id_tok in id_tokens:
            fdecl = fproc.get_fdecl(id_tok=id_tok)
            if fdecl is None:
                continue
            # Only bump use_count for local variables; skip dummies, COMMON, etc.
            if not fdecl.is_local():
                continue
            if getattr(fdecl, "use_count", 0) == 0:
                # Debug print; comment out when no longer needed
                print(
                    f"[CALL_USED] bump use_count for {id_tok.value}", file=sys.stderr)
                fdecl.use_count = 1


def _sig_kind_requires_mutable_actual(kind: str) -> bool:
    """Return True if a callee parameter kind requires a mutable (non-const) actual.

    This is used to propagate OUT/INOUT-ness across calls so we don't
    generate "const" or by-value dummy arguments that cannot be forwarded.

    The signature kinds come from mplapack_signatures.FUNCTION_SIGNATURES.
    We keep this intentionally conservative: any non-const pointer or
    non-const reference requires a mutable actual.
    """
    if kind is None:
        return False
    k = str(kind).strip().upper()

    # Scalar passed by non-const reference.
    if k == "REF_SCALAR":
        return True

    # Character pointers: distinguish input vs output/inout explicitly.
    # - PTR_CHAR_IN  : const char*  (input-only) -> does NOT require mutable actual
    # - PTR_CHAR_OUT : char*        (output/inout) -> requires mutable actual
    if k == "PTR_CHAR_IN":
        return False
    if k == "PTR_CHAR_OUT":
        return True

    # Legacy/ambiguous char pointer kind: do NOT force mutability.
    # This avoids incorrectly promoting input-only flags like TRANS.
    if k == "PTR_CHAR":
        return False

    # Numeric pointers are always treated as mutable in our generated interfaces.
    if k == "PTR_NUMERIC":
        return True

    # Other pointer-like kinds: keep conservative behavior, honoring explicit CONST markers.
    if k.startswith("PTR_"):
        if "CONST" in k:
            return False
        return True

    return False


def _split_call_actuals_tokens(arg_tokens):
    """Split a CALL actual-argument token list into per-argument token lists."""
    if arg_tokens is None:
        return []

    # In fable's token stream, comma-separated lists are represented as
    # 'seq' tokens. This mirrors convert_tokens(..., commas=True).
    try:
        from fable.tokenization import group_power
    except Exception:
        return [arg_tokens] if arg_tokens else []

    actuals = []
    for tok in group_power(tokens=arg_tokens):
        if tok.is_seq():
            actuals.append(tok.value)

    # Fallback: if we did not see any seq token, treat the whole list as one.
    if not actuals and arg_tokens:
        actuals = [arg_tokens]

    return actuals


def _actual_base_identifier_if_definable(fproc, actual_tokens):
    """Return base identifier name if the actual is a definable variable.

    We recognize:
      - a plain identifier: X
      - an array element or substring-like reference: A(i), A(i,j), C(1:1)
        (only if A/C is an array/CHARACTER variable, not a function call)

    For propagation we only need the base variable name.
    """
    if fproc is None or not actual_tokens:
        return None

    # Plain identifier
    if len(actual_tokens) == 1 and actual_tokens[0].is_identifier():
        return actual_tokens[0].value

    # Identifier followed by parentheses -> array element / substring.
    if (len(actual_tokens) >= 2
            and actual_tokens[0].is_identifier()
            and actual_tokens[1].is_parentheses()):
        id_tok = actual_tokens[0]
        try:
            fdecl = fproc.get_fdecl(id_tok=id_tok)
        except Exception:
            fdecl = None
        if fdecl is None:
            return None
        # Reject callable identifiers: foo(x) is not a definable variable.
        if getattr(fdecl, "is_user_defined_callable", lambda: False)():
            return None
        if getattr(fdecl, "dim_tokens", None) is not None:
            return id_tok.value
        # CHARACTER substring on a scalar CHARACTER is also definable.
        dt = getattr(fdecl, "data_type", None)
        if dt is not None and getattr(dt, "value", None) == "character":
            return id_tok.value

    return None


def _callee_mutable_mask_from_fproc(callee_fproc):
    """Return a list[bool] saying whether each callee argument needs a mutable actual."""
    mask = []
    if callee_fproc is None:
        return mask

    for id_tok in getattr(callee_fproc, "args", []) or []:
        if getattr(id_tok, "value", None) == "*":
            # Alternate return marker in old Fortran. Not an OUT/INOUT variable.
            mask.append(False)
            continue

        try:
            fdecl = callee_fproc.get_fdecl(id_tok=id_tok)
        except Exception:
            fdecl = None
        if fdecl is None:
            mask.append(False)
            continue

        # EXTERNAL callable dummies are function pointers.
        if getattr(fdecl, "is_user_defined_callable", lambda: False)():
            mask.append(False)
            continue

        dt = getattr(fdecl, "data_type", None)
        dt_code = getattr(dt, "value", None) if dt is not None else None
        is_char = (dt_code == "character")
        is_array = (getattr(fdecl, "dim_tokens", None) is not None)

        # --- Arrays ---
        if is_array:
            # Non-character arrays are emitted as raw pointers (TYPE*), i.e. non-const.
            # To keep the generated C++ compilable, treat them as requiring mutable
            # actuals (even if the Fortran INTENT would be IN).
            if not is_char:
                mask.append(True)
            else:
                # CHARACTER arrays are emitted as str_arr_cref/ref depending on is_modified.
                mask.append(bool(getattr(fdecl, "is_modified", False)))
            continue

        # --- Scalars ---
        # Numeric/logical scalars: non-const reference iff modified.
        # CHARACTER scalars: char* iff modified.
        mask.append(bool(getattr(fdecl, "is_modified", False)))

    return mask


def _propagate_out_inout_through_calls(topological_fprocs, *, max_rounds: int = 1000) -> None:
    """Propagate OUT/INOUT-ness from callees to callers.

    Rule:
      If a CALL passes an actual variable into a callee parameter that
      requires a mutable (non-const) lvalue in the generated C++ signature,
      then that actual variable must be treated as modified in the caller.

    This is crucial because we intentionally generate IN scalar dummies as
    pass-by-value 'TYPE const' to match MPLAPACK style; but once such a dummy
    is forwarded to an OUT/INOUT position, it must become non-const and
    non-value (i.e. '&' or pointer) at the caller boundary.

    This pass is monotone (False -> True only) and converges.
    """
    if topological_fprocs is None:
        return

    # Map callee name -> fproc (case-insensitive)
    fprocs_by_name = getattr(
        getattr(topological_fprocs, "all_fprocs", None), "fprocs_by_name", None)
    if callable(fprocs_by_name):
        fproc_map = fprocs_by_name()
    else:
        fproc_map = {}
    fproc_map_lower = {k.lower(): v for (k, v) in (fproc_map or {}).items()}

    rounds = 0
    changed = True

    while changed:
        rounds += 1
        if rounds > max_rounds:
            # Defensive: should never happen because we only flip bits to True.
            break

        changed = False

        for caller in getattr(topological_fprocs, "bottom_up_list", []) or []:
            execs = getattr(caller, "executable", None) or []
            for ei in execs:
                if getattr(ei, "key", None) != "call":
                    continue

                callee_name = getattr(
                    getattr(ei, "subroutine_name", None), "value", None)
                if not callee_name:
                    continue

                # Determine which callee argument positions require mutable actuals.
                callee_fproc = fproc_map_lower.get(str(callee_name).lower())

                if callee_fproc is not None:
                    mutable_mask = _callee_mutable_mask_from_fproc(
                        callee_fproc)
                else:
                    # External routine: consult mplapack_signatures if available.
                    sig = _lookup_routine_signature(str(callee_name))
                    if sig is None:
                        continue
                    mutable_mask = [
                        _sig_kind_requires_mutable_actual(k) for k in sig]

                # Collect actual argument token lists.
                arg_tok = getattr(ei, "arg_token", None)
                if arg_tok is None:
                    actuals = []
                else:
                    actuals = _split_call_actuals_tokens(
                        getattr(arg_tok, "value", None))

                if len(actuals) != len(mutable_mask):
                    # Be conservative: don't guess when counts mismatch.
                    continue

                # Mark caller-side variables as modified where needed.
                for idx, needs_mutable in enumerate(mutable_mask):
                    if not needs_mutable:
                        continue
                    base = _actual_base_identifier_if_definable(
                        caller, actuals[idx])
                    if base is None:
                        continue

                    # Look up the declaration in the caller.
                    fdecl = getattr(caller, "fdecl_by_identifier", {}).get(
                        base.lower())
                    if fdecl is None:
                        fdecl = getattr(
                            caller, "fdecl_by_identifier", {}).get(base)
                    if fdecl is None:
                        continue

                    if not getattr(fdecl, "is_modified", False):
                        fdecl.is_modified = True
                        changed = True


def _infer_user_defined_callable_signatures(conv_info):
    """Infer signatures for user-defined EXTERNAL callable dummy arguments.

    This pass is intentionally conservative:

      - It only infers a signature when the callable is actually invoked inside
        the current procedure (e.g. IF (SELECT(W)) or LOG = SELCTG(...)).
      - It requires the return type and all argument types to be inferred
        consistently. If anything is ambiguous or conflicting, the callable
        remains unresolved and will keep the existing UNHANDLED fallback.

    This function does NOT parse comments.
    """
    fproc = getattr(conv_info, "fproc", None)
    if fproc is None:
        return

    # Cache per procedure. fproc may use __slots__, so we cannot attach new attributes.
    key = id(fproc)
    if key in _INFERRED_CALLABLE_SIGNATURES:
        return

    # Collect dummy arguments declared as user-defined callables (EXTERNAL).
    targets = set()
    for id_tok in (getattr(fproc, "args", None) or []):
        if getattr(id_tok, "value", None) == "*":
            continue
        try:
            fdecl = fproc.get_fdecl(id_tok=id_tok)
        except Exception:
            fdecl = None
        if fdecl is not None and fdecl.is_user_defined_callable():
            targets.add(id_tok.value.lower())

    if not targets:
        _INFERRED_CALLABLE_SIGNATURES[key] = {}
        return

    from fable.tokenization import group_power, extract_identifiers

    def _mpl_type_from_fdecl(fdecl):
        try:
            ctype = convert_data_type(
                conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
        except Exception:
            return None
        return convert_to_mplapack_type(ctype)

    # Best-effort type inference from a C++ expression string.
    num_re = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
    int_re = re.compile(r"^[+-]?\d+$")
    id_head_re = re.compile(r"^([A-Za-z_][A-Za-z0-9_:]*)(?:\b|\[|\.|\()")

    def _mpl_type_from_expr(expr: str):
        if expr is None:
            return None
        s = expr.strip()
        if not s:
            return None
        if s.startswith("&"):
            s = s[1:].strip()
        s = _strip_outer_parens_balanced(s)

        if s in ("true", "false"):
            return "bool"

        if s.startswith("COMPLEX(") or "std::complex" in s:
            return "COMPLEX"

        if int_re.fullmatch(s):
            return "INTEGER"
        if num_re.fullmatch(s) and (("." in s) or ("e" in s) or ("E" in s)):
            return "REAL"

        # Explicit real/imag projections imply REAL.
        if ".real()" in s or ".imag()" in s or "castREAL(" in s:
            return "REAL"

        m = id_head_re.match(s)
        if not m:
            return None
        name = m.group(1).replace(" ", "")
        name = name.split("::")[-1]
        name = name.split(".")[-1]

        fd = fproc.fdecl_by_identifier.get(name.lower())
        if fd is None:
            fd = fproc.fdecl_by_identifier.get(name)
        if fd is None:
            return None
        return _mpl_type_from_fdecl(fd)

    def _find_invocations(tokens):
        """Return [(callee_lower, [arg_expr0, ...]), ...] for calls to targets."""
        found = []
        if tokens is None:
            return found
        prev = None
        for tok in group_power(tokens=tokens):
            if tok.is_seq():
                found.extend(_find_invocations(tok.value))
            elif tok.is_parentheses():
                # Scan nested expressions too.
                found.extend(_find_invocations(tok.value))
                if prev is not None and prev.is_identifier():
                    name = prev.value.lower()
                    if name in targets:
                        arg_str = convert_tokens(
                            conv_info=conv_info,
                            tokens=tok.value,
                            commas=True,
                        )
                        arg_exprs = _split_actuals(arg_str)
                        found.append((name, arg_exprs))
            prev = tok
        return found

    def _update(state, name, arg_types, ret_type):
        st = state.get(name)
        if st is None:
            state[name] = {"ret": ret_type, "args": list(
                arg_types) if arg_types is not None else None, "conflict": False}
            return

        if ret_type is not None:
            if st["ret"] is None:
                st["ret"] = ret_type
            elif st["ret"] != ret_type:
                st["conflict"] = True

        if arg_types is not None:
            if st["args"] is None:
                st["args"] = list(arg_types)
            else:
                if len(st["args"]) != len(arg_types):
                    st["conflict"] = True
                else:
                    for i, (old, new) in enumerate(zip(st["args"], arg_types)):
                        if new is None:
                            continue
                        if old is None:
                            st["args"][i] = new
                        elif old != new:
                            st["conflict"] = True

    state = {}

    execs = getattr(fproc, "executable", None) or []
    for ei in execs:
        # Condition contexts (IF/ELSEIF/DO WHILE/...) imply LOGICAL return.
        cond_tokens = getattr(ei, "cond_tokens", None)
        if cond_tokens is not None:
            for name, arg_exprs in _find_invocations(cond_tokens):
                arg_types = [_mpl_type_from_expr(a) for a in arg_exprs]
                _update(state, name, arg_types, "bool")

        # Assignments can pin the return type to the LHS variable type.
        if getattr(ei, "key", None) == "assignment":
            lhs_ids = extract_identifiers(tokens=ei.lhs_tokens)
            lhs_type = None
            if lhs_ids:
                lhs_fd = fproc.get_fdecl(id_tok=lhs_ids[0])
                if lhs_fd is not None:
                    lhs_type = _mpl_type_from_fdecl(lhs_fd)
            rhs_tokens = getattr(ei, "rhs_tokens", None)
            for name, arg_exprs in _find_invocations(rhs_tokens):
                arg_types = [_mpl_type_from_expr(a) for a in arg_exprs]
                _update(state, name, arg_types, lhs_type)

        # Scan nested invocations inside CALL actual arguments.
        if getattr(ei, "key", None) == "call":
            arg_tok = getattr(ei, "arg_token", None)
            if arg_tok is not None and getattr(arg_tok, "value", None) is not None:
                for name, arg_exprs in _find_invocations(arg_tok.value):
                    arg_types = [_mpl_type_from_expr(a) for a in arg_exprs]
                    _update(state, name, arg_types, None)

    final = {}
    for name, st in state.items():
        if st.get("conflict"):
            continue
        if st.get("ret") is None:
            continue
        args = st.get("args")
        if args is None or any(a is None for a in args):
            continue
        final[name] = (st["ret"], args)
        _INFERRED_CALLABLE_SIGNATURES[key] = final


def convert_to_cpp_function(
        cpp_callback,
        hpp_callback,
        conv_info,
        declaration_only=False,
        force_not_implemented=False):
    # Reset COMPLEX identifier tracking for this procedure
    global complex_identifiers, complex_pointer_identifiers
    complex_identifiers = set()
    complex_pointer_identifiers = set()

    # Ensure that any variable used as a CALL actual argument is treated
    # as "used" at least once, so its declaration is not dropped.
    _mark_call_actuals_used(conv_info)
    _infer_user_defined_callable_signatures(conv_info)

    if (not declaration_only):
        export_save_struct(callback=cpp_callback, conv_info=conv_info)
    fptr = []
    cargs = []

    def cargs_append(ctype, name):
        fptr.append(ctype)
        cargs.append(ctype + " " + name)
    args_fdecl_with_dim = []
    for id_tok in conv_info.fproc.args:
        if (id_tok.value == "*"):
            cargs_append("fem::star_type const&",
                         "/* UNHANDLED: star argument */")
            continue
        assert id_tok.value not in conv_info.vmap
        fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
        conv_info.set_vmap_from_fdecl(fdecl=fdecl)
        assert fdecl.parameter_assignment_tokens is None
        if (fdecl.use_count == 0):
            arg_name = "/* %s */" % prepend_identifier_if_necessary(
                id_tok.value)
        else:
            arg_name = prepend_identifier_if_necessary(id_tok.value)
        if (fdecl.data_type is not None
                and fdecl.data_type.value == "character"):
            # Mark this as a CHARACTER dummy argument so we can
            # skip its local declaration later.
            _mark_dummy_character_arg(conv_info, id_tok.value)
            if (fdecl.dim_tokens is None):
                # Scalar CHARACTER argument.
                #   - Default: classic LAPACK-style (const) char*
                #   - View mode (FABLE_SMALL_CHAR=0): str_cref / str_ref
                if FABLE_SMALL_CHAR_VIEW:
                    if fdecl.is_modified:
                        cargs_append("fem::str_ref", arg_name)
                    else:
                        cargs_append("fem::str_cref", arg_name)
                else:
                    if fdecl.is_modified:
                        cargs_append("char *", arg_name)
                    else:
                        cargs_append("const char *", arg_name)
            else:
                # CHARACTER arrays:
                #   - CHARACTER*1 arrays are modeled as a plain char pointer (char*/const char*)
                #     to match classic LAPACK interfaces (e.g. EI(*) in DLATME).
                #   - Arrays of longer strings keep str_arr_*ref.
                is_len1 = (
                    fdecl.size_tokens is None
                    or (len(fdecl.size_tokens) == 1
                        and fdecl.size_tokens[0].is_integer()
                        and fdecl.size_tokens[0].value == "1")
                )
                if is_len1:
                    if fdecl.is_modified:
                        cargs_append("char *", arg_name)
                    else:
                        cargs_append("const char *", arg_name)
                else:
                    if (len(fdecl.dim_tokens) == 1):
                        cdim = ""
                    else:
                        cdim = "%d" % len(fdecl.dim_tokens)
                    cargs_append(
                        "str_arr_%sref<%s>" % (
                            cconst(fdecl=fdecl, short=True), cdim),
                        arg_name,
                    )
        elif (not fdecl.is_user_defined_callable()):
            # Non-character, non-user-defined argument (typical BLAS/LAPACK args)
            ctype = convert_data_type(
                conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
            mplapack_type = convert_to_mplapack_type(ctype)

            # Determine original Fortran data type code (integer, real, complex, ...)
            dt_code = None
            if fdecl.data_type is not None:
                if isinstance(fdecl.data_type, str):
                    dt_code = fdecl.data_type
                else:
                    dt_code = getattr(fdecl.data_type, "value", None)
                if dt_code is not None:
                    dt_code = dt_code.lower()

            if (fdecl.dim_tokens is None):
                # Scalar argument.
                #   - OUT / INOUT (is_modified=True)      -> <TYPE>&
                #   - IN  scalar (any numeric kind)       -> <TYPE> const   (by value)
                name = prepend_identifier_if_necessary(arg_name)

                if fdecl.is_modified:
                    # OUT / INOUT scalar
                    cargs_append("%s &" % mplapack_type, name)
                else:
                    # IN scalar: pass by value for MPLAPACK/BLAS-style interfaces
                    cargs_append("%s const" % mplapack_type, name)

                # Track COMPLEX scalars
                if dt_code in ("complex", "doublecomplex"):
                    complex_identifiers.add(name)

            else:
                # Array argument: use plain pointer (REAL*, INTEGER*, ...)
                cargs_append("%s *" % mplapack_type, arg_name)

                # Track COMPLEX dummy arrays (COMPLEX* a, x, y, ...)
                if dt_code in ("complex", "doublecomplex") and fdecl.use_count != 0:
                    complex_identifiers.add(arg_name)
                    complex_pointer_identifiers.add(arg_name)

        else:
            passed = conv_info.fproc.externals_passed_by_arg_identifier.get(
                fdecl.id_tok.value)
            if (passed is None or len(passed) == 0):
                ctype = "UNHANDLED"
            sigs = _INFERRED_CALLABLE_SIGNATURES.get(
                id(conv_info.fproc), {}) or {}
            sig = sigs.get(fdecl.id_tok.value.lower())
            if sig is None:
                sig = sigs.get(fdecl.id_tok.value)
            if sig is not None:
                ret_type, arg_types = sig
                args = ", ".join(arg_types)
                # Type-only form used for forward declarations.
                fptr.append(f"{ret_type} (*)({args})")
                # Full parameter form (omit name if it was commented out).
                if str(arg_name).lstrip().startswith("/*"):
                    cargs.append(f"{ret_type} (*)({args}) {arg_name}")
                else:
                    cargs.append(f"{ret_type} (*{arg_name})({args})")
            else:
                passed = conv_info.fproc.externals_passed_by_arg_identifier.get(
                    fdecl.id_tok.value)
                if (passed is None or len(passed) == 0):
                    ctype = "UNHANDLED"
                else:
                    ctype = sorted(passed)[0]
                cargs_append(
                    ctype=ctype + "_function_pointer",
                    name=arg_name)

        if (fdecl.dim_tokens is not None and fdecl.use_count != 0):
            # args_fdecl_with_dim.append(fdecl)
            pass
    cdecl = "void"
    if (conv_info.fproc.name is not None):
        fdecl = conv_info.fproc.get_fdecl(id_tok=conv_info.fproc.name)
        if (fdecl.data_type is not None):
            cdecl = convert_data_type(
                conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
            # Convert return type to MPLAPACK style
            cdecl = convert_to_mplapack_type(cdecl)
            conv_info.vmap[conv_info.fproc.name.value] = "return_value"
    if (declaration_only):
        cpp_callback("")
        cpp_callback("// forward declaration (dependency cycle)")
        if (conv_info.inline_all):
            cpp_callback("inline")
        cpp_callback("%s %s(%s);" % (
            cdecl,
            convert_function_name_to_mplapack(conv_info.fproc.name.value),
            ", ".join(fptr)))
        return
    if (conv_info.fproc.is_passed_as_external):
        if (hpp_callback is None):
            cb = cpp_callback
        else:
            cb = hpp_callback
        cb("")
        cb("typedef %s (*%s_function_pointer)(%s);" % (
            cdecl,
            convert_function_name_to_mplapack(conv_info.fproc.name.value),
            ", ".join(fptr)))
    for callback in [hpp_callback, cpp_callback]:
        if (callback is None):
            continue
        callback("")
        if (callback is cpp_callback):
            produce_leading_comments(callback=callback, fproc=conv_info.fproc)
            produce_fortran_file_comment(
                conv_info=conv_info, callback=callback)
        if (conv_info.inline_all):
            callback("inline")
        callback(cdecl)
        if (callback is hpp_callback):
            last = ";"
        else:
            last = ""
        cname = convert_function_name_to_mplapack(conv_info.fproc.name.value)
        if (len(cargs) == 0):
            callback(cname+"()" + last)
        else:
            callback(cname + "(\n  " + ",\n  ".join(cargs) + ")" + last)
    cpp_callback("{")
    need_local_cmn = (
        getattr(conv_info.fproc, "needs_cmn", False)
        or getattr(conv_info.fproc, "uses_read", False)
        or getattr(conv_info.fproc, "uses_write", False)
    )
    if need_local_cmn:
        cpp_callback("  common cmn;")
    if (cdecl != "void"):
        cpp_callback("  %s %s = %s;" % (
            cdecl,
            conv_info.vmap[conv_info.fproc.name.value],
            zero_shortcut_if_possible(ctype=cdecl)))
    if (force_not_implemented):
        cpp_callback("  throw TBXX_NOT_IMPLEMENTED();")
    else:
        convert_executable(
            callback=cpp_callback,
            conv_info=conv_info,
            args_fdecl_with_dim=args_fdecl_with_dim)
    cpp_callback("}")
    produce_trailing_comments(callback=callback, fproc=conv_info.fproc)


def convert_to_struct(
        callback,
        separate_cmn_hpp,
        fproc,
        struct_type,
        struct_name,
        equivalence_simple,
        id_tok_list):
    assert struct_type in ["common", "save"]
    need_dynamic_parameters = False
    conv_info = conversion_info(fproc=fproc)
    callback("")
    callback("struct %s" % struct_name)
    callback("{")
    sve_equivalences = {}
    cmn_equivalences = {}
    have_variant_block = False
    if (struct_type == "save"
            and conv_info.fproc.conv_hook.needs_is_called_first_time):
        for common_name in sorted(conv_info.fproc.conv_hook.variant_common_names):
            callback("  fem::variant_bindings %s_bindings;" % common_name)
            have_variant_block = True
        #
        cei = conv_info.fproc.classified_equivalence_info()
        sve_equivalences = cei.save.equiv_tok_cluster_by_identifier
        if (len(sve_equivalences) != 0):
            callback("  fem::variant_core_and_bindings save_equivalences;")
            have_variant_block = True
        cmn_equivalences = cei.common.equiv_tok_cluster_by_identifier
    #
    from fable.tokenization import extract_identifiers
    const_identifiers = {}
    const_id_toks = []
    remaining_id_tok_list = []
    for id_tok in id_tok_list:
        if (id_tok.value in sve_equivalences):
            continue
        if (id_tok.value in cmn_equivalences):
            continue
        remaining_id_tok_list.append(id_tok)
        fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
        for tokens in [fdecl.size_tokens, fdecl.dim_tokens]:
            if (tokens is None):
                continue

            def parameter_recursion(tokens):
                have_dynamic_dependency = False
                for id_tok in extract_identifiers(tokens=tokens):
                    if (id_tok.value in const_identifiers):
                        continue
                    const_identifiers[id_tok.value] = None
                    fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
                    hdp = parameter_recursion(
                        tokens=fdecl.required_parameter_assignment_tokens())
                    if (hdp or id_tok.value in conv_info.fproc.dynamic_parameters):
                        have_dynamic_dependency = True
                    const_identifiers[id_tok.value] = have_dynamic_dependency
                    const_id_toks.append(id_tok)
                return have_dynamic_dependency
            parameter_recursion(tokens=tokens)
    initializers = []
    const_definitions = []
    # ISO C++ 9.4.2-4:
    #   "The member shall still be defined in a namespace scope if it is used
    #   in the program and the namespace scope definition shall not contain
    #   an initializer."
    if (len(const_id_toks) != 0):
        if (have_variant_block):
            callback("")
        append_empty_line = False
        for id_tok in const_id_toks:
            fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
            ctype = convert_data_type(
                conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
            if (const_identifiers[id_tok.value]):
                need_dynamic_parameters = True
                callback("  const %s %s;" % (
                    ctype, prepend_identifier_if_necessary(id_tok.value)))
                if (id_tok.value in conv_info.fproc.dynamic_parameters):
                    crhs = "dynamic_params." + prepend_identifier_if_necessary(
                        id_tok.value)
                else:
                    crhs = convert_tokens(
                        conv_info=conv_info, tokens=fdecl.parameter_assignment_tokens)
                initializers.append((id_tok.value, crhs))
            else:
                crhs = convert_tokens(
                    conv_info=conv_info, tokens=fdecl.parameter_assignment_tokens)
                callback("  static const %s %s = %s;" % (
                    ctype, prepend_identifier_if_necessary(id_tok.value), crhs))
                const_definitions.append(
                    "const %s %s::%s;" % (
                        ctype, struct_name, prepend_identifier_if_necessary(id_tok.value)))
                append_empty_line = True
        if (append_empty_line):
            callback("")
    #
    deferred_arr_members = []
    deferred_arr_initializers = []
    for id_tok in remaining_id_tok_list:
        fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
        if (fdecl.id_tok.value in const_identifiers):
            continue
        ctype, cdims, crhs = convert_data_type_and_dims(
            conv_info=conv_info, fdecl=fdecl, crhs=None, force_arr=True)[:3]
        if (cdims is None):
            callback("  %s %s;" % (
                ctype, prepend_identifier_if_necessary(id_tok.value)))
            if (crhs is None):
                crhs = zero_shortcut_if_possible(ctype=ctype)
            initializers.append((id_tok.value, crhs))
        elif (not equivalence_simple or need_dynamic_parameters):
            callback("  %s %s;" % (
                ctype, prepend_identifier_if_necessary(id_tok.value)))
            initializers.append((id_tok.value, "%s, fem::fill0" % cdims))
        else:
            ctype_core = convert_data_type(
                conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
            cstatic_size = convert_dims_to_static_size(
                conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
            callback("  %s %s_memory[%s];" % (
                ctype_core,
                prepend_identifier_if_necessary(id_tok.value),
                cstatic_size))
            if (fdecl.data_type.value == "character"):
                deferred_arr_members.append("  str_arr_ref<%d> %s;" % (
                    len(fdecl.dim_tokens),
                    prepend_identifier_if_necessary(id_tok.value)))
            else:
                deferred_arr_members.append("  %s %s;" % (
                    ad_hoc_change_arr_to_arr_ref(ctype=ctype),
                    prepend_identifier_if_necessary(id_tok.value)))
            deferred_arr_initializers.append((
                id_tok.value,
                "*%s_memory, %s, fem::fill0" % (
                    prepend_identifier_if_necessary(id_tok.value), cdims)))
    if (len(deferred_arr_members) != 0):
        callback("")
        for line in deferred_arr_members:
            callback(line)
        initializers.extend(deferred_arr_initializers)
    n = len(initializers)
    if (n != 0):
        callback("")
        if (not need_dynamic_parameters):
            callback("  %s() :" % struct_name)
        else:
            callback("  %s(" % struct_name)
            callback("    dynamic_parameters const& dynamic_params)")
            callback("  :")
        for i in range(n):
            ii = initializers[i]
            if (i+1 == n):
                comma = ""
            else:
                comma = ","
            callback("    %s(%s)%s" % (
                prepend_identifier_if_necessary(ii[0]), ii[1], comma))
        callback("  {}")
    callback("};")
    #
    if (len(const_definitions) != 0):
        callback("")
        need_ifdef = (separate_cmn_hpp and struct_type == "common")
        if (need_ifdef):
            callback("#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN")
        for cd in const_definitions:
            callback(cd)
        if (need_ifdef):
            callback("#endif")
    #
    return group_args(
        need_dynamic_parameters=need_dynamic_parameters)


def generate_common_report(
        common_fdecl_list_sizes,
        common_equiv_tok_seqs,
        ccode_registry,
        member_registry,
        variant_due_to_equivalence_common_names,
        stringio):
    variant_common_names = set()
    if (stringio is None):
        report = StringIO()
    else:
        report = stringio
    for common_name, fproc_cpp_pairs in ccode_registry.items():
        fprocs_by_cpp = {}
        for fproc, cpp in fproc_cpp_pairs:
            fprocs_by_cpp.setdefault("\n".join(cpp), []).append(fproc)
        if (len(fprocs_by_cpp) != 1):
            variant_common_names.add(common_name)
            fprocs_by_cpp_items = list(fprocs_by_cpp.items())

            def size_key(a):
                return len(a[0])
            fprocs_by_cpp_items.sort(key=size_key, reverse=True)
            import difflib
            diff_function = getattr(difflib, "unified_diff", difflib.ndiff)

            def show_fprocs(label, cpp_fprocs):
                print("procedures %s:" % label,
                      " ".join(sorted([fproc.name.value for fproc in cpp_fprocs[1]])), file=report)
            main_cpp_fprocs = fprocs_by_cpp_items[0]
            print("common name:", common_name, file=report)
            print("number of variants:", len(fprocs_by_cpp_items), file=report)
            print("total number of procedures using the common block:",
                  sum([len(fprocs) for cpp, fprocs in fprocs_by_cpp_items]), file=report)
            show_fprocs("first", main_cpp_fprocs)
            for other_cpp_fprocs in fprocs_by_cpp_items[1:]:
                show_fprocs("second", other_cpp_fprocs)
                print(" ".join([line for line in diff_function(
                              (main_cpp_fprocs[0]+"\n").splitlines(1),
                              (other_cpp_fprocs[0]+"\n").splitlines(1))]), file=report)
    #
    need_empty_line = False
    for identifier in sorted(member_registry.keys()):
        common_names = member_registry[identifier]
        if (len(common_names) != 1):
            print("Name clash: %s in COMMONs: %s" % (
                identifier, ", ".join(sorted(common_names))), file=report)
            need_empty_line = True
    if (need_empty_line):
        print(file=report)
    #
    vv = list(variant_due_to_equivalence_common_names - variant_common_names)
    if (len(vv) != 0):
        print("common variants due to equivalence:", len(vv), file=report)
        size_sums = {}
        for common_name, sizes in common_fdecl_list_sizes.items():
            size_sums[common_name] = sum(sizes)

        vv.sort(key=lambda element: (-size_sums[element], element))
        print("  %-20s   procedures    sum of members" %
              "common name", file=report)
        for common_name in vv:
            print("  %-20s   %8d         %8d" % (
                common_name,
                len(common_fdecl_list_sizes[common_name]),
                size_sums[common_name]), file=report)
        print(file=report)
        print("Locations of equivalence statements:", file=report)
        reported_already = set()
        for common_name in vv:
            print("  %s" % common_name, file=report)
            prev_loc = ""
            tab = []
            max_len_col1 = 6
            for tok_seq in common_equiv_tok_seqs[common_name]:
                sl, i = tok_seq.stmt_location()
                tag = (sl.file_name, sl.line_number, i)
                if (tag in reported_already):
                    break
                reported_already.add(tag)
                vn = tok_seq.value[0].value
                dn, bn = os.path.split(sl.file_name)
                loc = ("%s(%s) %s" % (bn, sl.line_number, dn)).rstrip()
                if (loc == prev_loc):
                    loc = ""
                else:
                    prev_loc = loc
                tab.append((vn, loc))
                max_len_col1 = max(max_len_col1, len(vn))
            if (len(tab) != 0):
                fmt = "    %%-%ds %%s" % max_len_col1
                for row in tab:
                    print(fmt % row, file=report)
    #
    if (len(report.getvalue()) != 0 and stringio is None):
        import sys
        report_file_name = "fable_cout_common_report"
        from libtbx.str_utils import show_string
        print("Writing file:", show_string(report_file_name), file=sys.stderr)
        open(report_file_name, "w").write(report.getvalue())
    #
    return variant_common_names


def convert_commons(
        callback,
        separate_cmn_hpp,
        topological_fprocs,
        dynamic_parameters,
        common_equivalence_simple,
        common_report_stringio):
    if (dynamic_parameters is not None):
        callback("")
        callback("struct dynamic_parameters")
        callback("{")
        for dp_props in dynamic_parameters:
            callback("  %s %s;" % (dp_props.ctype, dp_props.name))
        callback("""
  dynamic_parameters(
    fem::command_line_arguments const& command_line_args)
  :""")
        for dp_props in dynamic_parameters:
            if (dp_props is not dynamic_parameters[-1]):
                c = ","
            else:
                c = ""
            callback("    %s(%s)%s" % (
                prepend_identifier_if_necessary(dp_props.name),
                str(dp_props.default),
                c))
        callback("""\
  {
    fem::dynamic_parameters_from(command_line_args, %d)"""
                 % len(dynamic_parameters))
        for dp_props in dynamic_parameters:
            callback("      .reset_if_given(%s)"
                     % prepend_identifier_if_necessary(dp_props.name))
        callback("    ;")
        callback("  }")
        callback("};")
        callback("")
        callback("typedef")
        callback("  fem::dynamic_parameters_capsule<dynamic_parameters>")
        callback("    dynamic_parameters_capsule;")
    #
    common_fdecl_list_sizes = {}
    common_equiv_tok_seqs = {}
    common_ccode_registry = {}
    member_registry = {}

    variant_common_names = set()
    bottom_up_filtered = []
    for fproc in topological_fprocs.bottom_up_list:
        # Ensure conv_hook exists for all fprocs seen here
        ch = getattr(fproc, "conv_hook", None)
        if ch is None:
            ch = conv_hook_info()
            ch.ignore_common_and_save = False
            fproc.conv_hook = ch
        if not ch.ignore_common_and_save:
            bottom_up_filtered.append(fproc)

    struct_commons_need_dynamic_parameters = set()
    for fproc in bottom_up_filtered:
        fproc.conv_hook.needs_variant_bind = False
        for common_name, common_fdecl_list in fproc.common.items():
            common_fdecl_list_sizes.setdefault(common_name, []).append(
                len(common_fdecl_list))
            id_tok_list = []
            for common_fdecl in common_fdecl_list:
                assert common_fdecl.size_tokens is None
                id_tok_list.append(common_fdecl.id_tok)
                member_registry.setdefault(
                    common_fdecl.id_tok.value, set()).add(common_name)
                if (common_name not in common_equivalence_simple):
                    equiv_tok_cluster = (
                        fproc.equivalence_info()
                        .equiv_tok_cluster_by_identifier.get(common_fdecl.id_tok.value)
                    )
                    if (equiv_tok_cluster is not None):
                        fproc.conv_hook.needs_variant_bind = True
                        variant_common_names.add(common_name)
                        for equiv_tok in equiv_tok_cluster:
                            for tok_seq in equiv_tok.value:
                                common_equiv_tok_seqs.setdefault(common_name, []).append(
                                    tok_seq)
            struct_name = "common_" + common_name
            buffer = []
            info = convert_to_struct(
                callback=buffer.append,
                separate_cmn_hpp=separate_cmn_hpp,
                fproc=fproc,
                struct_type="common",
                struct_name=struct_name,
                equivalence_simple=(common_name in common_equivalence_simple),
                id_tok_list=id_tok_list)
            if (info.need_dynamic_parameters):
                struct_commons_need_dynamic_parameters.add(struct_name)
            common_ccode_registry.setdefault(common_name, []).append(
                (fproc, buffer))
    variant_common_names.update(generate_common_report(
        common_fdecl_list_sizes=common_fdecl_list_sizes,
        common_equiv_tok_seqs=common_equiv_tok_seqs,
        ccode_registry=common_ccode_registry,
        member_registry=member_registry,
        variant_due_to_equivalence_common_names=variant_common_names,
        stringio=common_report_stringio))
    commons_defined_already = set()
    struct_commons = []
    variant_commons = []
    for fproc in bottom_up_filtered:
        fproc.conv_hook.variant_common_names = set()
        for common_name, common_fdecl_list in fproc.common.items():
            if (common_name in variant_common_names):
                fproc.conv_hook.variant_common_names.add(common_name)
                if (common_name not in commons_defined_already):
                    commons_defined_already.add(common_name)
                    variant_commons.append(common_name)
            else:
                if (common_name not in commons_defined_already):
                    commons_defined_already.add(common_name)
                    struct_commons.append("common_"+common_name)
                    for line in common_ccode_registry[common_name][0][1]:
                        callback(line)
    #
    for fproc in bottom_up_filtered:
        if (not fproc.conv_hook.needs_variant_bind):
            fproc.conv_hook.needs_variant_bind = (
                len(fproc.conv_hook.variant_common_names) != 0
                or fproc.classified_equivalence_info().has_save())
        fproc.conv_hook.needs_is_called_first_time = (
            fproc.conv_hook.needs_variant_bind
            or len(fproc.data) != 0)
        fproc.conv_hook.data_init_after_variant_bind = (
            fproc.conv_hook.needs_variant_bind
            and len(fproc.data) != 0)
        if (fproc.conv_hook.needs_is_called_first_time):
            fproc.uses_save = True
    topological_fprocs.each_fproc_update_is_modified()
    topological_fprocs.each_fproc_update_needs_cmn()
    #
    save_struct_buffers = {}
    save_struct_names = []
    for fproc in bottom_up_filtered:
        id_tok_list = []
        for fdecl in fproc.fdecl_by_identifier.values():
            if (fdecl.is_save()):
                id_tok_list.append(fdecl.id_tok)
        if (len(id_tok_list) == 0
                and not fproc.conv_hook.needs_is_called_first_time):
            continue
        id_tok_list.sort(key=lambda token: token.value)
        struct_name = "%s_save" % fproc.name.value
        buffer = []
        info = convert_to_struct(
            callback=buffer.append,
            separate_cmn_hpp=separate_cmn_hpp,
            fproc=fproc,
            struct_type="save",
            struct_name=struct_name,
            equivalence_simple=False,
            id_tok_list=id_tok_list)
        save_struct_buffers[fproc.name.value] = buffer
        if (info.need_dynamic_parameters):
            fproc.conv_hook.needs_sve_dynamic_parameters = True
        save_struct_names.append(struct_name)
    if (len(commons_defined_already) == 0
            and len(save_struct_names) == 0
            and dynamic_parameters is None):
        # Disabled: callback("")
        # Disabled: callback("using fem::common;")
        return
    callback("")
    callback("struct common :")
    leading_bases = ["fem::common"]
    if (dynamic_parameters is not None):
        leading_bases.append("dynamic_parameters_capsule")
    callback("  " + ",\n  ".join(leading_bases + struct_commons))
    callback("{")
    need_empty_line = False
    for common_name in variant_commons:
        callback("  fem::variant_core common_%s;" % common_name)
        need_empty_line = True

    def save_as_sve(struct_name): return struct_name[:-3]+"ve"
    for struct_name in save_struct_names:
        callback("  fem::cmn_sve %s;" % save_as_sve(struct_name))
        need_empty_line = True
    if (need_empty_line):
        callback("")
    initializations = ["fem::common(argc, argv)"]
    if (dynamic_parameters is not None):
        initializations.append(
            "dynamic_parameters_capsule(command_line_args)")
        for struct_name in struct_commons:
            if (struct_name in struct_commons_need_dynamic_parameters):
                initializations.append("%s(dynamic_params)" % struct_name)
    callback("""\
  common(
    int argc,
    char const* argv[])
  :
    %s
  {}""" % ",\n    ".join(initializations))
    callback("};")
    #
    return group_args(
        member_registry=member_registry,
        save_struct_buffers=save_struct_buffers)


include_fem_hpp = ""


def include_guard(callback, namespace, suffix):
    if namespace:
        s = namespace.upper().replace("::", "_") + suffix
    else:
        s = "GUARD" + suffix
    callback("#ifndef %s" % s)
    callback("#define %s" % s)
    callback("")


def open_namespace(callback, namespace, using_namespace_major_types=True):
    # Disabled: no namespace wrapping
    # ns = namespace.split("::")
    # for component in ns:
    #     callback("namespace %s {" % component)
    # if (using_namespace_major_types):
    #     callback("""
    # using namespace fem::major_types;""")
    if namespace:
        ns = namespace.split("::")
    else:
        ns = []
    return ns


def close_namespace(callback, namespace, hpp_guard):
    # Disabled: no namespace wrapping
    # callback("")
    # ns = namespace.split("::")
    # callback("%s // namespace %s" % ("}"*len(ns), namespace))
    if (hpp_guard):
        callback("")
        callback("#endif // GUARD")
    if namespace:
        ns = namespace.split("::")
    else:
        ns = []
    return ns


class hpp_cpp_buffers(object):

    __slots__ = ["hpp", "cpp"]

    def __init__(O):
        O.hpp = []
        O.cpp = []


def convert_program(callback, global_conv_info, namespace, hpp_guard, debug):
    main_calls = []
    for fproc in global_conv_info.topological_fprocs.bottom_up_list:
        if (not fproc.is_program()):
            continue
        conv_info = global_conv_info.specialized(fproc=fproc)
        export_save_struct(callback=callback, conv_info=conv_info)
        cname = fproc.name.value
        main_calls.append(cname)
        callback("")
        produce_leading_comments(callback=callback, fproc=fproc)
        produce_fortran_file_comment(conv_info=conv_info, callback=callback)
        callback("""\
void
%s(
  int argc,
  char const* argv[])
{""" % cname)
        if (not fproc.needs_cmn):
            callback("""\
  if (argc != 1) {
    throw std::runtime_error("Unexpected command-line arguments.");
  }""")
        result_buffer = []
        try:
            convert_executable(
                callback=result_buffer.append,
                conv_info=conv_info,
                blockdata=global_conv_info.topological_fprocs.all_fprocs.blockdata)
        except Exception:
            if (not debug):
                raise
            show_traceback()
        else:
            need_cmn_obj = (
                fproc.needs_cmn
                or getattr(fproc, "uses_read", False)
                or getattr(fproc, "uses_write", False)
                or bool(getattr(global_conv_info.topological_fprocs.all_fprocs, "blockdata", None))
            )
            if need_cmn_obj:
                callback("  common cmn(argc, argv);")
            for line in result_buffer:
                callback(line)
            callback("}")
        produce_trailing_comments(callback=callback, fproc=fproc)
    #
    ns = close_namespace(
        callback=callback, namespace=namespace, hpp_guard=hpp_guard)
    #
    if (len(main_calls) != 0):
        callback("")
        callback("""\
int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    %s);
}""" % "::".join(ns + [main_calls[0]]))


def get_missing_external_return_type(fdecls):
    for fdecl in fdecls:
        if (fdecl.data_type is not None):
            return convert_data_type(
                conv_info=conversion_info(), fdecl=fdecl, crhs=None)[0]
    return "void"


default_arr_nd_size_max = 256


def _postprocess_mplapack_labels_and_comments(lines):
    """MPLAPACK-specific postprocessing for labels, comments, and trivial zero offsets."""
    import re

    new_lines = []
    for line in lines:
        # 1) Fix Mxerbla("XXXX ", info) labels using the MPLAPACK name map
        m = re.search(r'Mxerbla\("([^"]+)"', line)
        if m:
            label = m.group(1)      # e.g. "ZTRSV "
            core = label.strip()    # "ZTRSV"
            lower = core.lower()    # "ztrsv"
            mapped = _MPLAPACK_NAME_MAP.get(lower)
            if mapped is None:
                mapped = _mplapack_default_name(core)
            # preserve trailing spaces
            suffix = label[len(core):]
            repl = f'Mxerbla("{mapped}{suffix}"'
            line = line[:m.start()] + repl + line[m.end():]

        # 2) Fix "//     End of XXXX" comments using the MPLAPACK name map
        m = re.search(r'(End of )([A-Za-z][A-Za-z0-9_]*)', line)
        if m:
            core = m.group(2)
            lower = core.lower()
            mapped = _MPLAPACK_NAME_MAP.get(lower)
            if mapped is None:
                mapped = _mplapack_default_name(core)
            line = line[:m.start(2)] + mapped + line[m.end(2):]

        # 3) Remove trivial zero row offset: (1 - 1) + ...
        #    Example: a[(1 - 1) + (j - 1) * lda] -> a[(j - 1) * lda]
        line = re.sub(r'\(\s*1\s*-\s*1\s*\)\s*\+\s*', '', line)

        # 4) Remove trivial zero column offset: + (1 - 1) * lda
        #    Example: a[(i - 1) + (1 - 1) * lda] -> a[(i - 1)]
        line = re.sub(
            r'\+\s*\(\s*1\s*-\s*1\s*\)\s*\*\s*[A-Za-z_][A-Za-z0-9_]*',
            '',
            line,
        )
        #    Example: a[(1 - 1) * lda + (i - 1)] -> a[(i - 1)]
        line = re.sub(
            r'\(\s*1\s*-\s*1\s*\)\s*\*\s*[A-Za-z_][A-Za-z0-9_]*\s*\+\s*',
            '',
            line,
        )

        # NOTE: Do not touch general index expressions here.
        # 1-based to 0-based conversion must be handled in the main converter logic.

        new_lines.append(line)
    return new_lines


def _normalize_fortran_comment_prefix(lines):
    """Normalize Fortran-derived comments: //C... -> // ..."""
    normalized = []
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("//C"):
            # Leading whitespace before //C
            leading = line[:len(line) - len(stripped)]
            # Drop the 'C'
            rest = stripped[3:]  # after "//C"
            # Ensure there is exactly one space after //
            # rest may be "" or like "     foo"
            rest = rest.lstrip()
            if rest:
                new_line = f"{leading}// {rest}"
            else:
                new_line = f"{leading}//"
            normalized.append(new_line)
        else:
            normalized.append(line)
    return normalized


def _postprocess_complex_initializers(lines):
    """Normalize COMPLEX(a, b) initializers.

    Rewrite:
      COMPLEX z = (a, b);        -> COMPLEX z = COMPLEX(a, b);
      const COMPLEX z = (a, b);  -> const COMPLEX z = COMPLEX(a, b);
    """
    pat_var = re.compile(
        r'^(\s*COMPLEX\s+[A-Za-z_][A-Za-z0-9_]*\s*=\s*)\(\s*([^,]+?)\s*,\s*([^,]+?)\s*\);'
    )
    pat_const = re.compile(
        r'^(\s*const\s+COMPLEX\s+[A-Za-z_][A-Za-z0-9_]*\s*=\s*)\(\s*([^,]+?)\s*,\s*([^,]+?)\s*\);'
    )
    out = []
    for line in lines:
        line = pat_const.sub(r'\1COMPLEX(\2, \3);', line)
        line = pat_var.sub(r'\1COMPLEX(\2, \3);', line)
        out.append(line)
    return out


def _postprocess_complex_constant_assignments(lines):
    """Rewrite assignments of real literal pairs (a, b) into COMPLEX(a, b).

    Example:
      return_value = (0.0, 0.0);    ->  return_value = COMPLEX(0.0, 0.0);
      return_value = (1.0, 2.377);  ->  return_value = COMPLEX(1.0, 2.377);
    """
    # Real literal: optional sign, digits with optional decimal point, optional exponent
    real_lit = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?'
    pat = re.compile(
        rf'(\=\s*)\(\s*({real_lit})\s*,\s*({real_lit})\s*\);'
    )
    out = []
    for line in lines:
        line = pat.sub(r'\1COMPLEX(\2, \3);', line)
        out.append(line)
    return out


def _postprocess_intrinsic_aliases(lines):
    """Final cleanup for intrinsic helper names.

    Ensure that fem::abs, fem::dabs, fem::cdabs, fem::pow2
    are printed without the fem:: namespace.
    """
    out = []
    for line in lines:
        # abs family
        line = line.replace("fem::cdabs", "abs")
        line = line.replace("fem::dabs", "abs")
        line = line.replace("fem::abs", "abs")

        # pow2 helper
        line = line.replace("fem::pow2", "pow2")

        out.append(line)
    return out


def _postprocess_comment_name_map(lines):
    """Apply MPLAPACK name mapping inside C++ comment lines.

    For each line that starts with '//' (after optional whitespace),
    replace occurrences of uppercased Fortran routine names with
    their mapped C++ names from _MPLAPACK_NAME_MAP.
    """
    out = []
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("//"):
            for f77_name, cpp_name in _MPLAPACK_NAME_MAP.items():
                u = f77_name.upper()
                if u in line:
                    pattern = r"\b" + re.escape(u) + r"\b"
                    line = re.sub(pattern, cpp_name, line)
        out.append(line)
    return out


def _postprocess_math_intrinsics_upper(lines):
    """Rewrite math intrinsic calls (atan2, cos, sin, tan, log, exp, max, min, abs) to uppercase names."""
    out = []
    for line in lines:
        # ----------------------------------------------------------
        # 1) Canonicalize fem:: and double-precision variants
        #    to plain lowercase intrinsics
        # ----------------------------------------------------------

        # log / log10 (fem::LOG, fem::DLOG, LOG, DLOG, dlog, dlog10, etc.)
        line = re.sub(r'\bfem::d?log10\s*\(',
                      'log10(', line, flags=re.IGNORECASE)
        line = re.sub(r'\bfem::d?log\s*\(',
                      'log(',   line, flags=re.IGNORECASE)
        line = re.sub(r'\bdlog10\s*\(',
                      'log10(', line, flags=re.IGNORECASE)
        line = re.sub(r'\bdlog\s*\(',
                      'log(',   line, flags=re.IGNORECASE)

        # exp / pow / mod (fem::EXP, fem::exp)
        line = re.sub(r'\bfem::exp\s*\(',
                      'exp(',   line, flags=re.IGNORECASE)
        line = re.sub(r'\bfem::pow\s*\(',
                      'pow(',   line, flags=re.IGNORECASE)
        line = re.sub(r'\bfem::mod\s*\(',
                      'mod(',   line, flags=re.IGNORECASE)

        # NINT-like intrinsics (fem::nint, fem::idnint)
        line = re.sub(r'\bfem::nint\s*\(',
                      'nint(',  line, flags=re.IGNORECASE)
        line = re.sub(r'\bfem::idnint\s*\(',
                      'nint(',  line, flags=re.IGNORECASE)

        # ----------------------------------------------------------
        # 2) Trigonometric intrinsics
        # ----------------------------------------------------------
        line = re.sub(r'\bstd::atan2\s*\(', 'atan2(', line)
        line = re.sub(r'\batan2\s*\(',      'atan2(', line)
        line = re.sub(r'\bstd::cos\s*\(',   'cos(',   line)
        line = re.sub(r'\bcos\s*\(',        'cos(',   line)
        line = re.sub(r'\bstd::sin\s*\(',   'sin(',   line)
        line = re.sub(r'\bsin\s*\(',        'sin(',   line)
        line = re.sub(r'\bstd::tan\s*\(',   'tan(',   line)
        line = re.sub(r'\btan\s*\(',        'tan(',   line)

        # ----------------------------------------------------------
        # 3) Logarithm and exponential intrinsics (std:: + plain)
        # ----------------------------------------------------------
        line = re.sub(r'\bstd::log10\s*\(', 'log10(', line)
        line = re.sub(r'\blog10\s*\(',      'log10(', line)
        line = re.sub(r'\bstd::log\s*\(',   'log(',   line)
        line = re.sub(r'\blog\s*\(',        'log(',   line)
        line = re.sub(r'\bstd::exp\s*\(',   'exp(',   line)
        line = re.sub(r'\bexp\s*\(',        'exp(',   line)

        # ----------------------------------------------------------
        # 4) Extremum intrinsics
        # ----------------------------------------------------------
        line = re.sub(r'\bstd::max\s*\(',   'max(',   line)
        line = re.sub(r'\bmax\s*\(',        'max(',   line)
        line = re.sub(r'\bstd::min\s*\(',   'min(',   line)
        line = re.sub(r'\bmin\s*\(',        'min(',   line)

        # ----------------------------------------------------------
        # 5) Absolute value intrinsics
        # ----------------------------------------------------------
        line = re.sub(r'\bstd::abs\s*\(',   'abs(',   line)
        line = re.sub(r'\babs\s*\(',        'abs(',   line)

        # ----------------------------------------------------------
        # 6) Other intrinsics
        # ----------------------------------------------------------
        line = re.sub(r'\bstd::mod\s*\(',   'mod(',   line)
        line = re.sub(r'\bmod\s*\(',        'mod(',   line)
        line = re.sub(r'\bstd::sqrt\s*\(',  'sqrt(',  line)
        line = re.sub(r'\bsqrt\s*\(',       'sqrt(',  line)
        line = re.sub(r'\bpow2\s*\(',       'pow2(',  line)
        line = re.sub(r'\bstd::pow\s*\(',   'pow(',   line)
        line = re.sub(r'\bpow\s*\(',        'pow(',   line)

        # ----------------------------------------------------------
        # 7) NINT-like rounding intrinsics (already canonicalized)
        # ----------------------------------------------------------
        line = re.sub(r'\bnint\s*\(',       'nint(',  line)

        # ----------------------------------------------------------
        # 8) conjg
        # ----------------------------------------------------------
        line = re.sub(r'\bfem::dconjg\s*\(', 'conj(',   line)

        out.append(line)
    return out


def _postprocess_minmax_parens(lines):
    """Drop one level of redundant parentheses around MIN/MAX in index shifts.

    Example:
      ((MIN(k1 + 1, m)) - 1) -> (MIN(k1 + 1, m) - 1)
      ((MAX(i, j)) - 1)      -> (MAX(i, j) - 1)
    """
    pat_min = re.compile(r'\(\((min\([^()]*\))\)\s*-\s*1\)')
    pat_max = re.compile(r'\(\((max\([^()]*\))\)\s*-\s*1\)')
    out = []
    for line in lines:
        line = pat_min.sub(r'(\1 - 1)', line)
        line = pat_max.sub(r'(\1 - 1)', line)
        out.append(line)
    return out


def _postprocess_ilaenv_name_map(lines):
    """Apply MPLAPACK name mapping inside iMlaenv calls.

    Example:
      iMlaenv(3, "DGERQF", " ", m, n, -1, -1)
        -> iMlaenv(3, "Rgerqf", " ", m, n, -1, -1)
    """
    # Match: iMlaenv(..., "NAME", ...)
    pat = re.compile(
        r'(\biMlaenv\s*\(\s*[^,]*,\s*")([A-Za-z0-9_]+)("\s*,)'
    )

    def replace_name(name: str) -> str:
        """Map Fortran routine name to C++ name using _MPLAPACK_NAME_MAP."""
        for f77_name, cpp_name in _MPLAPACK_NAME_MAP.items():
            if name.upper() == f77_name.upper():
                return cpp_name
        return name

    out = []
    for line in lines:
        if "iMlaenv" in line:
            def _repl(m):
                old = m.group(2)
                new = replace_name(old)
                return m.group(1) + new + m.group(3)
            line = pat.sub(_repl, line)
        out.append(line)
    return out


def _postprocess_ilaenv_char_concat(lines):
    """Convert character concatenation in iMlaenv's 3rd argument to CHAR2/CHAR3.

    FORTRAN: ILAENV(1, 'DORMQR', SIDE // TRANS, M-1, N, M-1, -1)
    After FABLE: iMlaenv(1, "Rormqr", side + trans, m - 1, n, m - 1, -1)
    Expected:    iMlaenv(1, "Rormqr", CHAR2(side, trans), m - 1, n, m - 1, -1)

    Also handles 3-way concatenation:
    FORTRAN: SIDE // TRANS // DIAG
    After FABLE: side + trans + diag
    Expected:    CHAR3(side, trans, diag)

    Requires header definition:
        #define CHAR2(a, b) ((char[]){(a)[0], (b)[0], '\\0'})
        #define CHAR3(a, b, c) ((char[]){(a)[0], (b)[0], (c)[0], '\\0'})
    """
    import re

    def find_matching_paren(s, start):
        """Find index of closing paren matching the one at 'start'."""
        depth = 0
        for i in range(start, len(s)):
            if s[i] == '(':
                depth += 1
            elif s[i] == ')':
                depth -= 1
                if depth == 0:
                    return i
        return -1

    def split_args(s):
        """Split argument string by top-level commas."""
        parts = []
        current = []
        depth = 0
        for ch in s:
            if ch == '(':
                depth += 1
                current.append(ch)
            elif ch == ')':
                depth -= 1
                current.append(ch)
            elif ch == ',' and depth == 0:
                parts.append(''.join(current).strip())
                current = []
            else:
                current.append(ch)
        if current:
            parts.append(''.join(current).strip())
        return parts

    def is_char_concat(arg):
        """Check if argument looks like 'var1 + var2' or 'var1 + var2 + var3'.

        We support simple concatenations of (possibly substringed) CHARACTER
        variables, e.g.
            JOB(1:1)   // COMPZ(1:1)
            JOB(:1)    // COMPZ(:1)
            JOBU       // JOBVT
        and normalize them to base identifiers:
            ['job', 'compz'], ['jobu', 'jobvt'], etc.
        """
        stripped = arg.strip()

        # Split by '+' at top level (no nesting across parentheses)
        parts = []
        current = []
        depth = 0
        for ch in stripped:
            if ch == '(':
                depth += 1
                current.append(ch)
            elif ch == ')':
                depth -= 1
                current.append(ch)
            elif ch == '+' and depth == 0:
                parts.append(''.join(current).strip())
                current = []
            else:
                current.append(ch)
        if current:
            parts.append(''.join(current).strip())

        # We only care about 2 or 3-term concatenations.
        if len(parts) < 2 or len(parts) > 3:
            return None

        # Normalize each part:
        #   JOB(1,1)  -> "job"
        #   JOB(:1)   -> "job"
        #   JOB(1:1)  -> "job"
        #   JOB       -> "job"
        ident_pat = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')
        substr_pat1 = re.compile(
            r'^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*1\s*,\s*1\s*\)$')
        substr_pat2 = re.compile(
            r'^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*:?\s*1\s*\)$')
        substr_pat3 = re.compile(
            r'^([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*1\s*:\s*1\s*\)$')
        # Also accept C++-style element forms produced by earlier rewrites:
        #   job[0], &job[0], job[(1) - 1], &job[(1) - 1]
        arr0_pat = re.compile(
            r'^\s*&?\s*([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*0\s*\]\s*$')
        arr1m1_pat = re.compile(
            r'^\s*&?\s*([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*\(?\s*1\s*\)?\s*-\s*1\s*\]\s*$')

        norm_parts = []
        for p in parts:
            p0 = p.strip()
            m = (substr_pat1.match(p0) or substr_pat2.match(p0) or substr_pat3.match(p0)
                 or arr0_pat.match(p0) or arr1m1_pat.match(p0))
            base = m.group(1) if m else p0
            if not ident_pat.match(base):
                return None
            norm_parts.append(base.lower())

        return norm_parts

    def rewrite_line(line):
        # Find iMlaenv( calls
        pattern = r'\biMlaenv\s*\('
        result = []
        pos = 0

        for m in re.finditer(pattern, line):
            result.append(line[pos:m.start()])
            call_start = m.start()
            paren_start = m.end() - 1

            paren_end = find_matching_paren(line, paren_start)
            if paren_end < 0:
                result.append(line[call_start:m.end()])
                pos = m.end()
                continue

            # Extract arguments
            args_str = line[paren_start + 1:paren_end]
            args = split_args(args_str)

            # iMlaenv has at least 7 arguments, 3rd is the character string
            if len(args) >= 3:
                third_arg = args[2]
                char_parts = is_char_concat(third_arg)
                if char_parts:
                    # Replace with CHAR2 or CHAR3
                    if len(char_parts) == 2:
                        args[2] = f"CHAR2({char_parts[0]}, {char_parts[1]})"
                    elif len(char_parts) == 3:
                        args[2] = f"CHAR3({char_parts[0]}, {char_parts[1]}, {char_parts[2]})"
                    new_call = f"iMlaenv({', '.join(args)})"
                    result.append(new_call)
                    pos = paren_end + 1
                    continue

            # No transformation needed
            result.append(line[call_start:paren_end + 1])
            pos = paren_end + 1

        result.append(line[pos:])
        return ''.join(result)

    if isinstance(lines, list):
        return [rewrite_line(line) for line in lines]
    elif isinstance(lines, str):
        return '\n'.join(rewrite_line(l) for l in lines.split('\n'))
    return lines


def _postprocess_strip_float_suffix(lines):
    """Remove 'f' suffix from floating literals like 1.0f, 0.5f, 3.14e-1f."""
    # Match patterns like:
    #   1.0f, 0.5f, 3.14e-1f, 1.F, 1.e+0F, etc.
    pat = re.compile(r'(\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)[fF]\b')
    out = []
    for line in lines:
        out.append(pat.sub(r'\1', line))
    return out


def _postprocess_strip_wp_kind_suffix(lines):
    """Remove leftover Fortran kind suffixes like '_wp' from numeric literals.

    The tokenizer can split `0.0_wp` into `0.0` + `_wp`, and the C++ printer
    may additionally append a float suffix, producing `0.0f_wp`. These are
    not valid C++ literals. Normalize them to plain numeric literals.

    Examples:
      0.0f_wp  -> 0.0
      1.0_wp   -> 1.0
      1_wp     -> 1
      1.e-3_wp -> 1.e-3
    """
    pat = re.compile(
        r'(?<![A-Za-z0-9_])'
        r'(?P<num>(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+\-]?\d+)?)'
        r'(?:[fF])?'
        r'_wp\b',
        flags=re.IGNORECASE,
    )
    out = []
    for line in lines:
        idx = line.find("//")
        if idx >= 0:
            code, comment = line[:idx], line[idx:]
        else:
            code, comment = line, ""
        code = pat.sub(r'\g<num>', code)
        out.append(code + comment)
    return out


def _postprocess_index_zero_simplify(text):
    """
    Postprocess index expressions to remove redundant arithmetic in array
    subscripts, such as (1 - 1), +0, *0, and double parentheses around
    MIN/MAX(...)-1.

    This function is deliberately conservative:
    - It MUST NOT rewrite `name[0]` into `name`.
      That would change scalar arguments (REAL&) into pointers (REAL*),
      which breaks many LAPACK calls.
    - It accepts either a single string or a list of strings.
      In the latter case it returns a list of processed lines.
    - It does not touch comments starting from '//' on each line.
    """
    import re

    def _simplify_zero_mul_parenthesized(e: str) -> str:
        """Replace '0 * ( ... )' with '0' for a balanced parenthesis group."""
        pat = re.compile(r"\b0\s*\*\s*\(")
        while True:
            m = pat.search(e)
            if not m:
                break
            paren_start = m.end() - 1  # points to '('
            depth = 0
            i = paren_start
            while i < len(e):
                ch = e[i]
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        break
                i += 1
            if depth != 0:
                # Unbalanced parentheses: stop trying to simplify.
                break
            e = e[:m.start()] + "0" + e[i + 1:]
        return e

    # Simplify only the code part (before //), leave comments untouched.
    def simplify_code(code: str) -> str:
        # 0) Remove redundant parentheses for an array element immediately after '['.
        #    Example:
        #       sva[(iwork[p - 1]) - 1]  ->  sva[iwork[p - 1] - 1]
        code = re.sub(
            r"\[\s*\(\s*([A-Za-z_][A-Za-z0-9_]*\s*\[[^\[\]]*\])\s*\)",
            r"[\1",
            code,
        )

        # 1) Simplify expressions inside [...]
        bracket_re = re.compile(r"\[(.*?)\]")

        def simplify_expr(expr: str) -> str:
            e = expr

            # (i) -> i, (0) -> 0 when the whole index is just a single
            # parenthesized identifier or integer literal.
            m_simple = re.fullmatch(
                r"\(\s*([A-Za-z_][A-Za-z0-9_]*|\d+)\s*\)", e)
            if m_simple:
                e = m_simple.group(1)

            # (1 - 1) -> 0
            e = re.sub(r"\(\s*1\s*-\s*1\s*\)", "0", e)

            # 0 * ( ... ) -> 0  (balanced parentheses)
            e = _simplify_zero_mul_parenthesized(e)

            # 0 * ATOM -> 0   (ATOM: identifier or integer literal)
            e = re.sub(r"\b0\s*\*\s*(?:[A-Za-z_][A-Za-z0-9_]*|\d+)\b", "0", e)
            # ATOM * 0 -> 0
            e = re.sub(r"\b(?:[A-Za-z_][A-Za-z0-9_]*|\d+)\s*\*\s*0\b", "0", e)

            # Remove leading "0 + ..." safely (only at start).
            e = re.sub(r"^\s*0\s*\+\s*", "", e)

            # Remove "+ 0" and "- 0" only when the 0 is a standalone term.
            # Do NOT remove when it is followed by '*' or '/' (e.g., '+ 0 * ldb').
            e = re.sub(r"\+\s*0\b(?!\s*[*\/])", "", e)
            e = re.sub(r"\-\s*0\b(?!\s*[*\/])", "", e)

            # Collapse whitespace inside the index expression
            e = re.sub(r"\s+", " ", e).strip()
            return e

        def repl_brackets(m: re.Match) -> str:
            inner = m.group(1)
            return "[" + simplify_expr(inner) + "]"

        code2 = bracket_re.sub(repl_brackets, code)

        # 2) Simplify double-parenthesized MIN/MAX(...)-1 in the code part
        pattern_minmax = re.compile(
            r"\(\s*\(\s*("
            r"(?:[Mm][Ii][Nn]|[Mm][Aa][Xx])"   # MIN or MAX, case-insensitive
            r"\("
            # allow one level of nested parentheses
            r"(?:[^()]*|\([^()]*\))*"
            r"\)"
            r")\s*\)\s*-\s*1\s*\)",
            re.DOTALL,
        )
        code2 = pattern_minmax.sub(r"(\1 - 1)", code2)

        return code2

    def simplify_line(line: str) -> str:
        idx = line.find("//")
        if idx < 0:
            return simplify_code(line)
        code = line[:idx]
        comment = line[idx:]
        return simplify_code(code) + comment

    # If we get a list of lines, process each line and return a list again.
    if isinstance(text, list):
        return [simplify_line(line) for line in text]

    # Normal case: a single string containing the whole file.
    if isinstance(text, str):
        lines = text.split("\n")
        new_lines = [simplify_line(line) for line in lines]
        return "\n".join(new_lines)

    # Fallback: unknown type, do not try to be clever.
    return text


def _postprocess_array_slice_intrinsics(lines):
    """Rewrite array slice intrinsics (Mmaxloc, Mmaxval, Mminval), including stride.

    Mmaxloc: Mmaxloc(array[__SLICE__(start, end)], dim) -> Mmaxloc(array, start, end, dim)
    Mmaxval: Mmaxval(array[__SLICE__(start, end)])          -> Mmaxval(array, start, end, 1)
    Mminval: Mminval(array[__SLICE__(start, end)])          -> Mminval(array, start, end, 1)
    Mmaxval: Mmaxval(array[__SLICE__(start, end, step)])    -> Mmaxval(array, start, end, step)
    Mminval: Mminval(array[__SLICE__(start, end, step)])    -> Mminval(array, start, end, step)
    Mmaxloc: Mmaxloc(array[__SLICE__(start, end, step)], d) -> Mmaxloc(array, start, end, step, d)

    FORTRAN examples:
        itemp = maxloc( work( (n+j):(2*n) ), 1 )  -> Mmaxloc(work, (n+j), (2*n), 1)
        emin = abs(maxval(s(isbeg:isbeg+nsl-1)))  -> abs(Mmaxval(s, isbeg, isbeg+nsl-1, 1))
    """
    import re

    def find_matching_paren(s, start):
        """Find the index of the closing paren matching the one at 'start'."""
        depth = 0
        for i in range(start, len(s)):
            if s[i] == '(':
                depth += 1
            elif s[i] == ')':
                depth -= 1
                if depth == 0:
                    return i
        return -1

    def split_args(s):
        """Split by top-level comma."""
        parts = []
        current = []
        depth = 0
        for ch in s:
            if ch == '(':
                depth += 1
                current.append(ch)
            elif ch == ')':
                depth -= 1
                current.append(ch)
            elif ch == ',' and depth == 0:
                part = ''.join(current).strip()
                if part:
                    parts.append(part)
                current = []
            else:
                current.append(ch)
        if current:
            part = ''.join(current).strip()
            if part:
                parts.append(part)
        return parts

    # Pattern to match: FuncName(identifier[__SLICE__(args)]...)
    # FuncName is one of: (Mmaxloc|maxloc), (Mmaxval|maxval), (Mminval|minval)
    pattern = re.compile(
        r'\b(Mmaxloc|maxloc|Mmaxval|maxval|Mminval|minval)\s*\(\s*'
        r'([A-Za-z_][A-Za-z0-9_]*)'
        r'\s*\[\s*__SLICE__\s*\('
    )

    def rewrite_line(line):
        result = []
        pos = 0

        for m in pattern.finditer(line):
            result.append(line[pos:m.start()])
            func_name = m.group(1)
            array_name = m.group(2)

            # Normalize function name to MPLAPACK helper name.
            base = func_name.lower()
            if base in ("mmaxloc", "maxloc"):
                out_name = "Mmaxloc"
            elif base in ("mmaxval", "maxval"):
                out_name = "Mmaxval"
            elif base in ("mminval", "minval"):
                out_name = "Mminval"
            else:
                out_name = func_name

            # Find the end of __SLICE__(...)
            slice_paren_start = m.end() - 1  # position of '(' after __SLICE__
            slice_paren_end = find_matching_paren(line, slice_paren_start)
            if slice_paren_end < 0:
                # Fallback: keep original text
                result.append(line[m.start():m.end()])
                pos = m.end()
                continue

            # Extract slice args: (start, end) or (start, end, step)
            slice_args_str = line[slice_paren_start + 1:slice_paren_end]
            slice_arg_list = split_args(slice_args_str)

            if len(slice_arg_list) not in (2, 3):
                result.append(line[m.start():m.end()])
                pos = m.end()
                continue

            start_arg = slice_arg_list[0].strip()
            end_arg = slice_arg_list[1].strip()
            incx_arg = "1"
            if len(slice_arg_list) == 3:
                incx_arg = slice_arg_list[2].strip()

            # After __SLICE__(...) should be "]" followed by optional ", dim" then ")"
            rest_start = slice_paren_end + 1
            rest_pos = rest_start
            while rest_pos < len(line) and line[rest_pos] in ' \t':
                rest_pos += 1
            if rest_pos >= len(line) or line[rest_pos] != ']':
                result.append(line[m.start():m.end()])
                pos = m.end()
                continue
            rest_pos += 1  # skip ']'

            # Skip whitespace
            while rest_pos < len(line) and line[rest_pos] in ' \t':
                rest_pos += 1

            if rest_pos >= len(line):
                result.append(line[m.start():m.end()])
                pos = m.end()
                continue

            # Check what comes next: ',' (has extra args) or ')' (no extra args)
            if line[rest_pos] == ')':
                # No extra arguments (Mmaxval, Mminval style)
                if out_name in ("Mmaxval", "Mminval"):
                    new_call = f"{out_name}({array_name}, {start_arg}, {end_arg}, {incx_arg})"
                else:
                    # Fallback (shouldn't happen for the targeted intrinsics)
                    new_call = f"{out_name}({array_name}, {start_arg}, {end_arg})"
                result.append(new_call)
                pos = rest_pos + 1
            elif line[rest_pos] == ',':
                # Has extra arguments (Mmaxloc style with dim)
                rest_pos += 1  # skip ','

                # Find the closing ')' of the function call
                paren_depth = 1
                extra_args_start = rest_pos
                close_paren_pos = -1
                for i in range(rest_pos, len(line)):
                    if line[i] == '(':
                        paren_depth += 1
                    elif line[i] == ')':
                        paren_depth -= 1
                        if paren_depth == 0:
                            close_paren_pos = i
                            break

                if close_paren_pos < 0:
                    result.append(line[m.start():m.end()])
                    pos = m.end()
                    continue

                extra_args = line[extra_args_start:close_paren_pos].strip()
                if out_name == "Mmaxloc" and len(slice_arg_list) == 3:
                    # Strided slice: keep dim as the last argument.
                    #   maxloc(a(i:j:k), dim) -> Mmaxloc(a, i, j, k, dim)
                    new_call = f"{out_name}({array_name}, {start_arg}, {end_arg}, {incx_arg}, {extra_args})"
                else:
                    # Backward-compatible: keep the existing 4-arg form.
                    new_call = f"{out_name}({array_name}, {start_arg}, {end_arg}, {extra_args})"
                result.append(new_call)
                pos = close_paren_pos + 1
            else:
                result.append(line[m.start():m.end()])
                pos = m.end()
                continue

        result.append(line[pos:])
        return ''.join(result)

    if isinstance(lines, list):
        return [rewrite_line(line) for line in lines]
    elif isinstance(lines, str):
        return '\n'.join(rewrite_line(l) for l in lines.split('\n'))
    return lines


# Keep old name as alias for compatibility
_postprocess_mmaxloc = _postprocess_array_slice_intrinsics


def _postprocess_slice_assignment(lines):
    """Convert array slice assignment to for loop.

    Supported patterns:

      1D slice with __SLICE__ macro:
        A[__SLICE__(i_start, i_end)] op= expr;
        A[__SLICE__(i_start, i_end, i_step)] op= expr;
        -> step=1:  for (INTEGER i = i_start; i <= i_end; i++) { A[i - 1] op expr; }
        -> step!=1: for (INTEGER i = i_start; ...; i += i_step) { A[i - 1] op expr; }

      2D slice with __SLICE2D__ macro (first dimension only):
        A[__SLICE2D__(i_start, i_end, j, ldA)] op= expr;
        -> for (INTEGER i = i_start; i <= i_end; i++) {
               A[(i - 1) + (j - 1) * ldA] op expr';
           }
        where expr' has its own __SLICE2D__(...) expanded in the same way.

      Pseudo-slice shorthand (no macros):

        A(i_start, i_end) = zero;
        -> 1D zero-fill loop

        A(i_start, i_end, j_start, j_end) = zero;
        -> 2D zero-fill loop with ldA inferred as "ld" + array_name.lower().
    """
    import re

    def find_matching_paren(s, start):
        """Find the index of the closing paren matching the one at 'start'."""
        depth = 0
        for i in range(start, len(s)):
            if s[i] == "(":
                depth += 1
            elif s[i] == ")":
                depth -= 1
                if depth == 0:
                    return i
        return -1

    def split_args(s):
        """Split a comma-separated argument string at top level."""
        parts = []
        current = []
        depth = 0
        for ch in s:
            if ch == "(":
                depth += 1
                current.append(ch)
            elif ch == ")":
                depth -= 1
                current.append(ch)
            elif ch == "," and depth == 0:
                part = "".join(current).strip()
                if part:
                    parts.append(part)
                current = []
            else:
                current.append(ch)
        if current:
            part = "".join(current).strip()
            if part:
                parts.append(part)
        return parts

    # Helper to build "(i - 1) + (col - 1) * ldname"
    def make_index_expr(col_expr, ldname, loop_var="i"):
        """Build index expression: (loop_var - 1) + (col_expr - 1) * ldname.

        If col_expr has top-level arithmetic operators (+, -, *, /),
        wrap it in an extra pair of parentheses before subtracting 1.

        Examples:
          col_expr = "j"
            -> (j - 1)
          col_expr = "n + 1"
            -> ((n + 1) - 1)
          col_expr = "icolz + nsl - 2"
            -> ((icolz + nsl - 2) - 1)
        """
        col_expr = col_expr.strip()
        if col_expr == "1":
            # Special case: column 1 -> (1 - 1)
            col_term = "(1 - 1)"
        else:
            core = col_expr
            # If the column expression has top-level +, -, * or /,
            # protect it with an extra level of parentheses so that
            # we get ((expr) - 1) instead of expr - 1 gluing together.
            if _has_top_level_arith_op(core):
                core = f"({core})"
            col_term = f"({core} - 1)"
        return f"({loop_var} - 1) + {col_term} * {ldname}"

    def replace_slice2d_with_index(expr, loop_var="i"):
        """Replace A[__SLICE2D__(start, end, col, ld)] with A[(loop_var - 1) + (col - 1) * ld]."""
        pattern = re.compile(
            r"([A-Za-z_][A-Za-z0-9_]*)"
            r"\s*\[\s*__SLICE2D__\s*\("
        )
        result = []
        pos = 0

        while True:
            m = pattern.search(expr, pos)
            if not m:
                result.append(expr[pos:])
                break

            result.append(expr[pos:m.start()])
            array_name = m.group(1)

            paren_start = m.end() - 1
            paren_end = find_matching_paren(expr, paren_start)
            if paren_end < 0:
                # Give up if parentheses are unbalanced
                result.append(expr[m.start():])
                break

            args_str = expr[paren_start + 1:paren_end]
            args = split_args(args_str)
            if len(args) != 4:
                # Not the expected form: leave as-is
                result.append(expr[m.start():paren_end + 1])
                pos = paren_end + 1
                continue

            # args: start, end, col, ldname
            col_expr = args[2]
            ldname = args[3]
            index_expr = make_index_expr(col_expr, ldname, loop_var=loop_var)

            # Skip past closing ")]"
            rest_pos = paren_end + 1
            while rest_pos < len(expr) and expr[rest_pos] in " \t":
                rest_pos += 1
            if rest_pos < len(expr) and expr[rest_pos] == "]":
                rest_pos += 1

            result.append(f"{array_name}[{index_expr}]")
            pos = rest_pos

        return "".join(result)

    # Pattern for 1D slice:   A[__SLICE__(...)]
    pattern_1d = re.compile(
        r"(\s*)"
        r"([A-Za-z_][A-Za-z0-9_]*)"
        r"\s*\[\s*__SLICE__\s*\("
    )

    # Pattern for 2D slice:   A[__SLICE2D__(...)]
    pattern_2d = re.compile(
        r"(\s*)"
        r"([A-Za-z_][A-Za-z0-9_]*)"
        r"\s*\[\s*__SLICE2D__\s*\("
    )

    def rewrite_line(line: str) -> str:
        # --------------------------------------------------------------
        # 0) Pseudo-slice shorthand:  A(i1, i2) = zero;  or  A(i1,i2,j1,j2) = zero;
        # --------------------------------------------------------------
        m_zero = re.match(
            r"(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)\s*=\s*zero\s*;",
            line,
        )
        if m_zero:
            leading_ws = m_zero.group(1)
            array_name = m_zero.group(2)
            args_str = m_zero.group(3)
            args = split_args(args_str)

            # 1D pseudo slice: A(i_start, i_end) = zero;
            if len(args) == 2:
                row_start, row_end = [a.strip() for a in args]
                for_loop = (
                    f"{leading_ws}for (INTEGER i_ = {row_start}; "
                    f"i_ <= {row_end}; i_++) {{\n"
                    f"{leading_ws}    {array_name}[(i_ - 1)] = zero;\n"
                    f"{leading_ws}}}"
                )
                return for_loop

            # 2D pseudo slice: A(i_start, i_end, j_start, j_end) = zero;
            if len(args) == 4:
                row_start, row_end, col_start, col_end = [
                    a.strip() for a in args]
                ldname = "ld" + array_name.lower()
                for_loop = (
                    f"{leading_ws}for (INTEGER l_ = {row_start}; "
                    f"l_ <= {row_end}; l_++) {{\n"
                    f"{leading_ws}    for (INTEGER m_ = {col_start}; "
                    f"m_ <= {col_end}; m_++) {{\n"
                    f"{leading_ws}        {array_name}[(l_ - 1) + "
                    f"(m_ - 1) * {ldname}] = zero;\n"
                    f"{leading_ws}    }}\n"
                    f"{leading_ws}}}"
                )
                return for_loop

        # --------------------------------------------------------------
        # 1) 2D slice with __SLICE2D__(...)
        # --------------------------------------------------------------
        m = pattern_2d.search(line)
        if m:
            leading_ws = m.group(1)
            array_name = m.group(2)

            slice_paren_start = m.end() - 1
            slice_paren_end = find_matching_paren(line, slice_paren_start)
            if slice_paren_end < 0:
                return line

            args_str = line[slice_paren_start + 1:slice_paren_end]
            args = split_args(args_str)
            if len(args) != 4:
                return line

            start_expr, end_expr, col_expr, ldname = args

            # After __SLICE2D__(...) we expect: "] op= value ;"
            rest_start = slice_paren_end + 1
            rest = line[rest_start:]

            assign_match = re.match(r"\s*\]\s*([+\-*/]?=)\s*", rest)
            if not assign_match:
                return line

            assign_op = assign_match.group(1)
            value_start = rest_start + assign_match.end()

            # Find semicolon at top level to terminate RHS.
            depth_paren = 0
            depth_bracket = 0
            semicolon_pos = -1
            for i in range(value_start, len(line)):
                ch = line[i]
                if ch == "(":
                    depth_paren += 1
                elif ch == ")":
                    depth_paren -= 1
                elif ch == "[":
                    depth_bracket += 1
                elif ch == "]":
                    depth_bracket -= 1
                elif ch == ";" and depth_paren == 0 and depth_bracket == 0:
                    semicolon_pos = i
                    break

            if semicolon_pos < 0:
                return line

            value_expr = line[value_start:semicolon_pos].strip()
            value_expr = replace_slice2d_with_index(value_expr, loop_var="i_")

            index_expr = make_index_expr(col_expr, ldname, loop_var="i_")

            for_loop = (
                f"{leading_ws}for (INTEGER i_ = {start_expr}; "
                f"i_ <= {end_expr}; i_++) {{ "
                f"{array_name}[{index_expr}] {assign_op} {value_expr}; }}"
            )

            return line[:m.start()] + for_loop + line[semicolon_pos + 1:]

        # --------------------------------------------------------------
        # 2) 1D slice with __SLICE__(...)
        # --------------------------------------------------------------
        m = pattern_1d.search(line)
        if m:
            leading_ws = m.group(1)
            array_name = m.group(2)

            slice_paren_start = m.end() - 1
            slice_paren_end = find_matching_paren(line, slice_paren_start)
            if slice_paren_end < 0:
                return line

            args_str = line[slice_paren_start + 1:slice_paren_end]
            args = split_args(args_str)
            if len(args) not in (2, 3):
                return line

            start_expr = args[0].strip()
            end_expr = args[1].strip()
            step_expr = "1"
            if len(args) == 3:
                step_expr = args[2].strip()

            rest_start = slice_paren_end + 1
            rest = line[rest_start:]

            assign_match = re.match(r"\s*\]\s*([+\-*/]?=)\s*", rest)
            if not assign_match:
                return line

            assign_op = assign_match.group(1)
            value_start = rest_start + assign_match.end()

            depth_paren = 0
            depth_bracket = 0
            semicolon_pos = -1
            for i in range(value_start, len(line)):
                ch = line[i]
                if ch == "(":
                    depth_paren += 1
                elif ch == ")":
                    depth_paren -= 1
                elif ch == "[":
                    depth_bracket += 1
                elif ch == "]":
                    depth_bracket -= 1
                elif ch == ";" and depth_paren == 0 and depth_bracket == 0:
                    semicolon_pos = i
                    break

            if semicolon_pos < 0:
                return line

            value_expr = line[value_start:semicolon_pos].strip()

            loop_header = None
            if len(args) == 2:
                loop_header = f"for (INTEGER i_ = {start_expr}; i_ <= {end_expr}; i_++)"
            else:
                # Strided slice: a(i:j:k) assignment. Handle constant and unknown step signs.
                m_int = re.fullmatch(r"[+-]?\d+", step_expr)
                if m_int:
                    step_val = int(step_expr)
                    if step_val == 0:
                        return line
                    if step_val > 0:
                        loop_header = f"for (INTEGER i_ = {start_expr}; i_ <= {end_expr}; i_ += {step_expr})"
                    else:
                        loop_header = f"for (INTEGER i_ = {start_expr}; i_ >= {end_expr}; i_ += {step_expr})"
                else:
                    loop_header = (
                        f"for (INTEGER i_ = {start_expr}; "
                        f"(({step_expr}) > 0 ? (i_ <= {end_expr}) : (i_ >= {end_expr})); "
                        f"i_ += {step_expr})"
                    )

            for_loop = (
                f"{leading_ws}{loop_header} {{ "
                f"{array_name}[i_ - 1] {assign_op} {value_expr}; }}"
            )
            return line[:m.start()] + for_loop + line[semicolon_pos + 1:]

        # No slice pattern matched: return line unchanged.
        return line

    if isinstance(lines, list):
        return [rewrite_line(line) for line in lines]
    elif isinstance(lines, str):
        return "\n".join(rewrite_line(l) for l in lines.split("\n"))
    return lines


def _postprocess_complex_zero_initializers(lines):
    """Rewrite COMPLEX zero initializers using COMPLEX(0.0, 0.0).

    Examples:
      COMPLEX ztemp = (0.0, 0.0);       -> COMPLEX ztemp = COMPLEX(0.0, 0.0);
      const COMPLEX zero = (0.0, 0.0);  -> const COMPLEX zero = COMPLEX(0.0, 0.0);
    """
    pat_var = re.compile(
        r'^(\s*COMPLEX\s+[A-Za-z_][A-Za-z0-9_]*\s*=\s*)\(\s*0\.0\s*,\s*0\.0\s*\);'
    )
    pat_const = re.compile(
        r'^(\s*const\s+COMPLEX\s+[A-Za-z_][A-Za-z0-9_]*\s*=\s*)\(\s*0\.0\s*,\s*0\.0\s*\);'
    )
    out = []
    for line in lines:
        line = pat_var.sub(r'\1COMPLEX(0.0, 0.0);', line)
        line = pat_const.sub(r'\1COMPLEX(0.0, 0.0);', line)
        out.append(line)
    return out

def _postprocess_fix_scalar_str_length_from_slices(lines):
    """Fix wrong fem::str<N> length for scalar locals when slices exceed N.

    Why this exists:
      Some Fortran forms like "CHARACTER(LEN=3) PATH" may not surface as
      fdecl.size_tokens depending on the upstream parser. When that happens we may
      accidentally emit fem::str<1> based on later scalar assignments, while the
      code still contains substring operations like PATH(2:3). That yields
      -Wstringop-overflow in optimized builds.

    Strategy:
      1) Collect scalar fem::str<N> declarations (locals).
      2) Collect maximum substring end index used for those identifiers.
      3) If max_end > N, widen the declaration to fem::str<max_end>.

    This is a conservative postprocess: it only widens and never shrinks.
    """

    # Match scalar fem::str<N> declarations (optionally initialized).
    decl_re = re.compile(
        r"\b(?P<type>fem::str<(?P<n>\d+)>)\s+"
        r"(?P<name>[A-Za-z_]\w*)\b(?!\s*\[)"
    )

    # Match substring usage: name(2, 3)
    slice_re = re.compile(
        r"\b(?P<name>[A-Za-z_]\w*)\s*\(\s*(?P<first>\d+)\s*,\s*(?P<last>\d+)\s*\)"
    )

    declared_len = {}
    for line in lines:
        m = decl_re.search(line)
        if not m:
            continue
        name = m.group("name")
        # Skip references and declarations in function parameter lists.
        if "&" in line[: m.start("name")]:
            continue
        declared_len[name] = int(m.group("n"))

    if not declared_len:
        return lines

    need_len = dict(declared_len)
    for line in lines:
        for m in slice_re.finditer(line):
            name = m.group("name")
            if name not in need_len:
                continue
            last = int(m.group("last"))
            if last > need_len[name]:
                need_len[name] = last

    out = []
    for line in lines:
        m = decl_re.search(line)
        if not m:
            out.append(line)
            continue
        name = m.group("name")
        old_n = int(m.group("n"))
        new_n = need_len.get(name, old_n)
        if new_n <= old_n:
            out.append(line)
            continue
        out.append(line.replace(f"fem::str<{old_n}>", f"fem::str<{new_n}>", 1))
    return out

def _fix_fortran_externals(src):
    """Downgrade simple Fortran 90 EXTERNAL declarations to F77 style.

    Examples:
      EXTERNAL :: XERBLA, DHGEQZ
        -> EXTERNAL XERBLA, DHGEQZ
      DOUBLE PRECISION, EXTERNAL :: DLAMCH, DLANHS
        -> DOUBLE PRECISION DLAMCH, DLANHS
    """
    import re
    lines = src.splitlines(True)
    out = []
    for line in lines:
        # 1) EXTERNAL :: ...  ->  EXTERNAL ...
        line2 = re.sub(r'(?i)\bEXTERNAL\s*::', 'EXTERNAL', line)

        # 2) <type>, EXTERNAL ...  ->  <type> ...
        #    Handles typical LAPACK patterns:
        #      DOUBLE PRECISION, EXTERNAL :: DLAMCH
        #      LOGICAL, EXTERNAL :: LSAME
        #      INTEGER, EXTERNAL :: ILAENV
        line2 = re.sub(
            r'(?i)(\bDOUBLE\s+PRECISION|\bINTEGER|\bLOGICAL|\bREAL|\bCOMPLEX|\bCHARACTER)\s*,\s*EXTERNAL\b',
            r'\1',
            line2,
        )
        out.append(line2)
    return ''.join(out)


def _fix_fortran_end_statements(src: str) -> str:
    """Downgrade F90-style typed END statements to a bare END.

    Some fixed-form parsers accept only 'END' to terminate a program unit.
    Example:
        END SUBROUTINE DSYSWAPR  ->  END
    """
    import re
    lines = src.splitlines(True)
    out = []
    end_pat = re.compile(
        r'^(?P<indent>\s*)end\s+'
        r'(?P<kind>subroutine|function|program|module|block\s+data)\b.*$',
        flags=re.IGNORECASE,
    )
    for line in lines:
        if not line:
            out.append(line)
            continue

        # Skip fixed-form comment lines (column-1 'c'/'C'/'*')
        if line[0] in ("c", "C", "*"):
            out.append(line)
            continue

        # Preserve inline comments started by '!' (common in LAPACK)
        bang = line.find("!")
        if bang >= 0:
            code, comment = line[:bang], line[bang:]
        else:
            code, comment = line, ""

        m = end_pat.match(code.rstrip("\r\n"))
        if not m:
            out.append(line)
            continue

        indent = m.group("indent")
        if comment:
            out.append(f"{indent}END{comment}")
        else:
            eol = "\r\n" if line.endswith("\r\n") else (
                "\n" if line.endswith("\n") else "")
            out.append(f"{indent}END{eol}")
    return "".join(out)


def _drop_fortran_intrinsic_statements(src: str, *, source_form: str = "fixed") -> str:
    """Drop Fortran INTRINSIC statements (and their continuation lines).

    Example that breaks some parsers:
        INTRINSIC :: abs, sign, sqrt

    For C++ translation we do not need INTRINSIC markers, so we remove them.
    """
    intrinsic_re = re.compile(r"^\s*intrinsic\b", flags=re.IGNORECASE)

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    def is_full_line_comment_or_blank(raw: str) -> bool:
        if raw.strip() == "":
            return True
        s = raw.lstrip()
        if s.startswith("!"):
            return True
        # Fixed-form comment in column 1
        if raw and raw[0] in ("c", "C", "*", "!"):
            return True
        return False

    def split_code_comment(raw: str):
        idx = raw.find("!")
        if idx >= 0:
            return raw[:idx], raw[idx:]
        return raw, ""

    def is_fixed_form_continuation(raw: str) -> bool:
        # Fixed-form continuation marker is any non-blank in column 6.
        if len(raw) < 6:
            return False
        if raw and raw[0] in ("c", "C", "*", "!"):
            return False
        # Guard: treat as fixed-form only if columns 1-5 are a label field
        # (spaces/digits/tabs). This prevents false positives on free-form lines like
        # "   if(...)" where raw[5] may be '(' and would otherwise be misread as a
        # fixed-form continuation.
        prefix = raw[:5]
        for ch in prefix:
            if ch not in " \t0123456789":
                return False
        return raw[5] not in (" ", "\t")

    lines = src.splitlines(True)
    out = []
    i = 0
    skipping = False
    prev_trailing_amp = False

    while i < len(lines):
        raw, eol = split_eol(lines[i])

        if skipping:
            # Keep unrelated comment/blank lines as-is.
            if is_full_line_comment_or_blank(raw):
                out.append(raw + eol)
                i += 1
                continue

            lstr = raw.lstrip()
            is_free_cont = lstr.startswith("&")
            sf = (source_form or "fixed").lower()
            use_fixed = (sf != "free")
            is_cont = prev_trailing_amp or is_free_cont or (
                use_fixed and is_fixed_form_continuation(raw))

            if is_cont:
                code, _comment = split_code_comment(raw)
                prev_trailing_amp = code.rstrip().endswith("&")
                i += 1
                continue

            # First non-continuation line: stop skipping and re-process it normally.
            skipping = False
            prev_trailing_amp = False
            continue

        # Normal mode
        if is_full_line_comment_or_blank(raw):
            out.append(raw + eol)
            i += 1
            continue

        code, _comment = split_code_comment(raw)
        if intrinsic_re.match(code):
            prev_trailing_amp = code.rstrip().endswith("&")
            skipping = True
            i += 1
            continue

        out.append(raw + eol)
        i += 1

    return "".join(out)


def _strip_openmp_blocks(src: str) -> str:
    """Strip OpenMP-only code guarded by preprocessor directives.

    LAPACK test sources sometimes contain FPP directives like:

        #if defined(_OPENMP)
          use omp_lib
          ...
        #endif

    The FABLE reader is not a full Fortran preprocessor, and fixed-form parsing
    cannot handle free-form "USE" statements. For MPLAPACK conversion we treat
    OpenMP as disabled and keep the non-OpenMP branch (if any).

    This is a lightweight OpenMP-focused filter; it is not a general C preprocessor.
    """

    pp_re = re.compile(
        r"^\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b(.*)$", re.IGNORECASE)
    use_omp_re = re.compile(
        r"^\s*use\s+(?:,\s*intrinsic\s*::\s*)?omp_lib\b", re.IGNORECASE)
    omp_call_re = re.compile(r"\bomp_(get|set)_[a-z0-9_]+\b", re.IGNORECASE)

    def eval_openmp_cond(kind: str, rest: str) -> bool:
        """Evaluate a preprocessor condition that mentions _OPENMP.

        Policy: treat _OPENMP as NOT defined.
        """
        txt = (rest or "").strip()
        low = txt.lower()

        if kind == "ifdef":
            return False
        if kind == "ifndef":
            return True

        # #if / #elif cases
        if re.search(r"!\s*defined\s*\(\s*_openmp\s*\)", low) or re.search(r"!\s*defined\s+_openmp\b", low):
            return True
        if re.search(r"defined\s*\(\s*_openmp\s*\)", low) or re.search(r"defined\s+_openmp\b", low):
            return False
        if re.search(r"!\s*_openmp\b", low):
            return True

        # Fallback: any remaining expression that mentions _OPENMP is treated as false.
        return False

    lines = src.splitlines(True)
    out = []

    # Stack of OpenMP conditional blocks.
    # Each frame: {"depth": int, "keep": bool, "selected": bool}
    stack = []

    def currently_keeping() -> bool:
        return all(fr["keep"] for fr in stack)

    for line in lines:
        # Always drop USE omp_lib; it breaks fixed-form parsing.
        if use_omp_re.match(line):
            continue

        m = pp_re.match(line)
        if m:
            kw = m.group(1).lower()
            rest = (m.group(2) or "").strip()
            mentions_openmp = ("_openmp" in rest.lower())

            # Enter an OpenMP conditional block.
            if kw in ("if", "ifdef", "ifndef") and mentions_openmp:
                cond = eval_openmp_cond(kw, rest)
                stack.append({"depth": 1, "keep": bool(
                    cond), "selected": bool(cond)})
                continue

            if stack:
                fr = stack[-1]

                # Track nested conditionals inside an OpenMP block.
                if kw in ("if", "ifdef", "ifndef"):
                    fr["depth"] += 1
                    # Keep nested directives only if we keep this branch.
                    if currently_keeping():
                        out.append(line)
                    continue

                if kw == "endif":
                    fr["depth"] -= 1
                    if fr["depth"] == 0:
                        stack.pop()
                    continue

                # Switch branches for the top level of the OpenMP conditional.
                if fr["depth"] == 1 and kw in ("elif", "else"):
                    if kw == "else":
                        if not fr["selected"]:
                            fr["keep"] = True
                            fr["selected"] = True
                        else:
                            fr["keep"] = False
                    else:
                        if not fr["selected"]:
                            # If the elif condition still mentions _OPENMP, evaluate it.
                            # Otherwise assume it is the intended non-OpenMP fallback.
                            take = eval_openmp_cond("elif", rest) if (
                                "_openmp" in rest.lower()) else True
                            fr["keep"] = bool(take)
                            if take:
                                fr["selected"] = True
                        else:
                            fr["keep"] = False
                    continue

                # Other directives inside OpenMP blocks: keep only if active.
                if currently_keeping():
                    out.append(line)
                continue

            # Not in an OpenMP block: keep directives as-is.
            out.append(line)
            continue

        # Regular Fortran line
        if stack and not currently_keeping():
            continue

        # Optionally drop direct omp_* calls even if not guarded.
        if omp_call_re.search(line):
            continue

        out.append(line)

    return "".join(out)


def _detect_fortran_source_form(file_name: str, src: str) -> str:
    """Detect Fortran source form: 'free' (F90+) or 'fixed' (F77-style).

    Priority:
      1) File extension (strong hint)
      2) Lightweight content heuristics (tie-break / unknown extensions)

    This is intentionally conservative: when unsure, prefer 'fixed' because
    the historical FABLE reader behavior is typically fixed-form oriented.
    """
    ext = os.path.splitext(file_name)[1].lower()

    free_exts = {".f90", ".f95", ".f03", ".f08", ".f18"}
    fixed_exts = {".f", ".for", ".ftn", ".f77"}

    form_hint = None
    if ext in free_exts:
        form_hint = "free"
    elif ext in fixed_exts:
        form_hint = "fixed"

    # Scan only the head of the file; we want this to be cheap.
    head_lines = src.splitlines()[:200]

    free_hits = 0
    fixed_hits = 0

    for line in head_lines:
        if not line:
            continue

        # Fixed-form comment in column 1
        if line[0] in ("c", "C", "*"):
            fixed_hits += 2

        # Fixed-form label pattern (1-5 digits, then whitespace)
        # Example: " 123  CONTINUE" or "10    FORMAT(...)"
        if re.match(r"^\s{0,5}\d{1,5}\s+", line):
            fixed_hits += 1

        s = line.strip().lower()
        if not s:
            continue

        # Strong free-form tokens
        if "::" in s:
            free_hits += 2
        if s.startswith("use ") or s.startswith("use,"):
            free_hits += 2
        if s.startswith("module ") or s.startswith("contains"):
            free_hits += 2
        if s.startswith("select case"):
            free_hits += 2
        if s.startswith("end subroutine") or s.startswith("end function") or s.startswith("end module"):
            free_hits += 2

    # Decision logic:
    # - Extension wins unless the content strongly contradicts it.
    # - For unknown extensions, pick the higher score; default to fixed.
    if form_hint == "free":
        if fixed_hits >= free_hits + 6 and fixed_hits >= 6:
            return "fixed"
        return "free"
    if form_hint == "fixed":
        if free_hits >= fixed_hits + 3 and free_hits >= 3:
            return "free"
        return "fixed"

    if free_hits > fixed_hits:
        return "free"
    return "fixed"


def _fix_fortran_use_la_constants(src: str) -> str:
    """Replace 'use LA_CONSTANTS, only: ...' with local declarations.

    This is a pragmatic shim for FABLE's limited F90 parsing.
    We keep the original USE lines as comments and inject local constants.
    A marker comment '### MUST BE FIXED' is added intentionally.
    """
    import re

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    def code_part(raw: str) -> str:
        return raw.split("!", 1)[0]

    USE_HEAD_RE = re.compile(r"^\s*use\s+la_constants\b", flags=re.IGNORECASE)

    lines = src.splitlines(True)
    out = []
    i = 0
    found = False
    wp_kind = None  # "dp" or "sp"
    imported = set()

    while i < len(lines):
        raw, eol = split_eol(lines[i])
        if not USE_HEAD_RE.match(code_part(raw).strip()):
            out.append(raw + eol)
            i += 1
            continue

        # Collect the whole USE statement (free-form & continuations).
        found = True
        block = []
        while i < len(lines):
            r, e = split_eol(lines[i])
            block.append(r)
            cp = code_part(r).rstrip()
            i += 1
            if cp.endswith("&"):
                continue
            # Next line can be a continuation starting with '&' or 'only:'
            if i < len(lines):
                nxt = code_part(lines[i]).lstrip()
                if nxt.startswith("&") or nxt.lower().startswith("only:"):
                    continue
            break

        # Parse "only:" list if present.
        flat = " ".join([code_part(b).replace("&", " ") for b in block])
        m_only = re.search(r"\bonly\s*:\s*(.*)$", flat, flags=re.IGNORECASE)
        if m_only:
            only_list = m_only.group(1)
            # Split by commas (good enough here).
            items = [x.strip() for x in only_list.split(",") if x.strip()]
            for it in items:
                # Support renaming: lhs=>rhs
                if "=>" in it:
                    lhs, rhs = [x.strip().lower() for x in it.split("=>", 1)]
                    imported.add(lhs)
                    if lhs == "wp":
                        if rhs in ("dp", "sp"):
                            wp_kind = rhs
                else:
                    imported.add(it.strip().lower())

        # Emit replacement block.
        # Use free-form comments ('!'); if later lowered to fixed-form,
        # your layout normalizer can convert them to 'C' comments.
        out.append(
            "! ### MUST BE FIXED: replaced 'use LA_CONSTANTS, only: ...' with local stubs\n")
        for b in block:
            out.append("! original: " + b.strip() + "\n")

        # Decide precision if not found explicitly.
        if wp_kind is None:
            # Fallback heuristic: if we imported dp-like names, assume dp.
            wp_kind = "dp" if ("dp" in flat.lower(
            ) or "dzero" in flat.lower() or "done" in flat.lower()) else "sp"

        # IMPORTANT:
        #   Do NOT inject Rlamch() here.
        #   This preprocessing runs before fable.read builds declaration tables,
        #   and an undeclared identifier inside PARAMETER assignments can crash
        #   get_fdecl() with KeyError (e.g. 'rlamch').
        #
        # Use standard LAPACK typed functions DLAMCH / SLAMCH with explicit
        # declarations, then let the C++ name-map rename them (e.g. dlamch -> Rlamch).
        if wp_kind == "dp":
            out.append("DOUBLE PRECISION zero, half, one, safmin, safmax\n")
            out.append("PARAMETER (zero = 0.0D0)\n")
            out.append("PARAMETER (half = 0.5D0)\n")
            out.append("PARAMETER (one  = 1.0D0)\n")
            # Use intrinsics to avoid external/specification functions in PARAMETER.
            out.append("PARAMETER (safmin = one //UNHANDLED )\n")
            out.append("PARAMETER (safmax = one //UNHANDLED)\n")
        else:
            out.append("REAL zero, half, one, safmin, safmax\n")
            out.append("PARAMETER (zero = 0.0E0)\n")
            out.append("PARAMETER (half = 0.5E0)\n")
            out.append("PARAMETER (one  = 1.0E0)\n")
            out.append("PARAMETER (safmin = one //UNHANDLED)\n")
            out.append("PARAMETER (safmax = one //UNHANDLED)\n")

    src2 = "".join(out)
    if not found:
        return src2

    # Rewrite REAL(KIND=wp) / REAL(wp) and numeric literals with _wp suffix.
    if wp_kind == "dp":
        src2 = re.sub(r"\breal\s*\(\s*kind\s*=\s*wp\s*\)",
                      "DOUBLE PRECISION", src2, flags=re.IGNORECASE)
        src2 = re.sub(r"\breal\s*\(\s*wp\s*\)",
                      "DOUBLE PRECISION", src2, flags=re.IGNORECASE)
        # 1.0_wp -> 1.0D0, 1e-3_wp -> 1D-3

        def _wp_lit(m):
            num = m.group("num")
            num = re.sub(r"[eE]([+\-]?\d+)", r"D\1", num)
            return num if ("D" in num or "d" in num) else (num + "D0")
        src2 = re.sub(r"(?P<num>(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+\-]?\d+)?)\s*_wp\b",
                      _wp_lit, src2, flags=re.IGNORECASE)
    else:
        src2 = re.sub(r"\breal\s*\(\s*kind\s*=\s*wp\s*\)",
                      "REAL", src2, flags=re.IGNORECASE)
        src2 = re.sub(r"\breal\s*\(\s*wp\s*\)", "REAL",
                      src2, flags=re.IGNORECASE)

        def _wp_lit_sp(m):
            num = m.group("num")
            num = re.sub(r"[dD]([+\-]?\d+)", r"E\1", num)
            return num if ("E" in num or "e" in num) else (num + "E0")
        src2 = re.sub(r"(?P<num>(?:\d+\.\d*|\d*\.\d+|\d+)(?:[dD][+\-]?\d+)?)\s*_wp\b",
                      _wp_lit_sp, src2, flags=re.IGNORECASE)

    return src2


def _fix_fortran_use_la_xisnan(src: str) -> str:
    """Drop/Comment-out 'USE LA_XISNAN' statements.

    Some LAPACK/BLAS F90 sources contain lines like:
        USE LA_XISNAN, ONLY: LA_ISNAN

    FABLE's fixed-form reader does not understand 'USE' statements, so we
    comment them out during preprocessing. The actual NaN check is handled
    later in C++ by mapping LA_ISNAN(...) into Mla_isnan(...).
    """
    import re

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    def code_part(raw: str) -> str:
        return raw.split("!", 1)[0]

    USE_HEAD_RE = re.compile(
        r"^\s*use\s*(?:,\s*intrinsic\s*::\s*)?la_xisnan\b",
        flags=re.IGNORECASE,
    )

    lines = src.splitlines(True)
    out = []
    i = 0

    while i < len(lines):
        raw, eol = split_eol(lines[i])

        # Keep fixed-form whole-line comments intact.
        if raw and raw[0] in ("c", "C", "*", "!"):
            out.append(raw + eol)
            i += 1
            continue

        if not USE_HEAD_RE.match(code_part(raw).strip()):
            out.append(raw + eol)
            i += 1
            continue

        # Collect the whole USE statement (free-form '&' continuations).
        block = []
        while i < len(lines):
            r, e = split_eol(lines[i])
            block.append(r)
            cp = code_part(r).rstrip()
            i += 1
            if cp.endswith("&"):
                continue
            if i < len(lines):
                nxt = code_part(lines[i]).lstrip()
                if nxt.startswith("&") or nxt.lower().startswith("only:"):
                    continue
            break

        out.append(
            "! ### MUST BE FIXED: removed USE LA_XISNAN (NaN check is handled in C++)\n")
        for b in block:
            out.append("! original: " + b.strip() + "\n")

    return "".join(out)


def _should_lower_free_to_fixed(src: str) -> bool:
    """Heuristic: decide whether to lower free-form source into fixed-form."""
    low = src.lower()

    # Strong F90 markers
    if "::" in low:
        return True
    if re.search(r"(?m)^\s*use\b", low):
        return True
    if re.search(r"(?m)^\s*select\s+case\b", low):
        return True

    # Free-form continuation: trailing '&' before any '!' comment
    for line in src.splitlines():
        s = line.lstrip()
        if not s or s.startswith("!"):
            continue
        code = line.split("!", 1)[0]
        if code.rstrip().endswith("&"):
            return True

    return False


def _fix_fortran_intrinsic_kind_keyword_arguments(src: str) -> str:
    """Drop F90 keyword arguments like KIND=... in intrinsic calls.

    FABLE's legacy expression tokenizer does not accept keyword arguments
    inside calls, e.g.:
        CMPLX(x, y, KIND=WP)
        REAL(z, KIND=WP)
        INT(i, KIND=KIND(i))

    For the C++ translation we do not need the KIND selector at the call site,
    because target types are fixed by the translated declarations.

    This pass removes top-level arguments matching /^KIND\s*=/ (case-insensitive)
    from a conservative set of intrinsics that commonly appear in LAPACK F90.
    Nested calls are handled recursively within a single physical line.
    """
    import re

    intrinsic_with_kind = {
        "cmplx",
        "real",
        "int",
        "nint",
        "aint",
        "anint",
        "floor",
        "ceiling",
    }

    kind_arg_re = re.compile(r"^kind\s*=", flags=re.IGNORECASE)

    def find_matching_paren(s: str, start: int) -> int:
        depth = 0
        for i in range(start, len(s)):
            ch = s[i]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    return i
        return -1

    call_re = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\(")

    def rewrite_segment(code: str) -> str:
        out = []
        i = 0
        while i < len(code):
            m = call_re.search(code, i)
            if not m:
                out.append(code[i:])
                break

            name = m.group(1)
            paren_start = m.end() - 1
            paren_end = find_matching_paren(code, paren_start)
            if paren_end < 0:
                out.append(code[i:])
                break

            out.append(code[i:m.start(1)])
            head = code[m.start(1):paren_start]

            inner = code[paren_start + 1:paren_end]
            inner = rewrite_segment(inner)

            if name.lower() in intrinsic_with_kind and "kind" in inner.lower():
                args = _split_top_level_commas(inner)
                args = [a for a in args if not kind_arg_re.match(a.strip())]
                inner = ", ".join(args)

            out.append(f"{head}({inner})")
            i = paren_end + 1

        return "".join(out)

    def is_full_line_comment(raw: str) -> bool:
        if not raw:
            return False
        # Fixed-form whole-line comment (column 1)
        if raw[0] in ("c", "C", "*", "!"):
            return True
        # Free-form whole-line comment
        if raw.lstrip().startswith("!"):
            return True
        return False

    out_lines = []
    for line in src.splitlines(True):
        if is_full_line_comment(line):
            out_lines.append(line)
            continue

        # Preserve inline comments started by '!'
        bang = line.find("!")
        if bang >= 0:
            code_part = line[:bang]
            comment_part = line[bang:]
        else:
            code_part = line
            comment_part = ""

        out_lines.append(rewrite_segment(code_part) + comment_part)

    return "".join(out_lines)


def _preprocess_fortran_fixed_form(src: str) -> typing.Tuple[str, str]:
    """Preprocess a fixed-form Fortran source.

    Returns: (new_src, emit_form)
      emit_form is 'fixed' or 'free' and controls the temp file suffix.
    """
    new_src = src
    new_src = _fix_fortran_externals(new_src)
    new_src = _fix_fortran_use_la_constants(new_src)
    new_src = _fix_fortran_use_la_xisnan(new_src)
    new_src = _drop_fortran_intrinsic_statements(
        new_src, source_form="fixed")  # drop unconditionally
    new_src = _fix_fortran_iso_fortran_env_real64(
        new_src)  # harmless even if not present
    new_src = _fix_fortran_intrinsic_kind_keyword_arguments(new_src)
    new_src = _fix_fortran_f90_decl_syntax(new_src)
    new_src = _fix_fortran_select_case_to_if(new_src)
    new_src = _fix_fortran_end_statements(new_src)
    return new_src, "fixed"


def _preprocess_fortran_free_form(src: str) -> typing.Tuple[str, str]:
    """Preprocess a free-form Fortran source.

    Returns: (new_src, emit_form)
      emit_form is 'fixed' or 'free' and controls the temp file suffix.
    """
    new_src = src
    new_src = _fix_fortran_externals(new_src)
    new_src = _fix_fortran_use_la_constants(new_src)
    new_src = _fix_fortran_use_la_xisnan(new_src)
    new_src = _drop_fortran_intrinsic_statements(
        new_src, source_form="free")  # drop unconditionally
    new_src = _fix_fortran_intrinsic_kind_keyword_arguments(new_src)
    lower_to_fixed = _should_lower_free_to_fixed(new_src)

    if lower_to_fixed:
        # Lowering pipeline: free-form -> fixed-form
        new_src = _fix_fortran_iso_fortran_env_real64(new_src)
        new_src = _fix_fortran_f90_decl_syntax(new_src)
        new_src = _fix_fortran_select_case_to_if(new_src)
        new_src = _fix_fortran_free_form_ampersand_continuations(new_src)
        new_src = _normalize_free_form_to_fixed_form_layout(new_src)
        new_src = _fix_fortran_end_statements(new_src)
        return new_src, "fixed"

    # If we keep it free-form, do only safe transforms and emit ".f90"
    new_src = _fix_fortran_end_statements(new_src)
    return new_src, "free"


def _split_top_level_commas(s: str):
    """Split by commas, ignoring commas inside parentheses."""
    items = []
    buf = []
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            items.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    tail = "".join(buf).strip()
    if tail:
        items.append(tail)
    return items


def _fix_fortran_f90_decl_syntax(src: str) -> str:
    """
    Rewrite a subset of Fortran 90 declaration syntax into a form that fable's
    fixed-form reader accepts.

    Constraints / policy:
      - Do NOT emit standalone DIMENSION statements.
        Always embed dimensions in the type declaration: "DOUBLE PRECISION Q(M,M)".
      - Handle fixed-form continuation lines (col 6 non-blank, non-'0') as one statement.
      - For ALLOCATABLE deferred-shape arrays (:) or (:,:), try to substitute the
        explicit shape from a later ALLOCATE(...) statement (RAII-friendly).
      - If an ALLOCATE becomes redundant due to explicit-shape substitution, remove the
        entire continued ALLOCATE statement and set STAT=0.
      - Leave DEALLOCATE untouched (expected to remain UNHANDLED downstream).
    """
    import re

    TYPE_RE = re.compile(
        r'^\s*(DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|REAL|INTEGER|LOGICAL|CHARACTER|COMPLEX)\b',
        flags=re.IGNORECASE
    )

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        return line, ""

    def is_full_line_comment(raw: str) -> bool:
        return bool(raw) and raw[0] in ("c", "C", "*", "!")

    def is_fixed_form_continuation(raw: str) -> bool:
        # Fixed-form continuation if column 6 is non-blank and not '0'
        # (col index 5 in 0-based).
        if is_full_line_comment(raw):
            return False
        if len(raw) < 6:
            return False
        cc = raw[5]
        return (cc != " ") and (cc != "0") and (cc != "\t")

    def split_inline_comment(raw: str):
        """Split raw line into (code, comment) at first '!' outside quotes."""
        in_str = False
        quote = None
        for i, ch in enumerate(raw):
            if ch in ("'", '"'):
                if not in_str:
                    in_str = True
                    quote = ch
                elif quote == ch:
                    in_str = False
                    quote = None
            if (not in_str) and ch == "!":
                return raw[:i], raw[i:]
        return raw, ""

    def code_field_fixed_form(code: str) -> str:
        # Columns 7-72 are statement field in fixed form (0-based index 6:).
        # If the line is shorter, fall back to lstrip to avoid empty strings.
        return code[6:] if len(code) >= 6 else code.lstrip()

    def gather_fixed_statement(lines, i: int):
        """
        Gather a full fixed-form statement starting at line i, consuming continuation
        lines (col 6 continuation). Returns:
          (stmt_text, first_inline_comment, first_eol, j_next)
        where j_next is the first line index after the statement.
        """
        raw0, eol0 = split_eol(lines[i])
        code0, cmt0 = split_inline_comment(raw0)
        stmt = code_field_fixed_form(code0).rstrip()

        j = i + 1
        while j < len(lines):
            rawj, _ = split_eol(lines[j])
            if not is_fixed_form_continuation(rawj):
                break
            codej, _ = split_inline_comment(rawj)
            part = code_field_fixed_form(codej).strip()
            if part:
                stmt += " " + part
            j += 1
        return stmt, cmt0, eol0, j

    def find_top_level_eq(text: str) -> int:
        """Find '=' at nesting level 0 (ignore strings and '=>')."""
        depth = 0
        in_str = False
        quote = None
        i = 0
        while i < len(text):
            ch = text[i]
            if ch in ("'", '"'):
                if not in_str:
                    in_str = True
                    quote = ch
                elif quote == ch:
                    in_str = False
                    quote = None
                i += 1
                continue
            if not in_str:
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth = max(0, depth - 1)
                elif ch == "=" and depth == 0:
                    if i + 1 < len(text) and text[i + 1] == ">":
                        i += 2
                        continue
                    return i
            i += 1
        return -1

    def extract_paren_content(text: str, lpar: int):
        """Given text[lpar]=='(', return (content, rpar_index)."""
        depth = 0
        for i in range(lpar, len(text)):
            if text[i] == "(":
                depth += 1
            elif text[i] == ")":
                depth -= 1
                if depth == 0:
                    return text[lpar + 1:i], i
        return text[lpar + 1:], len(text) - 1

    def normalize_base_type(bt: str) -> str:
        low = " ".join(bt.strip().split()).lower()
        if low == "double precision":
            return "DOUBLE PRECISION"
        if low == "double complex":
            return "DOUBLE COMPLEX"
        return bt.strip().upper()

    def normalize_dim_spec(dim: str):
        """
        Replace deferred-shape ':' with a minimal placeholder '1'
        so it becomes explicit-shape for the parser.
        """
        if dim is None:
            return None
        dim = dim.strip()
        if not dim:
            return None
        parts = _split_top_level_commas(dim)
        if not parts:
            return None
        out_parts = []
        changed = False
        for p in parts:
            ps = p.strip()
            if ps == ":":
                out_parts.append("1")
                changed = True
            else:
                out_parts.append(ps)
        return ", ".join(out_parts) if changed else dim

    def parse_allocate_statement(stmt: str):
        """
        Parse ALLOCATE(...) statement; return ([(name, shape)], stat_var) or None.
        """
        low = stmt.lower().lstrip()
        if not low.startswith("allocate"):
            return None
        lpar = stmt.lower().find("allocate")
        lpar = stmt.find("(", lpar)
        if lpar < 0:
            return None
        inside, _ = extract_paren_content(stmt, lpar)
        items = _split_top_level_commas(inside)
        objs = []
        stat_var = None
        for it in items:
            it = it.strip()
            if not it:
                continue
            eq = find_top_level_eq(it)
            if eq != -1:
                key = it[:eq].strip().lower()
                val = it[eq + 1:].strip()
                if key == "stat":
                    stat_var = val
                continue
            lpar2 = it.find("(")
            if lpar2 != -1:
                shape, _ = extract_paren_content(it, lpar2)
                name = it[:lpar2].strip()
                if name:
                    objs.append((name, shape.strip()))
            else:
                name = it.strip()
                if name:
                    objs.append((name, None))
        return objs, stat_var

    lines = src.splitlines(True)

    # Pass 1: collect shapes from continued ALLOCATE statements.
    allocate_shapes = {}
    i = 0
    while i < len(lines):
        rawi, _ = split_eol(lines[i])
        if is_full_line_comment(rawi):
            i += 1
            continue
        stmt, _, _, j = gather_fixed_statement(lines, i)
        parsed = parse_allocate_statement(stmt)
        if parsed:
            objs, _ = parsed
            for name, shape in objs:
                if shape:
                    allocate_shapes.setdefault(name.lower(), shape)
        i = j

    # Pass 2: rewrite continued declarations and optionally remove continued ALLOCATE.
    out_lines = []
    raii_alloc_names = set()  # all ALLOCATABLE vars we rewrite into explicit-shape decls
    indent = "      "

    i = 0
    while i < len(lines):
        rawi, eoli = split_eol(lines[i])

        if is_full_line_comment(rawi):
            out_lines.append(rawi + eoli)
            i += 1
            continue

        stmt, first_cmt, eol0, j = gather_fixed_statement(lines, i)

        # Remove ALLOCATE(...) if it only targets RAII variables.
        parsed_alloc = parse_allocate_statement(stmt)
        if parsed_alloc:
            objs, stat_var = parsed_alloc
            if objs and all((name.lower() in raii_alloc_names) for name, _ in objs):
                out_lines.append(
                    f"!FABLE: ALLOCATE removed (RAII in C++){eol0}")
                if stat_var:
                    out_lines.append(f"{indent}{stat_var} = 0{eol0}")
                i = j
                continue
            # Keep original lines unchanged.
            out_lines.extend(lines[i:j])
            i = j
            continue

        # Rewrite only typed '::' declarations.
        if ("::" not in stmt) or (not TYPE_RE.match(stmt)):
            out_lines.extend(lines[i:j])
            i = j
            continue

        left, right = stmt.split("::", 1)
        left_parts = _split_top_level_commas(left.strip())
        if not left_parts:
            out_lines.extend(lines[i:j])
            i = j
            continue

        base_type = normalize_base_type(left_parts[0])
        attrs_raw = [p.strip() for p in left_parts[1:] if p.strip()]

        dim_attr = None
        is_parameter = False
        is_allocatable = False

        for attr in attrs_raw:
            low = attr.lower().strip()
            if low.startswith("dimension"):
                lpar = attr.find("(")
                if lpar != -1:
                    inside, _ = extract_paren_content(attr, lpar)
                    dim_attr = inside.strip()
            elif low.startswith("parameter"):
                is_parameter = True
            elif low.startswith("allocatable"):
                is_allocatable = True
            else:
                # Drop INTENT/OPTIONAL/VALUE/TARGET/etc.
                pass

        dim_attr_norm = normalize_dim_spec(dim_attr) if dim_attr else None
        decl_items = _split_top_level_commas(right.strip())

        if is_parameter:
            # Emit as:
            #   TYPE name1, name2
            #   PARAMETER (name1=..., name2=...)
            names = []
            for it in decl_items:
                eq = find_top_level_eq(it)
                lhs = it.strip() if eq == -1 else it[:eq].strip()
                if not lhs:
                    continue
                lpar = lhs.find("(")
                if lpar != -1:
                    lhs = lhs[:lpar].strip()
                if lhs:
                    names.append(lhs)

            if names:
                out_lines.append(
                    f"{indent}{base_type} {', '.join(names)}{eol0}")

            # Preserve existing behavior when available.
            if "_contains_machine_const_intrinsics" in globals() and "_FORCE_UNHANDLED_PARAMETER_NAMES" in globals():
                if _contains_machine_const_intrinsics(stmt):
                    for nm in names:
                        _FORCE_UNHANDLED_PARAMETER_NAMES.add(nm.lower())

            # If this PARAMETER statement was continued, we already gathered it into stmt,
            # so we can emit it in one line.
            out_lines.append(
                f"{indent}PARAMETER ({', '.join([it.strip() for it in decl_items if it.strip()])}){eol0}")
            i = j
            continue

        # Non-PARAMETER: embed dims in declarators; emit one decl per variable.
        has_init = any((find_top_level_eq(it.strip()) != -1)
                       for it in decl_items)
        if has_init:
            # Too risky to refactor initialized declarators; keep as-is but drop F90 attrs by removing '::'
            # The simplest safe action is to keep original statement unchanged.
            out_lines.extend(lines[i:j])
            i = j
            continue

        declarators = []
        for it in decl_items:
            it = it.strip()
            if not it:
                continue

            lhs = it
            dim = None

            lpar = lhs.find("(")
            if lpar != -1:
                inside, _ = extract_paren_content(lhs, lpar)
                dim = normalize_dim_spec(inside)
                name = lhs[:lpar].strip()
            else:
                name = lhs.strip()

            if not name:
                continue

            # Apply DIMENSION(...) attribute if present.
            if dim is None and dim_attr_norm is not None:
                dim = dim_attr_norm

            if is_allocatable:
                # Prefer explicit shape from ALLOCATE mapping (RAII-friendly).
                shape = allocate_shapes.get(name.lower())
                if shape:
                    dim = shape.strip()
                # Mark as RAII-converted regardless.
                raii_alloc_names.add(name.lower())
                # If still unknown, force minimal placeholder.
                if dim is None:
                    dim = "1"

            if dim:
                declarators.append(f"{name}({dim})")
            else:
                declarators.append(name)

        for k, decl in enumerate(declarators):
            cmt = first_cmt if k == 0 else ""
            out_lines.append(f"{indent}{base_type} {decl}{cmt}{eol0}")

        i = j

    return "".join(out_lines)


def _fix_fortran_select_case_to_if(src: str) -> str:
    """Rewrite SELECT CASE into IF/ELSE IF/END IF (F77-friendly).

    Supports patterns like:
      SELECT CASE (expr)
        CASE (1)
          ...
        CASE (2,3)
          ...
        CASE (1:5)
          ...
        CASE DEFAULT
          ...
      END SELECT
    Nested SELECT CASE is handled recursively.
    """
    import re

    SELECT_RE = re.compile(
        r'^\s*select\s+case\s*\(\s*(?P<expr>[^)]+?)\s*\)\s*$', flags=re.IGNORECASE)
    CASE_RE = re.compile(
        r'^\s*case\s*(?:(?P<default>default)|\(\s*(?P<sel>.+?)\s*\))\s*$', flags=re.IGNORECASE)
    ENDSEL_RE = re.compile(r'^\s*end\s*select\b.*$', flags=re.IGNORECASE)

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    def split_code_comment(raw: str):
        idx = raw.find("!")
        if idx >= 0:
            return raw[:idx], raw[idx:]
        return raw, ""

    def is_comment_or_blank(raw: str) -> bool:
        if raw.strip() == "":
            return True
        s = raw.lstrip()
        if s.startswith("!"):
            return True
        if raw and raw[0] in ("c", "C", "*", "!"):
            return True
        return False

    def split_top_level_commas(s: str):
        items = []
        buf = []
        depth = 0
        for ch in s:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth = max(0, depth - 1)
            if ch == "," and depth == 0:
                items.append("".join(buf).strip())
                buf = []
            else:
                buf.append(ch)
        tail = "".join(buf).strip()
        if tail:
            items.append(tail)
        return items

    def case_selector_to_cond(expr: str, sel: str) -> str:
        sel = sel.strip()
        # sel is already inside (...) by regex
        parts = split_top_level_commas(sel)
        conds = []
        for p in parts:
            p = p.strip()
            if ":" in p:
                lo, hi = p.split(":", 1)
                lo = lo.strip()
                hi = hi.strip()
                if lo == "" and hi != "":
                    conds.append(f"({expr} <= {hi})")
                elif hi == "" and lo != "":
                    conds.append(f"({expr} >= {lo})")
                else:
                    conds.append(f"({expr} >= {lo} .AND. {expr} <= {hi})")
            else:
                conds.append(f"({expr} == {p})")
        if len(conds) == 1:
            return conds[0]
        return "(" + " .OR. ".join(conds) + ")"

    lines = src.splitlines(True)

    def rewrite_from(i: int):
        raw0, eol0 = split_eol(lines[i])
        code0, _c0 = split_code_comment(raw0)
        msel = SELECT_RE.match(code0.strip())
        if not msel:
            return [lines[i]], i + 1

        indent = re.match(r'^(\s*)', raw0).group(1)
        expr = msel.group("expr").strip()
        nl = eol0 if eol0 else "\n"

        i += 1
        cases = []
        current_sel = None  # None means "not started yet"
        current_body = []

        while i < len(lines):
            raw, eol = split_eol(lines[i])
            code, _comment = split_code_comment(raw)
            s = code.strip()

            if is_comment_or_blank(raw):
                # Keep comments/blanks in the current body.
                current_body.append(raw + eol)
                i += 1
                continue

            # Nested SELECT CASE
            if SELECT_RE.match(s):
                nested_out, i = rewrite_from(i)
                current_body.extend(nested_out)
                continue

            mcase = CASE_RE.match(s)
            if mcase:
                # Flush previous case
                if current_sel is not None:
                    cases.append((current_sel, current_body))
                current_body = []
                if mcase.group("default"):
                    current_sel = "DEFAULT"
                else:
                    current_sel = mcase.group("sel").strip()
                i += 1
                continue

            if ENDSEL_RE.match(s):
                if current_sel is not None:
                    cases.append((current_sel, current_body))
                i += 1
                break

            current_body.append(raw + eol)
            i += 1

        # Build IF/ELSE IF chain
        out = []
        first = True
        for sel, body in cases:
            if sel == "DEFAULT":
                cond = None
            else:
                cond = case_selector_to_cond(expr, sel)

            if first:
                if cond is None:
                    out.append(f"{indent}IF (.TRUE.) THEN{nl}")
                else:
                    out.append(f"{indent}IF ({cond}) THEN{nl}")
                first = False
            else:
                if cond is None:
                    out.append(f"{indent}ELSE{nl}")
                else:
                    out.append(f"{indent}ELSE IF ({cond}) THEN{nl}")
            out.extend(body)

        out.append(f"{indent}END IF{nl}")
        return out, i

    out_all = []
    i = 0
    while i < len(lines):
        raw, _eol = split_eol(lines[i])
        code, _comment = split_code_comment(raw)
        if SELECT_RE.match(code.strip()):
            chunk, i = rewrite_from(i)
            out_all.extend(chunk)
        else:
            out_all.append(lines[i])
            i += 1

    return "".join(out_all)


def _normalize_free_form_to_fixed_form_layout(src: str) -> str:
    """Normalize free-form indentation so fixed-form continuation detection won't misfire.

    FABLE currently writes temporary sources with suffix '.f', so the reader
    applies fixed-form rules. For free-form sources (.f90 etc.), a normal
    indentation can accidentally put a non-space character in column 6, which
    fixed-form interprets as a continuation line. That triggers an assertion
    in the continuation combiner.

    We rewrite each non-comment line into a safe fixed-form layout:
      - label (if present) is right-justified in columns 1-5
      - column 6 is blank for normal lines
      - statement starts at column 7
      - a leading '&' (free-form continuation marker) becomes a fixed-form
        continuation line: '     &...'
    """
    import re

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    LABEL_RE = re.compile(r"^\s*(\d{1,5})\s+(.*)$")

    out = []
    for line in src.splitlines(True):
        raw, eol = split_eol(line)
        raw = raw.expandtabs(8)

        if raw.strip() == "":
            out.append(raw + eol)
            continue

        lstr = raw.lstrip()

        # Free-form comment line -> fixed-form comment line
        if lstr.startswith("!"):
            out.append("C " + lstr[1:] + eol)
            continue

        # Already a fixed-form comment line
        if raw and raw[0] in ("c", "C", "*", "!"):
            out.append(raw + eol)
            continue

        m = LABEL_RE.match(raw)
        if m:
            label = m.group(1)
            stmt = m.group(2).lstrip()
            out.append(label.rjust(5) + " " + stmt + eol)
            continue

        # Continuation line in free form often starts with '&'
        if lstr.startswith("&"):
            stmt = lstr[1:].lstrip()
            out.append("     &" + stmt + eol)
        else:
            out.append("      " + lstr + eol)

    return "".join(out)


def _fix_fortran_iso_fortran_env_real64(src: str) -> str:
    """Make some F90 'iso_fortran_env real64' constructs digestible for FABLE.

    FABLE's declaration parser may treat USE statements as if they were type
    declarations and fail with "Syntax error" at lines like:
        USE, INTRINSIC :: iso_fortran_env, only: real64

    In LAPACK drivers (e.g. DGEDMD), this is typically used only to define:
        INTEGER, PARAMETER :: WP = real64
    and then REAL(KIND=WP) ... / 1.0_wp literals.

    Strategy:
      1) Drop (comment out) the USE line (parser compatibility).
      2) If WP=real64 is present, rewrite REAL(KIND=WP) -> DOUBLE PRECISION
         and rewrite *_wp literals -> *d0 so the inferred types stay correct.
    """
    import re

    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    USE_REAL64_RE = re.compile(
        r'^\s*use\s*,\s*intrinsic\s*::\s*iso_fortran_env\s*,\s*only\s*:\s*real64\s*$',
        flags=re.IGNORECASE,
    )
    WP_REAL64_RE = re.compile(
        r'^\s*integer\s*,\s*parameter\s*::\s*wp\s*=\s*real64\s*$',
        flags=re.IGNORECASE,
    )

    lines = src.splitlines(True)
    out = []
    has_wp_real64 = False

    for line in lines:
        raw, eol = split_eol(line)

        if USE_REAL64_RE.match(raw):
            # Fixed-form safe comment (column 1)
            out.append("C " + raw.lstrip() + eol)
            continue

        if WP_REAL64_RE.match(raw):
            has_wp_real64 = True
            # Remove WP definition; we'll normalize types below.
            out.append("C " + raw.lstrip() + eol)
            continue

        out.append(raw + eol)

    src2 = "".join(out)
    if not has_wp_real64:
        return src2

    # Normalize REAL(KIND=WP/wp) -> DOUBLE PRECISION
    REAL_KIND_WP_RE = re.compile(
        r'\breal\s*\(\s*kind\s*=\s*wp\s*\)', flags=re.IGNORECASE)

    def rewrite_code_only(line: str) -> str:
        # Keep inline '!' comments untouched.
        idx = line.find("!")
        if idx >= 0:
            code, comment = line[:idx], line[idx:]
        else:
            code, comment = line, ""

        code = REAL_KIND_WP_RE.sub("DOUBLE PRECISION", code)

        # Rewrite numeric literals with _wp suffix to double literals.
        # Examples: 1.0_wp -> 1.0d0, 0.95_wp -> 0.95d0, 1e-3_wp -> 1d-3
        def _lit(m: re.Match) -> str:
            num = m.group("num")
            # convert E exponent to D exponent for double literals
            num = re.sub(r'[eE]([+\-]?\d+)', r'D\1', num)
            # if no exponent, append d0
            if "D" in num or "d" in num:
                return num
            return num + "d0"

        WP_LIT_RE = re.compile(
            r'(?P<num>(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+\-]?\d+)?)\s*_(?:wp)\b',
            flags=re.IGNORECASE,
        )
        code = WP_LIT_RE.sub(_lit, code)
        return code + comment

    src2_lines = src2.splitlines(True)
    src2_lines = [rewrite_code_only(ln) for ln in src2_lines]
    return "".join(src2_lines)


def _fix_fortran_free_form_ampersand_continuations(src: str) -> str:
    """Convert free-form '&' continuations into fixed-form continuation lines.

    We sometimes write patched temporary sources with suffix '.f' to force
    the fixed-form parser. Free-form continuation (a trailing '&') is not
    recognized in fixed form, so we rewrite it to a column-6 continuation.

    Rule of thumb implemented here:
      - If a non-comment line ends with '&' (before any inline '! comment'),
        drop that '&' and mark the next non-comment/non-blank line as a
        continuation line by putting '&' in column 6: '     &...'.
      - If the continuation line starts with a leading '&' (free-form style),
        remove it.
    This is intentionally conservative and does not attempt full Fortran parsing.
    """
    def split_eol(line: str):
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        if line.endswith("\r"):
            return line[:-1], "\r"
        return line, ""

    def is_comment_or_blank(raw: str) -> bool:
        if raw.strip() == "":
            return True
        s = raw.lstrip()
        # Free-form whole-line comment
        if s.startswith("!"):
            return True
        # Fixed-form whole-line comment in column 1
        if raw and raw[0] in ("c", "C", "*", "!"):
            return True
        return False

    out = []
    pending_cont = False

    for line in src.splitlines(True):
        raw, eol = split_eol(line)

        # If previous line had trailing '&', rewrite this line as a fixed-form continuation.
        if pending_cont:
            if is_comment_or_blank(raw):
                out.append(raw + eol)
                continue

            s = raw.lstrip()
            if s.startswith("&"):
                s = s[1:].lstrip()
            # Column-6 continuation marker
            raw = "     &" + s
            pending_cont = False

        # Do not treat pure comment lines as candidates for trailing '&'
        if is_comment_or_blank(raw):
            out.append(raw + eol)
            continue

        # Strip trailing '&' only in the code part (before inline '!' comment)
        bang = raw.find("!")
        if bang >= 0:
            code = raw[:bang]
            comment = raw[bang:]
        else:
            code = raw
            comment = ""

        code_r = code.rstrip()
        if code_r.endswith("&"):
            code_r = code_r[:-1].rstrip()
            pending_cont = True

        out.append(code_r + comment + eol)

    return "".join(out)


def _preprocess_fortran_files(file_names):
    """Return (patched_file_names, temp_files) for FABLE parsing.

    The input source form (free/fixed) is detected using the file extension
    plus a lightweight content heuristic. The preprocessing pipeline and the
    temporary file suffix are chosen accordingly.

    - fixed-form input  -> emit suffix ".f" (fixed)
    - free-form input   -> emit suffix ".f90" (free) by default

    If preprocessing does not change the source, the original filename is
    returned unless we need to force the suffix to match the detected form.
    """
    patched = []
    temp_files = []
    for fn in (file_names or []):
        try:
            with open(fn, "r") as f:
                src = f.read()
            original_src = src
            src = _strip_openmp_blocks(src)
        except OSError:
            patched.append(fn)
            continue

        src_form = _detect_fortran_source_form(fn, src)
        if src_form == "free":
            new_src, emit_form = _preprocess_fortran_free_form(src)
        else:
            new_src, emit_form = _preprocess_fortran_fixed_form(src)

        ext = os.path.splitext(fn)[1].lower()
        free_exts = {".f90", ".f95", ".f03", ".f08", ".f18"}
        fixed_exts = {".f", ".for", ".ftn", ".f77"}

        # If the content is unchanged but the extension does not match the
        # detected/desired form, still write a temporary copy with a matching
        # suffix so the downstream reader can pick the correct mode.
        need_temp = (new_src != original_src)
        if emit_form == "free" and ext not in free_exts:
            need_temp = True
        if emit_form == "fixed" and ext not in fixed_exts:
            need_temp = True

        if not need_temp:
            patched.append(fn)
            continue

        suffix = ".f90" if emit_form == "free" else ".f"

        fd, tmp = tempfile.mkstemp(prefix="fable_tmp_", suffix=suffix)
        os.close(fd)
        with open(tmp, "w") as g:
            g.write(new_src)
        patched.append(tmp)
        temp_files.append(tmp)

    return patched, temp_files


def process(
        file_names=None,
        all_fprocs=None,
        top_procedures=None,
        namespace="please_specify",
        include_prefix=None,
        include_guard_suffix=None,
        top_cpp_file_name=None,
        dynamic_parameters=None,
        fortran_file_comments=False,
        fem_do_safe=True,
        arr_nd_size_max=default_arr_nd_size_max,
        inline_all=False,
        common_equivalence_simple=set(),
        suppress_program=False,
        suppress_common=False,
        separate_cmn_hpp=False,
        number_of_function_files=None,
        separate_files_main_namespace={},
        write_separate_files_main_namespace="All",
        separate_files_separate_namespace={},
        write_separate_files_separate_namespace="All",
        ignore_common_and_save=set(),
        force_not_implemented=set(),
        ignore_missing=set(),
        suppress_functions=set(),
        suppress_function_definitions=set(),
        common_report_stringio=None,
        data_values_block_size=8,
        data_specializations=True,
        debug=False):
    assert [file_names, all_fprocs].count(None) == 1

    # Environment override: suppress emitting COMMON/SAVE boilerplate structs.
    if FABLE_SUPPRESS_COMMON:
        suppress_common = True

    # Reset per-run marker set (used to force UNHANDLED initializers for
    # machine-constant-style PARAMETER expressions).
    _FORCE_UNHANDLED_PARAMETER_NAMES.clear()

    if (namespace is None or namespace == "please_specify"):
        namespace = ""  # Disabled: was "placeholder_please_replace"

    import fable.read

    # Keep a copy of the original file names for output paths.
    orig_file_names = file_names
    temp_files = []

    if (all_fprocs is None):
        patched_file_names = None
        if file_names is not None:
            patched_file_names, temp_files = _preprocess_fortran_files(
                file_names)
        all_fprocs = fable.read.process(file_names=patched_file_names)

    if top_cpp_file_name is None and orig_file_names is not None and len(orig_file_names) == 1:
        # Derive output file name from the input file stem (basename without extension),
        # not from the first fproc name.
        #
        # Reason:
        #   For PROGRAM units, the internal procedure name is often prefixed with
        #   "program_" (e.g. PROGRAM DCHKAB -> program_dchkab). Using that name would
        #   produce program_dchkab.cpp, which is not desired here.
        src_path = orig_file_names[0]
        base_name = os.path.splitext(os.path.basename(src_path))[0]
        if not base_name:
            main_fproc = all_fprocs.all_in_input_order[0]
            base_name = main_fproc.name.value
        src_path = orig_file_names[0]
        stem = os.path.splitext(os.path.basename(src_path))[0]
        base_name = convert_function_name_to_mplapack(stem) if stem else None
        if not base_name:
            main_fproc = all_fprocs.all_in_input_order[0]
            base_name = convert_function_name_to_mplapack(
                main_fproc.name.value)
        src_dir = os.path.dirname(src_path)
        if src_dir:
            top_cpp_file_name = os.path.join(src_dir, base_name + ".cpp")
        else:
            top_cpp_file_name = base_name + ".cpp"

    result = []

    def callback(line):
        if (len(result) == 0):
            prev_line = None
        else:
            prev_line = result[-1]
        lines = break_lines(cpp_text=[line+"\n"], prev_line=prev_line)
        if (len(lines) != 0):
            if (debug):
                print("\n".join(lines))
            result.extend(lines)
    #
    need_function_hpp = False
    if (len(separate_files_main_namespace) != 0):
        need_function_hpp = True
    if (number_of_function_files is not None):
        assert number_of_function_files > 0
        need_function_hpp = True
    if (need_function_hpp):
        separate_cmn_hpp = True
    #
    if (include_guard_suffix is not None):
        include_guard(
            callback=callback, namespace=namespace, suffix=include_guard_suffix)
    #
    if (separate_cmn_hpp):
        callback("#define FEM_TRANSLATION_UNIT_WITH_MAIN")
        callback("")
    #

    def include_separate(callback):
        if (len(separate_files_separate_namespace) == 0):
            return False
        for name in sorted(separate_files_separate_namespace.keys()):
            callback('#include "%s.hpp"' % name)
        return True
    #

    def include_with_prefix(name):
        if (include_prefix is None):
            return '#include "%s.hpp"' % name
        return '#include <%s/%s.hpp>' % (include_prefix, name)
    #
    need_using_major_types = False
    if (need_function_hpp):
        callback(include_with_prefix("functions"))
    elif (separate_cmn_hpp):
        callback(include_with_prefix("cmn"))
    else:
        callback(include_fem_hpp)
        need_using_major_types = True
    callback("")
    if (not need_function_hpp):
        if (include_separate(callback=callback)):
            callback("")

    open_namespace(
        callback=callback,
        namespace=namespace,
        using_namespace_major_types=need_using_major_types)
    #
    topological_fprocs = all_fprocs.build_bottom_up_fproc_list_following_calls(
        top_procedures=top_procedures)

    # Mark procedures whose COMMON/SAVE handling should be ignored.
    # We will implement these routines manually (no auto-generated
    # *_save structs or cmn_sve members).
    hard_ignore_common = {
        "zlacon",  # add more names here if needed
        "drotm",
        "drotmg",
        "dlacon",
        "dlasy2",
        "dlaln2",
        # "slacon",
    }

    for fproc in topological_fprocs.bottom_up_list:
        if fproc.name is None:
            continue
        if fproc.name.value.lower() in hard_ignore_common:
            ch = getattr(fproc, "conv_hook", None)
            if ch is None:
                ch = conv_hook_info()
                fproc.conv_hook = ch
            ch.ignore_common_and_save = True

    missing = topological_fprocs.missing_external_fdecls_by_identifier

    # Do not emit stub implementations for missing external functions.
    # All such functions (lsame, xerbla, etc.) must be provided elsewhere.
    # If you ever want to see the list, you can add diagnostic prints here,
    # but no code should be generated into the output.

    #
    dep_cycles = topological_fprocs.dependency_cycles
    if (len(dep_cycles) != 0):
        callback("")
        callback("/* Dependency cycles: %d" % len(dep_cycles))
        for cycle in dep_cycles:
            callback("     " + " ".join(cycle))
        callback(" */")
    #
    if (dynamic_parameters is not None):
        assert len(dynamic_parameters) != 0
        for fproc in topological_fprocs.bottom_up_list:
            for dp_props in dynamic_parameters:
                fdecl = fproc.fdecl_by_identifier.get(dp_props.name)
                if (fdecl is not None):
                    fproc.dynamic_parameters.add(dp_props.name)
    #

    emit_cmn_hpp = (separate_cmn_hpp and not suppress_common)
    if (emit_cmn_hpp):
        cmn_buffer = []
        cmn_callback = cmn_buffer.append
        include_guard(
            callback=cmn_callback, namespace=namespace, suffix="_CMN_HPP")
        cmn_callback(include_fem_hpp)
        cmn_callback("")
        open_namespace(callback=cmn_callback, namespace=namespace)
    elif (suppress_common):
        def cmn_callback(line): pass
    else:
        cmn_callback = callback
    try:
        converted_commons_info = convert_commons(
            callback=cmn_callback,
            separate_cmn_hpp=separate_cmn_hpp,
            topological_fprocs=topological_fprocs,
            dynamic_parameters=dynamic_parameters,
            common_equivalence_simple=common_equivalence_simple,
            common_report_stringio=common_report_stringio)
    except Exception:
        if (not debug):
            raise
        show_traceback()
        converted_commons_info = None
    if (emit_cmn_hpp):
        close_namespace(callback=cmn_callback,
                        namespace=namespace, hpp_guard=True)
        with open("cmn.hpp", "w") as f:
            print("\n".join(break_lines(cpp_text=cmn_buffer)), file=f)
    # Propagate OUT/INOUT-ness through calls.
    # If a variable is forwarded into a callee parameter that requires a mutable
    # actual (non-const ref/pointer), treat it as modified in the caller too.
    _propagate_out_inout_through_calls(topological_fprocs)

    # Defensive: needs_cmn may depend on call graph properties. Recompute after
    # we changed is_modified flags.
    topological_fprocs.each_fproc_update_needs_cmn()
    separate_function_buffers = []
    separate_function_buffer_by_function_name = {}
    for name, identifiers in separate_files_main_namespace.items():
        if (len(identifiers) == 0):
            raise RuntimeError(
                "separate_files_main_namespace: empty list: %s" % name)
        buffer = []
        buffer.append(include_with_prefix("functions"))
        buffer.append("")
        separate_function_buffers.append((name, buffer))
        open_namespace(callback=buffer.append, namespace=namespace)
        for identifier in identifiers:
            if (identifier in separate_function_buffer_by_function_name):
                raise RuntimeError(
                    "separate_files_main_namespace:"
                    " ambiguous assignment: %s" % identifier)
            separate_function_buffer_by_function_name[identifier] = buffer
    #
    separate_namespaces = {}
    separate_namespaces_buffers = {}
    for name, identifiers in separate_files_separate_namespace.items():
        if (len(identifiers) == 0):
            raise RuntimeError(
                "separate_files_separate_namespace: empty list: %s" % name)
        buffers = hpp_cpp_buffers()
        for ext in ["hpp", "cpp"]:
            buffer = getattr(buffers, ext)
            if (ext == "hpp"):
                include_guard(callback=buffer.append,
                              namespace=name, suffix="_HPP")
                buffer.append(include_fem_hpp)
            else:
                buffer.append('#include "%s.hpp"' % name)
            buffer.append("")
            open_namespace(callback=buffer.append, namespace=name)
        for identifier in identifiers:
            if (identifier in separate_namespaces):
                raise RuntimeError(
                    "separate_files_separate_namespace:"
                    " ambiguous assignment: %s" % identifier)
            separate_namespaces[identifier] = name
            separate_namespaces_buffers[identifier] = buffers
    #
    if (not need_function_hpp):
        function_declarations = None
        function_definitions = None
    else:
        function_declarations = []
        function_definitions = []
    #
    global_conv_info = global_conversion_info(
        topological_fprocs=topological_fprocs,
        dynamic_parameters=dynamic_parameters,
        fortran_file_comments=fortran_file_comments,
        fem_do_safe=fem_do_safe,
        arr_nd_size_max=arr_nd_size_max,
        inline_all=inline_all,
        converted_commons_info=converted_commons_info,
        separate_namespaces=separate_namespaces,
        data_values_block_size=data_values_block_size,
        data_specializations=data_specializations)
    #
    for fproc in topological_fprocs.bottom_up_list:
        if (fproc.is_program()):
            continue
        if (fproc.name.value in suppress_functions):
            continue
        hpp_callback = None
        cpp_callback = None
        suppress_cpp = (fproc.name.value in suppress_function_definitions)
        buffers = separate_namespaces_buffers.get(fproc.name.value)
        if (buffers is None):
            if (not need_function_hpp):
                if (not suppress_cpp):
                    cpp_callback = callback
            else:
                function_hpp_buffer = []
                function_declarations.append(function_hpp_buffer)
                hpp_callback = function_hpp_buffer.append
                if (not suppress_cpp):
                    buffer = separate_function_buffer_by_function_name.get(
                        fproc.name.value)
                    if (buffer is None):
                        if (number_of_function_files is None):
                            cpp_callback = callback
                        else:
                            function_cpp_buffer = []
                            function_definitions.append(function_cpp_buffer)
                            cpp_callback = function_cpp_buffer.append
                    else:
                        cpp_callback = buffer.append
        else:
            hpp_callback = buffers.hpp.append
            if (not suppress_cpp):
                cpp_callback = buffers.cpp.append
        if (cpp_callback is None):
            cpp_diverted = []
            cpp_callback = cpp_diverted.append
            if (hpp_callback is None):
                hpp_callback = callback
        if (not need_function_hpp):
            fwds = topological_fprocs.forward_uses_by_identifier.get(
                fproc.name.value)
            if (fwds is not None):
                for fwd_identifier in fwds:
                    fwd_fproc = all_fprocs.fprocs_by_name()[fwd_identifier]
                    try:
                        convert_to_cpp_function(
                            hpp_callback=None,
                            cpp_callback=cpp_callback,
                            conv_info=global_conv_info.specialized(
                                fproc=fwd_fproc),
                            declaration_only=True)
                    except Exception:
                        if (not debug):
                            raise
                        show_traceback()
        try:
            convert_to_cpp_function(
                hpp_callback=hpp_callback,
                cpp_callback=cpp_callback,
                conv_info=global_conv_info.specialized(fproc=fproc),
                force_not_implemented=(fproc.name.value in force_not_implemented))
        except Exception:
            if (not debug):
                raise
            show_traceback()
    #
    for name, buffer in separate_function_buffers:
        close_namespace(
            callback=buffer.append, namespace=namespace, hpp_guard=False)
        if (write_separate_files_main_namespace == "All"
                or name in write_separate_files_main_namespace):
            with open(name+".cpp", "w") as f:
                print("\n".join(break_lines(cpp_text=buffer)), file=f)
    #
    for name, identifiers in separate_files_separate_namespace.items():
        buffers = separate_namespaces_buffers[identifiers[0]]
        for ext in ["hpp", "cpp"]:
            buffer = getattr(buffers, ext)
            close_namespace(
                callback=buffer.append, namespace=name, hpp_guard=(ext == "hpp"))
            if (write_separate_files_separate_namespace == "All"
                    or name in write_separate_files_separate_namespace):
                with open(name+"."+ext, "w") as f:
                    print("\n".join(break_lines(cpp_text=buffer)), file=f)
    #
    if (function_declarations is not None):
        def write_functions(buffers, serial=None):
            if (buffers is function_declarations):
                assert serial is None
                fn = "functions.hpp"
            elif (serial is None):
                fn = "functions.cpp"
            else:
                fn = "functions_%03d.cpp" % serial
            f = open(fn, "w")
            def fcb(line): print(line, file=f)
            if (buffers is function_declarations):
                include_guard(
                    callback=fcb, namespace=namespace, suffix="_FUNCTIONS_HPP")
                fcb(include_with_prefix("cmn"))
                include_separate(callback=fcb)
            else:
                fcb(include_with_prefix("functions"))
            fcb("")
            open_namespace(
                callback=fcb, namespace=namespace, using_namespace_major_types=False)
            for lines in buffers:
                for line in break_lines(cpp_text=lines):
                    fcb(line)
            close_namespace(
                callback=fcb,
                namespace=namespace,
                hpp_guard=(buffers is function_declarations))
            f.close()
        write_functions(function_declarations)
        if (function_definitions is not None and len(function_definitions) != 0):
            buffer_blocks = create_buffer_blocks(
                target_number_of_blocks=number_of_function_files,
                buffers=function_definitions)
            if (len(buffer_blocks) == 1):
                write_functions(buffers=buffer_blocks[0])
            else:
                serial = 0
                for buffers in buffer_blocks:
                    serial += 1
                    write_functions(buffers=buffers, serial=serial)
    #
    hpp_guard = (include_guard_suffix is not None)
    if (suppress_program):
        close_namespace(
            callback=callback, namespace=namespace, hpp_guard=hpp_guard)
    else:
        try:
            convert_program(
                callback=callback,
                global_conv_info=global_conv_info,
                namespace=namespace,
                hpp_guard=hpp_guard,
                debug=debug)
        except Exception:
            if (not debug):
                raise
            show_traceback()

    # Rewrite single-character string literals for CHARACTER*1 variables (char scalars).
    # In view mode (FABLE_SMALL_CHAR=0), scalar CHARACTER dummies/locals are fem::str, so
    # rewriting "A" -> 'A' would break assignments like: fem::str_ref dist; dist = "S";
    if not FABLE_SMALL_CHAR_VIEW:
        result = [rewrite_single_char_string_literals(line) for line in result]

    result = _postprocess_mplapack_labels_and_comments(result)
    result = _normalize_fortran_comment_prefix(result)
    result = _postprocess_complex_initializers(result)
    result = _postprocess_complex_constant_assignments(result)

    # Final intrinsic cleanup (abs / pow2 aliases).
    result = _postprocess_intrinsic_aliases(result)

    # Uppercase selected math intrinsics (ATAN2, COS, SIN, TAN, LOG, EXP, MAX, MIN, ABS).
    result = _postprocess_math_intrinsics_upper(result)

    # Drop redundant parentheses around MIN/MAX index shifts.
    result = _postprocess_minmax_parens(result)

    # Apply MPLAPACK name mapping inside comment lines.
    result = _postprocess_comment_name_map(result)

    # Apply MPLAPACK name mapping inside iMlaenv calls.
    result = _postprocess_ilaenv_name_map(result)

    # Convert character concatenation in iMlaenv to CHAR2/CHAR3.
    result = _postprocess_ilaenv_char_concat(result)

    # Strip C-style float suffixes from literals (1.0f -> 1.0, etc.).
    result = _postprocess_strip_float_suffix(result)

    # Strip leftover Fortran kind suffixes on literals (e.g. 0.0f_wp -> 0.0).
    result = _postprocess_strip_wp_kind_suffix(result)

    #
    result = _postprocess_index_zero_simplify(result)

    # Rewrite Mmaxloc(array(start, end), dim) to Mmaxloc(array, start, end, dim)
    result = _postprocess_mmaxloc(result)

    # Convert array slice assignments to for loops
    result = _postprocess_slice_assignment(result)

    # Widen fem::str<N> scalar locals when substring slices exceed N.
    # This fixes cases like:
    #   CHARACTER(LEN=3) PATH
    # where the upstream parser may drop LEN=3 and we would otherwise emit
    # fem::str<1> path but still generate path(2, 3) = ...
    result = _postprocess_fix_scalar_str_length_from_slices(result)
    
    # Clean up temporary Fortran files created for preprocessing.
    for tmp in temp_files:
        try:
            os.remove(tmp)
        except OSError:
            pass

    if (top_cpp_file_name is not None):
        with open(top_cpp_file_name, "w") as f:
            print("\n".join(result), file=f)

    return result
