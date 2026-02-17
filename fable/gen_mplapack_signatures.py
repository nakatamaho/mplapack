#!/usr/bin/env python3
import sys
import re
from typing import List, Dict, Tuple

# ---------------------------------------------------------------------------
# Prototype parsing helpers
# ---------------------------------------------------------------------------


def collect_prototypes(text: str) -> List[str]:
    """
    Collect C++ prototypes terminated by ';' from the given header text.
    Handles prototypes that span multiple lines.
    """
    protos: List[str] = []
    buf = ""
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        buf += " " + stripped
        if ";" in stripped:
            parts = buf.split(";")
            for p in parts[:-1]:
                stmt = p.strip()
                if stmt:
                    protos.append(stmt + ";")
            buf = parts[-1]
    return protos


def split_args(arg_str: str) -> List[str]:
    """
    Split an argument list string into individual arguments,
    respecting parentheses (for function pointers).
    """
    args: List[str] = []
    cur: List[str] = []
    depth = 0
    for ch in arg_str:
        if ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            arg = "".join(cur).strip()
            if arg:
                args.append(arg)
            cur = []
        else:
            cur.append(ch)
    if cur:
        arg = "".join(cur).strip()
        if arg:
            args.append(arg)
    return args


def classify_return_type(rtype: str) -> str:
    """
    Classify the return type into one of:
      - "void"
      - "real"
      - "complex"
      - "integer"
      - "other"

    This is intentionally coarse and tailored to mpblas/mplapack style
    prototypes (REAL, COMPLEX, INTEGER, void, etc.).
    """
    s = rtype.strip()
    if not s:
        return "void"

    # Remove common qualifiers / reference / pointer markers
    base = s.replace("const", " ")
    base = base.replace("&", " ")
    base = base.replace("*", " ")
    base = re.sub(r"\s+", " ", base).strip()
    lower = base.lower()
    upper = base.upper()

    # void (including things like "static inline void")
    if re.search(r"\bvoid\b", lower):
        return "void"

    # COMPLEX* types
    if "COMPLEX" in upper:
        return "complex"

    # REAL / DOUBLE PRECISION / double / float → treat as "real"
    if "REAL" in upper or "double" in lower or "float" in lower:
        return "real"

    # INTEGER / int
    if "INTEGER" in upper or re.search(r"\bint\b", lower):
        return "integer"

    return "other"


def classify_arg(arg: str) -> str:
    """
    Classify a single argument into one of:
      - "VAL":        value (integer/real/logical by value)
      - "REF_SCALAR": reference to scalar (REAL&, COMPLEX&, INTEGER&)
      - "REF_ARRAY4": fixed-length numeric array reference (ISEED(4))
      - "PTR_CHAR_IN":  pointer to const character data (const char*, char const*)
      - "PTR_CHAR_OUT": pointer to mutable character data (char*)
      - "PTR_NUMERIC":pointer to numeric data (REAL*, COMPLEX*, INTEGER*)
      - "PTR_OTHER":  other pointer types (function pointers, etc.)
    """
    arg = arg.strip()
    if not arg or arg == "void":
        return None

    # Function pointer: we never try to add '&' to these.
    if "(*" in arg:
        return "PTR_OTHER"

    # Fixed-length numeric array reference (ISEED(4)):
    #   INTEGER (&iseed)[4] / mplapackint (&iseed)[4]
    # Only classify this specific argument to avoid affecting other arrays.
    m = re.search(r"\(\s*&\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)\s*\[\s*(\d+)\s*\]", arg)
    if m:
        name = m.group(1)
        try:
            n = int(m.group(2))
        except Exception:
            n = None
        if name.lower() == "iseed" and n == 4:
            return "REF_ARRAY4"

    # Pointer types (check before reference since we look for '*')
    if "*" in arg:
        # Character pointers: distinguish constness of the pointee.
        #
        # We detect "const char*" and "char const*" as input-only.
        # Note: "char * const" is a const pointer to mutable data; treat it as OUT.
        if re.search(r"\bchar\b", arg):
            low = arg.lower()
            if (re.search(r"\bconst\s+char\b", low)
                    or re.search(r"\bchar\s+const\b", low)):
                return "PTR_CHAR_IN"
            return "PTR_CHAR_OUT"
        # Numeric pointers (REAL / COMPLEX / INTEGER)
        base = arg.replace("const", " ").replace("&", " ")
        if "REAL" in base or "COMPLEX" in base or "INTEGER" in base:
            return "PTR_NUMERIC"

        return "PTR_OTHER"

    # Reference types (REAL&, COMPLEX&, INTEGER&)
    # These are scalar references passed as name[0] in the converted code.
    if "&" in arg:
        # Check if it's a numeric reference
        base = arg.replace("const", " ").replace("&", " ")
        if "REAL" in base or "COMPLEX" in base or "INTEGER" in base:
            return "REF_SCALAR"
        # Other reference types (rare) - treat as value
        return "VAL"

    # Everything else is treated as a value argument.
    return "VAL"


def process_header(text: str) -> Dict[str, Tuple[str, List[str]]]:
    """
    Process a header text and return a mapping:
        lowercased_function_name -> (return_kind, [arg_kind, ...])
    """
    sigs: Dict[str, Tuple[str, List[str]]] = {}
    for proto in collect_prototypes(text):
        # Match: return_type name(arglist);
        m = re.match(r"\s*(.+?)\s+([A-Za-z0-9_]+)\s*\((.*)\)\s*;", proto)
        if not m:
            # Skip lines that do not look like prototypes
            continue
        rtype, name, arglist = m.groups()
        rkind = classify_return_type(rtype)
        args = split_args(arglist)
        kinds: List[str] = []
        for a in args:
            k = classify_arg(a)
            if k is not None:
                kinds.append(k)
        sigs[name.lower()] = (rkind, kinds)
    return sigs

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main(argv: List[str]) -> None:
    if len(argv) < 2:
        print(
            "Usage: gen_mplapack_signatures.py mpblas_generic.h mplapack_generic.h [...]",
            file=sys.stderr,
        )
        sys.exit(1)

    all_returns: Dict[str, str] = {}
    all_sigs: Dict[str, List[str]] = {}

    for path in argv[1:]:
        with open(path, "r") as f:
            text = f.read()
        sigs = process_header(text)
        # Later headers override earlier ones if names collide.
        for name, (rkind, kinds) in sigs.items():
            all_returns[name] = rkind
            all_sigs[name] = kinds

    # Emit a Python module that can be imported from cout.py
    print("# Auto-generated by gen_mplapack_signatures.py")
    print("# Do not edit by hand.")
    print()

    # Return types (coarse classification)
    print("FUNCTION_RETURNS = {")
    for name in sorted(all_returns.keys()):
        rkind = all_returns[name]
        print(f"    {name!r}: {rkind!r},")
    print("}")
    print()

    # Argument kinds (existing behavior)
    print("FUNCTION_SIGNATURES = {")
    for name in sorted(all_sigs.keys()):
        kinds = all_sigs[name]
        kinds_list = ", ".join(repr(k) for k in kinds)
        print(f"    {name!r}: [{kinds_list}],")
    print("}")
    print()


if __name__ == "__main__":
    main(sys.argv)
