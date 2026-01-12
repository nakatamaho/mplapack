# Standalone Python implementation of FABLE scanner utilities.
# C++ extension `fable_ext` is intentionally ignored here.

# ----------------------------------------------------------------------
# Unsigned integer scan
# ----------------------------------------------------------------------


def py_fem_utils_unsigned_integer_scan(code, start=0, stop=-1):
    """
    Scan digits starting at `start` up to `stop` (exclusive).
    Return index of first non-digit or -1 if no digits were found.
    """
    i = start
    # NOTE: `stop` is assumed to be already normalized to len(code).
    while i < stop:
        c = code[i]
        if not c.isdigit():
            break
        i += 1
    if i == start:
        return -1
    return i


def py_ext_get_code_stop(code, stop):
    """
    Normalize `stop` to a valid end index for `code`.
    If stop < 0, use len(code).
    """
    len_code = len(code)
    if stop < 0:
        return len_code
    assert stop <= len_code
    return stop


def py_unsigned_integer_scan(code, start=0, stop=-1):
    """
    Public wrapper that mimics the C++ unsigned_integer_scan behavior.
    """
    code_stop = py_ext_get_code_stop(code, stop)
    return py_fem_utils_unsigned_integer_scan(
        code=code, start=start, stop=code_stop
    )


# ----------------------------------------------------------------------
# Floating point scans
# ----------------------------------------------------------------------


def py_floating_point_scan_after_exponent_char(code, start=0, stop=-1):
    """
    Scan integer part after an exponent character (e or d).
    Handles optional + / - sign.
    Returns index after the last digit, or -1 if no digits found.
    """
    code_stop = py_ext_get_code_stop(code=code, stop=stop)
    i = start
    if i < code_stop:
        c = code[i]
        if c == "+" or c == "-":
            i += 1
        return py_unsigned_integer_scan(code=code, start=i, stop=stop)
    return -1


def py_floating_point_scan_after_dot(code, start=0, stop=-1):
    """
    Scan fractional part after dot and optional exponent.
    Returns index after the fractional / exponent part.
    """
    code_stop = py_ext_get_code_stop(code=code, stop=stop)
    i = py_unsigned_integer_scan(code=code, start=start, stop=stop)
    if i < 0:
        i = start
    if i < code_stop:
        c = code[i]
        if c == "e" or c == "d":
            return py_floating_point_scan_after_exponent_char(
                code=code, start=i + 1, stop=stop
            )
    return i


# ----------------------------------------------------------------------
# Identifier scan
# ----------------------------------------------------------------------


def py_identifier_scan(code, start=0, stop=-1):
    """
    Scan Fortran-like identifier: [a-z_][a-z0-9_]*
    Returns index of first non-identifier char, or -1 if invalid start.
    """
    code_stop = py_ext_get_code_stop(code=code, stop=stop)
    i = start
    if i < code_stop:
        c = code[i]
        i += 1
        if (c < "a" or c > "z") and c != "_":
            return -1
        while i < code_stop:
            c = code[i]
            i += 1
            if ((c < "a" or c > "z") and (c < "0" or c > "9") and c != "_"):
                return i - 1
        return i
    return -1


# ----------------------------------------------------------------------
# Parenthesis matching
# ----------------------------------------------------------------------


def py_find_closing_parenthesis(code, start=0, stop=-1):
    """
    Find matching ')' starting from index `start`.
    Handles nested parentheses.
    Returns index of the matching ')', or -1 if not found.
    """
    code_stop = py_ext_get_code_stop(code=code, stop=stop)
    n_inner = 0
    for i in range(start, code_stop):
        c = code[i]
        if c == ")":
            if n_inner == 0:
                return i
            n_inner -= 1
        elif c == "(":
            n_inner += 1
    return -1


# ----------------------------------------------------------------------
# Public API bindings (always Python implementation)
# ----------------------------------------------------------------------

unsigned_integer_scan = py_unsigned_integer_scan
floating_point_scan_after_exponent_char = py_floating_point_scan_after_exponent_char
floating_point_scan_after_dot = py_floating_point_scan_after_dot
identifier_scan = py_identifier_scan
find_closing_parenthesis = py_find_closing_parenthesis


class SemanticError(Exception):
    """Semantic error raised by FABLE front-end."""
    pass
