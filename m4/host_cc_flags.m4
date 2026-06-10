dnl ===================================================================
dnl M4 Macro: Detect stable C compiler ID/version tag for directory names
dnl
dnl Goals
dnl   * Strip wrapper commands from $CC (ccache/sccache/distcc)
dnl   * Normalize compiler family to a 1-word ID: gcc | clang | icx | icc | intel | cc
dnl   * Extract version as major_minor_patch (fill missing parts with 'x')
dnl   * Export Makefile-substitutable variables:
dnl       MPLAPACK_CC_ID
dnl       MPLAPACK_CC_VER
dnl       MPLAPACK_CC_TAG   (= <id>-<ver>)
dnl
dnl Portability notes
dnl   * Use POSIX-ish shell and utilities.
dnl   * Avoid sed/grep regex constructs that are known to vary widely.
dnl   * Inside AC_DEFUN bodies, M4 strips one bracket level.  Therefore,
dnl     regex bracket expressions must be written doubled (e.g. [[0-9]]).
dnl ===================================================================
AC_DEFUN([AX_MPLAPACK_CC_TAG],[
  AC_REQUIRE([AC_PROG_CC])

dnl -----------------------------------------------------------------
dnl 1. Strip wrapper commands from $CC to find the real compiler binary
dnl -----------------------------------------------------------------
  AC_MSG_CHECKING([C compiler id (ccache-stripped)])

  mplapack_cc_cmd="$CC"
  mplapack_cc_real=""

dnl Tokenize $CC on whitespace.  This targets the common cases:
dnl   CC="ccache gcc"
dnl   CC="sccache clang"
dnl   CC="distcc gcc"
dnl   CC="icx -qopenmp"
  set -- $mplapack_cc_cmd
  while test $# -gt 0; do
    mplapack_cc_tok="$1"
    mplapack_cc_tok_base=$(basename "$mplapack_cc_tok")
    case "$mplapack_cc_tok_base" in
      ccache|sccache|distcc)
        shift
        ;;
      *)
        mplapack_cc_real="$mplapack_cc_tok"
        break
        ;;
    esac
  done

  if test -z "$mplapack_cc_real"; then
    mplapack_cc_real="$CC"
  fi

  mplapack_cc_base=$(basename "$mplapack_cc_real")

dnl Force C locale where possible to reduce localization variance.
dnl Capture both stdout and stderr because some compilers print version to stderr.
  mplapack_cc_version_out=$(LC_ALL=C "$mplapack_cc_real" --version 2>&1 | head -n 5)
  if test -z "$mplapack_cc_version_out"; then
dnl Fallback: run through the full $CC command (may include wrappers/flags).
    mplapack_cc_version_out=$(LC_ALL=C $CC --version 2>&1 | head -n 5)
  fi

dnl -----------------------------------------------------------------
dnl 2. Identify compiler family (normalized ID)
dnl -----------------------------------------------------------------
  mplapack_cc_id="cc"

dnl Prefer basename heuristics first (gcc-13, x86_64-*-gcc, clang-17, etc.)
  case "$mplapack_cc_base" in
    *icx*|*icx.exe)   mplapack_cc_id="icx" ;;
    *icc*|*icc.exe)   mplapack_cc_id="icc" ;;
    *clang*|*clang-cl*|*clang.exe) mplapack_cc_id="clang" ;;
    *gcc*|*gnu*)      mplapack_cc_id="gcc" ;;
  esac

dnl Fall back to --version text matching (handles "cc" symlinks, Apple clang, etc.)
  if test "x$mplapack_cc_id" = "xcc"; then
    if printf '%s\n' "$mplapack_cc_version_out" | grep -Eiq 'clang'; then
      mplapack_cc_id="clang"
    elif printf '%s\n' "$mplapack_cc_version_out" | grep -Eiq 'gcc|GNU'; then
      mplapack_cc_id="gcc"
    elif printf '%s\n' "$mplapack_cc_version_out" | grep -Eiq 'Intel|oneAPI|ICC'; then
dnl Try to separate icx vs icc if the text hints it.
      if printf '%s\n' "$mplapack_cc_version_out" | grep -Eiq 'icx|oneAPI'; then
        mplapack_cc_id="icx"
      elif printf '%s\n' "$mplapack_cc_version_out" | grep -Eiq 'icc'; then
        mplapack_cc_id="icc"
      else
        mplapack_cc_id="intel"
      fi
    fi
  fi

mplapack_cc_version_out=$(LC_ALL=C $CC --version 2>&1 | head -n 5)

dnl Keep the ID filename-safe: [A-Za-z0-9_]
  MPLAPACK_CC_ID=$(printf '%s' "$mplapack_cc_id" | sed 's/[[^A-Za-z0-9_]]/_/g')

  AC_MSG_RESULT([$MPLAPACK_CC_ID])

dnl -----------------------------------------------------------------
dnl 3. Extract version: pick the first valid version from --version output
dnl -----------------------------------------------------------------
  AC_MSG_CHECKING([C compiler version])

  dnl Capture version output.
  mplapack_cc_ver_input=$(LC_ALL=C $CC --version 2>&1 | head -n 5)

  mplapack_cc_ver_triplet=""
  mplapack_cc_major="x"
  mplapack_cc_minor="x"
  mplapack_cc_patch="x"

  if test -n "$mplapack_cc_ver_input"; then
    dnl 1. Remove the compiler command name (e.g., x86_64-...)
    dnl 2. Remove everything inside parentheses (e.g., "(GCC)")
    dnl 3. Convert non-digit/non-dot to spaces
    mplapack_cc_ver_triplet=$(printf '%s\n' "$mplapack_cc_ver_input" \
      | sed "s|^[[^ ]]*||" \
      | sed 's/([[^)]]*)//g' \
      | sed 's/[[^0-9.]]/ /g' \
      | awk '{
          for (i = 1; i <= NF; i++) {
            if ($i == "") continue;
            # Case A: Dotted version (11.4.0, 2025.3.2)
            if ($i ~ /^[[0-9]]+[.][[0-9]]+/) {
                if ($i ~ /^[[0-9]]+[.][[0-9]]+[.][[0-9]]+[.][[0-9]]+$/) {
                  split($i, a, "."); print a[[1]] "." a[[2]] "." a[[3]]; exit
                }
                print $i; exit
            }
            # Case B: Single number (10 from 10-win32).
            if ($i ~ /^[[0-9]]+$/) {
                if (length($i) < 4 && $i < 1000) {
                  print $i; exit
                }
            }
          }
        }')
  fi

  dnl Split identified triplet into major, minor, patch
  if test -n "$mplapack_cc_ver_triplet"; then
    if echo "$mplapack_cc_ver_triplet" | grep "\." > /dev/null; then
      dnl If dots exist, use cut
      mplapack_cc_major=$(echo "$mplapack_cc_ver_triplet" | cut -d. -f1)
      mplapack_cc_minor=$(echo "$mplapack_cc_ver_triplet" | cut -d. -f2)
      mplapack_cc_patch=$(echo "$mplapack_cc_ver_triplet" | cut -d. -f3)
    else
      dnl No dots (e.g., "10"): only major is set
      mplapack_cc_major="$mplapack_cc_ver_triplet"
      mplapack_cc_minor=""
      mplapack_cc_patch=""
    fi
  fi

  dnl Final fallback to 'x' if empty
  test -z "$mplapack_cc_major" && mplapack_cc_major="x"
  test -z "$mplapack_cc_minor" && mplapack_cc_minor="x"
  test -z "$mplapack_cc_patch" && mplapack_cc_patch="x"

  MPLAPACK_CC_VER="${mplapack_cc_major}_${mplapack_cc_minor}_${mplapack_cc_patch}"
  MPLAPACK_CC_VER=$(printf '%s' "$MPLAPACK_CC_VER" | sed 's/[[^A-Za-z0-9_]]/_/g')

  AC_MSG_RESULT([$MPLAPACK_CC_VER])

dnl -----------------------------------------------------------------
dnl 4. Assemble tag (<id>-<ver>) and export to Makefiles
dnl -----------------------------------------------------------------
  AC_MSG_CHECKING([C compiler tag])

  MPLAPACK_CC_TAG="${MPLAPACK_CC_ID}-${MPLAPACK_CC_VER}"
  MPLAPACK_CC_TAG=$(printf '%s' "$MPLAPACK_CC_TAG" | sed 's/[[^A-Za-z0-9_-]]/_/g')

  AC_MSG_RESULT([$MPLAPACK_CC_TAG])

  AC_SUBST([MPLAPACK_CC_ID])
  AC_SUBST([MPLAPACK_CC_VER])
  AC_SUBST([MPLAPACK_CC_TAG])
])
