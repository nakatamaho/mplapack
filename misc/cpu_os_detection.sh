#!/bin/sh
# -----------------------------------------------------------------------------
# Deterministic platform identifier for CI/test output directories.
#
# Output format:
#   <cpu_model>-<os_name><os_version>
#
# Examples:
#   Core_i5-8500B-macos15
#   Ryzen_Threadripper_3970X-ubuntu2204
#   Xeon_Gold_6338-rocky93
#
# Design goals:
# - POSIX sh (no bashisms)
# - Works in minimal environments (no non-standard packages required)
# - Never exits without printing *something*
# - Deterministic string normalization for directory naming
# -----------------------------------------------------------------------------

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

# Make a token directory-safe:
# - spaces -> _
# - non [A-Za-z0-9._+-] -> _
# - collapse repeated _
# - trim leading/trailing _
dir_safe_token() {
    # Reads stdin, writes token to stdout.
    tr ' ' '_' |
    sed 's/[^A-Za-z0-9._+-]/_/g; s/__*/_/g; s/^_*//; s/_*$//'
}

detect_cpu_raw() {
    os_s=$(uname -s 2>/dev/null || echo "")

    raw=""

    case "$os_s" in
        Linux*)
            # Priority 1: lscpu
            if have_cmd lscpu; then
                raw=$(LC_ALL=C lscpu 2>/dev/null |
                    awk -F: '/Model name/ {sub(/^[ \t]+/, "", $2); print $2; exit}')
            fi

            # Priority 2: /proc/cpuinfo
            if [ -z "$raw" ] && [ -r /proc/cpuinfo ]; then
                raw=$(awk -F: '
                    {
                        k = tolower($1)
                        if (k ~ /^model name$/ || k ~ /^hardware$/ || k ~ /^processor$/) {
                            sub(/^[ \t]+/, "", $2)
                            print $2
                            exit
                        }
                    }' /proc/cpuinfo 2>/dev/null)
            fi
            ;;
        Darwin*)
            # Priority: sysctl machdep.cpu.brand_string, fallback hw.model
            if have_cmd sysctl; then
                raw=$(sysctl -n machdep.cpu.brand_string 2>/dev/null)
                if [ -z "$raw" ]; then
                    raw=$(sysctl -n hw.model 2>/dev/null)
                fi
            fi
            ;;
        CYGWIN*|MINGW*|MSYS*)
            # Priority: wmic cpu get Name, fallback: powershell.exe
            if have_cmd wmic; then
                raw=$(wmic cpu get Name 2>/dev/null |
                    tr -d '\r' |
                    awk 'NR>1 && $0 !~ /^[[:space:]]*$/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0); print; exit}')
            fi
            if [ -z "$raw" ] && have_cmd powershell.exe; then
                raw=$(powershell.exe -NoProfile -Command \
                    "(Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)" \
                    2>/dev/null |
                    tr -d '\r' |
                    awk 'NR==1 {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0); print; exit}')
            fi
            ;;
    esac

    # Final fallback: uname -m
    if [ -z "$raw" ]; then
        raw=$(uname -m 2>/dev/null || echo "")
    fi

    if [ -z "$raw" ]; then
        raw="unknown_cpu"
    fi

    printf '%s\n' "$raw"
}

normalize_cpu() {
    in=$1

    # Strip CR (common in Windows command output).
    in=$(printf '%s' "$in" | tr -d '\r')

    # Remove trademark markers and symbols.
    # (These are common in brand strings: Intel(R), Core(TM), etc.)
    in=$(printf '%s' "$in" | LC_ALL=C tr -cd '\000-\177' | LC_ALL=C sed 's/(R)//g; s/(TM)//g')

    # Make frequency parsing easier by removing '@'.
    in=$(printf '%s' "$in" | tr '@' ' ')

    # Token-level filtering to remove required noise while keeping family + model.
    filtered=$(printf '%s\n' "$in" | awk '
        function lc(s) { return tolower(s) }
        {
            gsub(/[,;]/, " ")
            n = split($0, a, /[[:space:]]+/)

            out = ""
            for (i = 1; i <= n; i++) {
                t = a[i]
                if (t == "") continue
                tl = lc(t)

                # Remove vendor noise.
                if (tl == "intel" || tl == "amd") continue

                # Remove generic words.
                if (tl == "cpu" || tl == "processor") continue

                # Remove core count phrases: "32-Core", "16-Core", etc.
                if (t ~ /^[0-9]+-[Cc]ore(s)?$/) continue

                # Remove "N Core" (two-token form).
                if (t ~ /^[0-9]+$/ && i < n && lc(a[i+1]) ~ /^core(s)?$/) { i++; continue }

                # Remove frequency info: "3.00GHz", "3000MHz", etc.
                if (t ~ /^[0-9]+(\.[0-9]+)?[GgMm][Hh][Zz]$/) continue
                if (tl == "mhz" || tl == "ghz") continue

                # Keep token.
                if (out == "") out = t
                else out = out " " t
            }

            print out
            exit
        }')

    safe=$(printf '%s' "$filtered" | dir_safe_token)

    if [ -z "$safe" ]; then
        safe="unknown_cpu"
    fi

    printf '%s\n' "$safe"
}

detect_os() {
    os_s=$(uname -s 2>/dev/null || echo "")

    case "$os_s" in
        Linux*)
            os_name=""
            os_ver=""

            if [ -r /etc/os-release ]; then
                # /etc/os-release is designed to be shell-sourceable.
                . /etc/os-release 2>/dev/null
                os_name=$ID
                os_ver=$VERSION_ID
            fi

            # Minimal fallback if /etc/os-release is missing or incomplete.
            if [ -z "$os_name" ]; then
                os_name="linux"
            fi

            # VERSION_ID: remove dots (e.g. 22.04 -> 2204).
            if [ -n "$os_ver" ]; then
                os_ver=$(printf '%s' "$os_ver" | tr '.' '_')
            fi

            # Make os_name directory-safe (rarely needed, but harmless).
            os_name=$(printf '%s' "$os_name" | dir_safe_token)

            printf '%s%s\n' "$os_name" "$os_ver"
            ;;
        Darwin*)
            # macOS: use major version only.
            major=""
            if have_cmd sw_vers; then
                major=$(sw_vers -productVersion 2>/dev/null | awk -F. 'NR==1 {print $1; exit}')
            fi

            # If we cannot read the version, still emit a stable identifier.
            if [ -z "$major" ]; then
                printf 'macos\n'
            else
                printf 'macos%s\n' "$major"
            fi
            ;;
        CYGWIN*|MINGW*|MSYS*)
            # Windows: version is optional by requirement; keep it stable.
            printf 'windows\n'
            ;;
        *)
            printf 'unknown_os\n'
            ;;
    esac
}

main() {
    cpu_raw=$(detect_cpu_raw)
    cpu_norm=$(normalize_cpu "$cpu_raw")
    os_id=$(detect_os)

    # Absolute last-resort guards (must never print empty).
    if [ -z "$cpu_norm" ]; then cpu_norm="unknown_cpu"; fi
    if [ -z "$os_id" ]; then os_id="unknown_os"; fi

    printf '%s|%s\n' "$cpu_norm" "$os_id"
}

main
