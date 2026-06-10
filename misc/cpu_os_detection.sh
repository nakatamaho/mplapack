# -----------------------------------------------------------------------------
# Deterministic platform identifier for CI/test output directories.
#
# Output format:
#   <cpu_model>|<os_name><os_version>|<cpu_raw>
#
# Examples:
#   Core_i5-8500B|macos15|Intel(R) Core(TM) i5-8500B CPU @ 3.00GHz
#   Ryzen_Threadripper_3970X|ubuntu2204|AMD Ryzen Threadripper 3970X 32-Core Processor
#   Xeon_Gold_6338|rocky93|Intel(R) Xeon(R) Gold 6338 CPU @ 2.00GHz
#
# Field 1: directory-safe CPU tag    (for OUTDIR)
# Field 2: OS identifier             (for OUTDIR)
# Field 3: human-readable CPU string (for gnuplot MODELNAME)
#
# Design goals:
# - POSIX sh (no bashisms)
# - Works in minimal environments (no non-standard packages required)
# - Never exits without printing *something*
# - Deterministic string normalization for directory naming
# - Field 3 added without breaking existing consumers of fields 1 and 2
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
    tr ' ' '_' |
    sed 's/[^A-Za-z0-9._+-]/_/g; s/__*/_/g; s/^_*//; s/_*$//'
}

detect_cpu_raw() {
    # Allow caller to override CPU identification entirely.
    # Useful inside containers where the virtualization layer hides CPU details.
    # Usage: CPU_MODEL_OVERRIDE="Apple M4" bash cpu_os_detection.sh
    if [ -n "${CPU_MODEL_OVERRIDE:-}" ]; then
        printf '%s\n' "$CPU_MODEL_OVERRIDE"
        return
    fi

    os_s=$(uname -s 2>/dev/null || echo "")

    raw=""

    case "$os_s" in
        Linux*)
            if have_cmd lscpu; then
                raw=$(LC_ALL=C lscpu 2>/dev/null |
                    awk -F: '/Model name/ {sub(/^[ \t]+/, "", $2); print $2; exit}')
            fi

            # Apple Silicon inside Linux containers (Podman/Docker on macOS):
            # lscpu reports Vendor ID "Apple" but Model name "-" or empty.
            # The virtualization layer (Apple VZ/QEMU) hides real CPU part codes,
            # so /proc/cpuinfo CPU part is unreliable (often 0x000).
            # Best-effort: infer generation from ISA extension flags.
            if [ -z "$raw" ] || [ "$raw" = "-" ]; then
                vendor=$(LC_ALL=C lscpu 2>/dev/null |
                    awk -F: '/Vendor ID/ {sub(/^[ \t]+/, "", $2); print $2; exit}')
                if [ "$vendor" = "Apple" ] && [ "$(uname -m 2>/dev/null)" = "aarch64" ]; then
                    flags=$(LC_ALL=C lscpu 2>/dev/null |
                        awk '/^[[:space:]]*Flags:/ || /^[[:space:]]*flags:/ {
                            sub(/^[^:]+:/, ""); print; exit}')
                    has_bf16=0; has_afp=0
                    echo "$flags" | grep -qw "bf16" && has_bf16=1
                    echo "$flags" | grep -qw "afp"  && has_afp=1
                    if   [ "$has_bf16" = "1" ] && [ "$has_afp" = "1" ]; then
                        raw="Apple M3 or M4"   # cannot distinguish inside container
                    elif [ "$has_bf16" = "1" ]; then
                        raw="Apple M2"
                    else
                        raw="Apple M1"
                    fi
                fi
            fi

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
            if have_cmd sysctl; then
                raw=$(sysctl -n machdep.cpu.brand_string 2>/dev/null)
                if [ -z "$raw" ]; then
                    raw=$(sysctl -n hw.model 2>/dev/null)
                fi
            fi
            ;;
        CYGWIN*|MINGW*|MSYS*)
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

    in=$(printf '%s' "$in" | tr -d '\r')
    in=$(printf '%s' "$in" | LC_ALL=C tr -cd '\000-\177' | LC_ALL=C sed 's/(R)//g; s/(TM)//g')
    in=$(printf '%s' "$in" | tr '@' ' ')

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

                if (tl == "intel" || tl == "amd") continue
                if (tl == "cpu" || tl == "processor") continue
                if (t ~ /^[0-9]+-[Cc]ore(s)?$/) continue
                if (t ~ /^[0-9]+$/ && i < n && lc(a[i+1]) ~ /^core(s)?$/) { i++; continue }
                if (t ~ /^[0-9]+(\.[0-9]+)?[GgMm][Hh][Zz]$/) continue
                if (tl == "mhz" || tl == "ghz") continue

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
                . /etc/os-release 2>/dev/null
                os_name=$ID
                os_ver=$VERSION_ID
            fi

            if [ -z "$os_name" ]; then
                os_name="linux"
            fi

            if [ -n "$os_ver" ]; then
                os_ver=$(printf '%s' "$os_ver" | tr '.' '_')
            fi

            os_name=$(printf '%s' "$os_name" | dir_safe_token)

            printf '%s%s\n' "$os_name" "$os_ver"
            ;;
        Darwin*)
            major=""
            if have_cmd sw_vers; then
                major=$(sw_vers -productVersion 2>/dev/null | awk -F. 'NR==1 {print $1; exit}')
            fi

            if [ -z "$major" ]; then
                printf 'macos\n'
            else
                printf 'macos%s\n' "$major"
            fi
            ;;
        CYGWIN*|MINGW*|MSYS*)
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

    if [ -z "$cpu_norm" ]; then cpu_norm="unknown_cpu"; fi
    if [ -z "$os_id" ];    then os_id="unknown_os";    fi
    if [ -z "$cpu_raw" ];  then cpu_raw="$cpu_norm";   fi

    # Field 1: directory-safe CPU tag      (e.g. Ryzen_Threadripper_3970X)
    # Field 2: OS identifier               (e.g. ubuntu2204)
    # Field 3: human-readable CPU string   (e.g. AMD Ryzen Threadripper 3970X 32-Core Processor)
    printf '%s|%s|%s\n' "$cpu_norm" "$os_id" "$cpu_raw"
}

main
