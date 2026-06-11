#!/bin/sh
# entrypoint.sh - Personality wrapper for 32-bit containers on 64-bit hosts
#
# Problem: Docker --platform linux/386 on x86_64 hosts shares the host kernel,
#          so uname -m returns "x86_64" despite 32-bit userland.  This causes
#          config.guess / GMP / autotools to select 64-bit code paths, which
#          then fail at build time (e.g., "mp_limb_t is 32 bits but assembler
#          expects 64 bits").
#
# Solution: Use linux32 (setarch i686) to set the correct personality so that
#           uname -m reports "i686".
#
# Usage in Dockerfile:
#   COPY common/entrypoint.sh /entrypoint.sh
#   RUN chmod +x /entrypoint.sh
#   ENTRYPOINT ["/entrypoint.sh"]
#   CMD ["/bin/bash", "-lc", "cd /work && ./configure && make -j$(nproc)"]

case "$(uname -m)" in
    x86_64)
        # Detect 32-bit userland: check if the default pointer size is 4 bytes,
        # or fall back to dpkg/rpm architecture detection.
        is_32bit=false
        if command -v dpkg >/dev/null 2>&1; then
            [ "$(dpkg --print-architecture 2>/dev/null)" = "i386" ] && is_32bit=true
        elif command -v rpm >/dev/null 2>&1; then
            arch="$(rpm --eval '%{_arch}' 2>/dev/null)"
            case "$arch" in i?86) is_32bit=true ;; esac
        elif [ -f /etc/apk/arch ]; then
            [ "$(cat /etc/apk/arch)" = "x86" ] && is_32bit=true
        fi

        if $is_32bit; then
            exec linux32 "$@"
        fi
        ;;
esac

exec "$@"
