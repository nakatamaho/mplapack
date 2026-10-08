# Shell access to backends.txt, the MPLAPACK backend table.
# Source this file; it defines no variables other than mplapack_backends_file.
#
#   mplapack_backend_names            backend names in table order
#   mplapack_backend_field NAME FIELD FIELD is real, complex, default or traits

mplapack_backends_file="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/backends.txt"

mplapack_backend_names() {
    awk '!/^#/ && NF { print $1 }' "$mplapack_backends_file"
}

mplapack_backend_field() {
    local column
    case "$2" in
        real) column=2 ;;
        complex) column=3 ;;
        default) column=4 ;;
        traits) column=5 ;;
        *) echo "mplapack_backend_field: unknown field '$2'" >&2; return 1 ;;
    esac
    awk -v name="$1" -v column="$column" '
        !/^#/ && NF && $1 == name { print $column; found = 1 }
        END { if (!found) exit 1 }' "$mplapack_backends_file" || {
        echo "mplapack_backend_field: no backend '$1' in $mplapack_backends_file" >&2
        return 1
    }
}
