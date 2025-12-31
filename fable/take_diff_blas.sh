ref_dir="$HOME/mplapack/mpblas/reference"

comm -12 \
  <(ls *.cpp | sort) \
  <(cd "$ref_dir" && ls *.cpp | sort) \
| while read -r cpp; do
    # Here, $cpp exists in both BLAS/SRC and reference.
    # Example actions:

    # 1) Show diff:
    diff -u "$ref_dir/$cpp" "$cpp" || true

    # 2) Copy new BLAS/SRC version over the reference version:
    # cp "$cpp" "$ref_dir/$cpp"

    echo "common: $cpp"
done
