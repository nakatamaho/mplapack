. "$(dirname "$0")/clang_format_common.sh"

for f in *.cpp *.cc *.h *.hpp *.h.in; do
    run_clang_format "$f"
done
