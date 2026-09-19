#!/bin/sh
set -eu

repo="${1:?repository or bundle path is required}"
ref="${2:?ref is required}"
destination="${3:?destination is required}"

if [ -f "$repo" ]; then
    git clone --no-checkout "$repo" "$destination"
    git -C "$destination" checkout "$ref"
elif printf '%s\n' "$ref" | grep -Eq '^[0-9a-fA-F]{7,40}$'; then
    git init "$destination"
    git -C "$destination" remote add origin "$repo"
    git -C "$destination" fetch --depth=1 origin "$ref"
    git -C "$destination" checkout --detach FETCH_HEAD
else
    git clone --depth 1 --branch "$ref" "$repo" "$destination"
fi
