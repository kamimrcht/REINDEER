#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
apply_patch() {
    dependency="$1"
    patch_file="$PWD/patches/$2"
    if git -C "$dependency" apply --reverse --check "$patch_file" 2>/dev/null; then
        return
    fi
    git -C "$dependency" apply --check "$patch_file"
    git -C "$dependency" apply "$patch_file"
}
apply_patch blight blight-gcc15.patch
apply_patch blight blight-parallel-files.patch
if [ "${1:-}" = --bcalm ]; then
    apply_patch bcalm2/gatb-core gatb-gcc15.patch
fi
