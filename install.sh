#!/bin/bash
set -eu
cd "$(dirname "$0")"

build_jobs="${JOBS:-4}"
mkdir -p bin
sh script/prepare_dependencies.sh --bcalm

# The pinned BCALM/GATB CMake files predate CMake 4's policy minimum.
cmake -S bcalm2 -B bcalm2/build -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build bcalm2/build --target bcalm --parallel "$build_jobs"
cp bcalm2/build/bcalm bin/bcalm
make -j"$build_jobs"
