#!/usr/bin/env bash
# Fetch the AD backends compared by this benchmark at the versions used for the results in README.md, and print the
# matching CMake arguments. Only sources are downloaded; nothing is built or installed.
# Usage: benchmarks/automatic-differentiation/fetchBackends.sh <destination directory> <directory containing Fastor/>
set -euo pipefail
DEST=$(realpath -m "${1:?destination directory}")
FASTOR_DIR=$(realpath "${2:?directory that contains Fastor/}")
mkdir -p "$DEST"
cd "$DEST"

clone() { # url tag dir
  [ -d "$3" ] || git -c advice.detachedHead=false clone -q --depth 1 --branch "$2" "$1" "$3"
}

# autodiff 1.1.2 with the Eigen 5 patch autodiff#397 (not merged upstream yet), as Marmot's CI uses it
if [ ! -d autodiff ]; then
  clone https://github.com/autodiff/autodiff.git v1.1.2 autodiff
  curl -sfL https://github.com/autodiff/autodiff/commit/a4401bdd338784952f6084b69960cf914ba4e6af.diff | git -C autodiff apply
fi
# Sacado: only packages/sacado of Trilinos is needed
if [ ! -d trilinos ]; then
  git -c advice.detachedHead=false clone -q --depth 1 --branch trilinos-release-17-3-0 --filter=blob:none --sparse \
    https://github.com/trilinos/Trilinos.git trilinos
  git -C trilinos sparse-checkout set packages/sacado
fi
clone https://github.com/ceres-solver/ceres-solver.git 2.2.0 ceres
clone https://github.com/stan-dev/math.git v5.4.0 stan-math
clone https://github.com/rjhogan/Adept-2.git v2.1.3 adept2
if [ ! -d tinyad ]; then # no releases; the commit the results were measured with
  git init -q tinyad
  git -C tinyad fetch -q --depth 1 https://github.com/patr-schm/TinyAD.git 4b48d1a1a588874556a692a3abbdecd0db4c23e1 2>/dev/null ||
    git -C tinyad fetch -q --depth 1 https://github.com/patr-schm/TinyAD.git HEAD
  git -C tinyad -c advice.detachedHead=false checkout -q FETCH_HEAD
fi

# Fastor alone, so that no other autodiff or Boost from the same prefix ends up on the include path
mkdir -p fastor-include
ln -sfn "$FASTOR_DIR/Fastor" fastor-include/Fastor

echo "-DAUTODIFF_INCLUDE_DIR=$DEST/autodiff -DFASTOR_INCLUDE_DIR=$DEST/fastor-include" \
  "-DSACADO_SOURCE_DIR=$DEST/trilinos/packages/sacado -DCERES_SOURCE_DIR=$DEST/ceres" \
  "-DSTAN_SOURCE_DIR=$DEST/stan-math -DTINYAD_SOURCE_DIR=$DEST/tinyad -DADEPT_SOURCE_DIR=$DEST/adept2"
