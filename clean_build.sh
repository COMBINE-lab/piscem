#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "cleaning external dependency directories of cuttlefish and piscem-cpp."
rm -fr ${SCRIPT_DIR}/cuttlefish/external/*
rm -fr ${SCRIPT_DIR}/piscem-cpp/external/zlib-cloudflare

echo "invoking cargo clean"
cargo clean --target-dir ${SCRIPT_DIR}/target
