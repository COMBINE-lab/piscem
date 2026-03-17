#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "invoking cargo clean"
cargo clean --target-dir ${SCRIPT_DIR}/target
