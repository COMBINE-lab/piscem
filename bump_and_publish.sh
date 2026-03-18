#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--dry-run]
  ./bump_and_publish.sh [--dry-run] <version>

Options:
  --dry-run  Show what would be done without modifying tracked files, creating commits or tags, pushing, or publishing
  -h, --help Show this help message
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            if [[ -n "$VERSION" ]]; then
                die "version specified more than once"
            fi
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || {
    usage
    exit 1
}

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
TAG="v${VERSION}"
CRATE_NAME="piscem"
TMP_TARGET_DIR=""
MANIFEST_BACKUP=""
LOCKFILE_BACKUP=""
MANIFEST_UPDATED=false
COMMIT_CREATED=false

cleanup() {
    local status=$?

    if [[ -n "$TMP_TARGET_DIR" && -d "$TMP_TARGET_DIR" ]]; then
        rm -rf "$TMP_TARGET_DIR"
    fi

    if [[ "$status" -ne 0 && "$DRY_RUN" == false && "$MANIFEST_UPDATED" == true && "$COMMIT_CREATED" == false ]]; then
        if [[ -n "$MANIFEST_BACKUP" && -f "$MANIFEST_BACKUP" ]]; then
            cp "$MANIFEST_BACKUP" "$ROOT_CARGO"
        fi
        if [[ -n "$LOCKFILE_BACKUP" && -f "$LOCKFILE_BACKUP" ]]; then
            cp "$LOCKFILE_BACKUP" "$LOCKFILE"
        fi
        echo "restored $ROOT_CARGO and $LOCKFILE after failure" >&2
    fi

    if [[ -n "$MANIFEST_BACKUP" && -f "$MANIFEST_BACKUP" ]]; then
        rm -f "$MANIFEST_BACKUP"
    fi
    if [[ -n "$LOCKFILE_BACKUP" && -f "$LOCKFILE_BACKUP" ]]; then
        rm -f "$LOCKFILE_BACKUP"
    fi

    return "$status"
}

trap cleanup EXIT

[[ -f "$ROOT_CARGO" ]] || die "not found: $ROOT_CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE"

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine current crate version from $ROOT_CARGO"

if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    die "crate version is already $VERSION"
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi

if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi

if ! git remote get-url origin >/dev/null 2>&1; then
    die "git remote 'origin' is not configured"
fi

echo "Current crate version : $CURRENT_VERSION"
echo "New crate version     : $VERSION"
echo "Tag                   : $TAG"
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run               : yes"
else
    echo "Dry-run               : no"
fi
echo

echo "Preflight checks before changing version"
cargo check -q
TMP_TARGET_DIR="$(mktemp -d "${TMPDIR:-/tmp}/piscem-release-check.XXXXXX")"
CARGO_TARGET_DIR="$TMP_TARGET_DIR" cargo package --offline --allow-dirty --no-verify >/dev/null
rm -rf "$TMP_TARGET_DIR"
TMP_TARGET_DIR=""

echo "Updating $ROOT_CARGO"
echo "  version: $CURRENT_VERSION -> $VERSION"
echo "Updating $LOCKFILE"
echo "  package entry version: $CURRENT_VERSION -> $VERSION"

if [[ "$DRY_RUN" == false ]]; then
    MANIFEST_BACKUP="$(mktemp "${TMPDIR:-/tmp}/piscem-Cargo.toml.XXXXXX")"
    LOCKFILE_BACKUP="$(mktemp "${TMPDIR:-/tmp}/piscem-Cargo.lock.XXXXXX")"
    cp "$ROOT_CARGO" "$MANIFEST_BACKUP"
    cp "$LOCKFILE" "$LOCKFILE_BACKUP"

    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$ROOT_CARGO"
    rm -f "${ROOT_CARGO}.bak"

    sed -i.bak "/^name = \"${CRATE_NAME}\"$/,/^dependencies = \\[$/s/^version = \".*\"/version = \"${VERSION}\"/" "$LOCKFILE"
    rm -f "${LOCKFILE}.bak"

    MANIFEST_UPDATED=true
else
    echo "Dry-run: would rewrite $ROOT_CARGO and $LOCKFILE"
fi

UPDATED_VERSION="$CURRENT_VERSION"
UPDATED_LOCK_VERSION="$CURRENT_VERSION"
if [[ "$DRY_RUN" == false ]]; then
    UPDATED_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
    UPDATED_LOCK_VERSION="$(sed -n "/^name = \"${CRATE_NAME}\"$/,/^dependencies = \\[$/s/^version = \"\(.*\)\"/\1/p" "$LOCKFILE" | head -1)"

    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "crate version update failed"
    [[ "$UPDATED_LOCK_VERSION" == "$VERSION" ]] || die "lockfile version update failed"
fi

echo
echo "Post-bump validation"
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run: would run cargo check and cargo package against the bumped version"
else
    cargo check -q
    TMP_TARGET_DIR="$(mktemp -d "${TMPDIR:-/tmp}/piscem-release-check.XXXXXX")"
    CARGO_TARGET_DIR="$TMP_TARGET_DIR" cargo package --offline --allow-dirty --no-verify >/dev/null
    rm -rf "$TMP_TARGET_DIR"
    TMP_TARGET_DIR=""
fi

run git add "$ROOT_CARGO" "$LOCKFILE"
run git commit -m "chore(release): bump ${CRATE_NAME} to v${VERSION}"

if [[ "$DRY_RUN" == false ]]; then
    COMMIT_CREATED=true
fi

run cargo publish
run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

echo
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run complete"
else
    echo "Release bump and publish complete for v${VERSION}"
fi
