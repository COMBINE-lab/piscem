#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <version>  (e.g. $0 0.15.0)" >&2
    exit 1
fi

VERSION="$1"
TAG="v${VERSION}"
CARGO_TOML="Cargo.toml"

# Validate version looks reasonable (digits and dots)
if ! echo "$VERSION" | grep -qE '^[0-9]+\.[0-9]+\.[0-9]+'; then
    echo "Error: version '$VERSION' doesn't look like a semver (expected X.Y.Z)" >&2
    exit 1
fi

# Check we're in the right directory
if [ ! -f "$CARGO_TOML" ]; then
    echo "Error: $CARGO_TOML not found. Run this from the piscem repo root." >&2
    exit 1
fi

# Check for uncommitted changes
if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "Error: working tree has uncommitted changes. Commit or stash first." >&2
    exit 1
fi

# Check tag doesn't already exist
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "Error: tag '$TAG' already exists." >&2
    exit 1
fi

# Update version in Cargo.toml
sed -i '' "s/^version = \".*\"/version = \"${VERSION}\"/" "$CARGO_TOML"

echo "Updated $CARGO_TOML to version $VERSION"

# Commit, tag, and push
git add "$CARGO_TOML"
git commit -m "release ${TAG}"
git tag "$TAG"
git push origin main "$TAG"

echo "Released ${TAG} and pushed to origin."
