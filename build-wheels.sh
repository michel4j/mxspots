#!/usr/bin/env bash
set -e

# Define directories
UNREPAIRED_DIR="pre-dist"
DIST_DIR="dist"

echo "🧹 Cleaning up old builds..."
rm -rf "$UNREPAIRED_DIR" "$DIST_DIR"
mkdir -p "$UNREPAIRED_DIR" "$DIST_DIR"

echo "🐳 Running build inside manylinux Docker container..."
podman run --rm \
  -v "$(pwd)":/io:z \
  -w /io \
  quay.io/pypa/manylinux_2_28_x86_64 \
  bash -c "
    echo '📦 Building raw wheel...'
    /opt/python/cp312-cp312/bin/pip wheel . -w $UNREPAIRED_DIR --no-deps
    /opt/python/cp313-cp313/bin/pip wheel . -w $UNREPAIRED_DIR --no-deps
    /opt/python/cp314-cp314/bin/pip wheel . -w $UNREPAIRED_DIR --no-deps

    echo '🔧 Repairing wheel with auditwheel...'
    auditwheel repair $UNREPAIRED_DIR/*.whl --only-plat -w $DIST_DIR/
  "

echo "🧹 Cleaning up temporary build artifacts..."
rm -rf "$UNREPAIRED_DIR" *.egg-info build/

echo "✅ Success! Your PyPI-compliant wheel is ready in the '$DIST_DIR/' directory."
