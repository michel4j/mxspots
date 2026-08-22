# 0004 - scikit-build-core and CMake for Shared Library Compilation

`mxspots` uses `scikit-build-core` as its build backend in `pyproject.toml` to orchestrate CMake compilation of `libmxspots`. This enables standard cross-platform C compilation (`.so`/`.dylib`/`.dll`), easy inclusion of compiler optimization flags (Release builds), and seamless integration with Python wheel building.
