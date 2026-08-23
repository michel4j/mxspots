# 0017 - Multi-Linux Wheel Distribution and Portable Binary Standards

## Context

`mxspots` contains a C shared library (`libmxspots`) compiled with OpenMP support. For seamless user installations without requiring a local C compiler, OpenMP development headers, or CMake toolchains, binary wheels must be pre-built and distributed for Linux platforms.

Distributing pre-compiled binary wheels across diverse Linux installations introduces two core portability constraints:
1. **Instruction Set Portability**: Host CPU microarchitectures differ across client systems. Hardcoded compiler flags such as `-march=native` emit host-specific vector instructions (AVX-512, AVX2, etc.) that cause fatal illegal instruction (`SIGILL`) faults on older machines.
2. **C Runtime and ABI Compatibility**: Linux distributions use differing C standard libraries (glibc vs musl) and glibc version baselines.

## Decision

1. **Portable Baseline Compilation**:
   - `CMakeLists.txt` defaults to `PORTABLE_BUILD=ON` and `ENABLE_NATIVE_TUNING=OFF`.
   - Wheel builds use generic compiler optimization (`-O3`) targeting baseline instruction sets (`x86-64` and generic `aarch64`).
   - Local builds retain `-DENABLE_NATIVE_TUNING=ON` to allow users to opt into host microarchitecture tuning.

2. **Standard Linux ABI Targets**:
   - **`manylinux_2_28`**: Targets standard glibc-based distributions (Rocky/Alma/RHEL 8+, Ubuntu 20.04+, Debian 11+, Fedora 28+).
   - **`musllinux_1_2`**: Targets Alpine Linux 3.19+ and other musl-based distributions.
   - **Architectures**: `x86_64` (Intel/AMD 64-bit) and `aarch64` (ARM 64-bit).
   - **Python Versions**: CPython 3.12 and 3.13 (`cp312-*`, `cp313-*`).

3. **cibuildwheel and GitHub Actions Integration**:
   - Build and test matrix automated via `cibuildwheel` on GitHub Actions runners with QEMU for ARM64 emulation.
   - Every wheel is verified by running the full test suite (`pytest`) inside the target container environment before wheel packaging and artifact upload.

## Consequences

- End users on modern Linux distributions (both glibc and musl) can install `mxspots` via `pip install mxspots` without compiling from source.
- Binaries execute reliably without CPU microarchitecture mismatch crashes.
