# 0002 - Standalone Shared Library Loaded via Ctypes

`mxspots` implements its core computational routines in a standalone C shared library (`libmxspots.so`/`.dylib`/`.dll`) rather than a Python C extension module. Python interacts with this C library through `ctypes`, passing contiguous raw array pointers from `mxio` data buffers. This keeps the C codebase simple, free from CPython ABI dependencies, and independently buildable/testable.
