# 0001 - Custom C Engine for Spot Finding and Analysis

For peak detection and frame quality analysis, `mxspots` uses a custom, lightweight C core library/extension rather than relying on pure Python or large external crystallographic software suites. Frames are loaded using `mxio` and passed directly as numerical pixel arrays to the C engine for fast processing.
