# 0007 - Streaming Two-Pass Run-Length Connected Component Labeling

For spot segmentation, `mxspots` uses a streaming row-by-row run-length encoding (RLE) connected component labeling algorithm with disjoint-set union-find. This bounds temporary memory to compact stack buffers for adjacent rows, prevents recursion stack overflow on high-resolution detector frames, and computes exact spot moments and peak intensity in two streaming passes.
