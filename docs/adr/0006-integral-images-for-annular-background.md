# 0006 - Integral Images for Annular Local Background Estimation

To compute local annular background mean and variance in $O(1)$ time per pixel regardless of filter radius, `mxspots` constructs 2D integral images (Summed-Area Tables for sum, sum-of-squares, and valid pixel counts). This avoids expensive sliding-window or convolutional passes and accurately handles unmeasured, zero, or masked pixels via the integral count table.
