# 0009 - Azimuthal Radial Profiling for Canonical Ice Ring Detection

To detect and mask hexagonal water ice contamination, `mxspots` computes a 1D azimuthal radial intensity profile binned from the beam center and evaluates signal-to-noise ratio in canonical Ice $I_h$ resolution shells. Active ice rings are detected in sub-10ms and masked during candidate pixel thresholding to prevent false-positive spot picking.
