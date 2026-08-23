# 0011 - 95th Percentile Resolution Limit

To ensure robust evaluation of frame diffraction resolution limits, `mxspots` uses a $d_{95}$ metric computed strictly over confirmed Bragg reflections (spots belonging to regular lattice recurrence components with $\ge 5$ members) rather than all detected spots or the absolute single highest-resolution spot.

This outlier-resistant approach guarantees:
1. Isolated random noise peaks and false positives at extreme detector edges cannot artificially inflate the reported frame resolution.
2. In the absence of regular Bragg reflections ($N_B = 0$), the resolution limit evaluates to unresolvable ($d_{\text{min}} = 999.0\,\text{Å}$).
3. The composite quality score accurately reflects the diffraction power of true crystalline lattice reflections.
