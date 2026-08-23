# 0011 - 95th Percentile Resolution Limit

To ensure robust evaluation of frame diffraction resolution limits, `mxspots` uses a $d_{95}$ metric (the highest resolution, i.e., smallest $d$-spacing, bounding 95% of detected spots) rather than the absolute maximum resolution spot ($d_{max}$). This outlier-resistant approach prevents isolated random noise peaks at the detector edges from artificially inflating the reported frame resolution and skewing the composite quality score.
